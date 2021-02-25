import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()
using LinearAlgebra
using SparseArrays
using ForwardDiff
using OSQP
using RobotDynamics
using RobotZoo: PlanarRocket
using RobotZoo
using StaticArrays
using Plots
include("rocket.jl")


# Planar Rocket model
model = PlanarRocket(max_roll=5.0)
n,m = state_dim(model), control_dim(model)
xeq = [zeros(6); model.m*model.g; 0]
ueq = zeros(2) 
norm(dynamics(model, xeq, ueq)) ≈ 0 # make sure it's an equilibrium point
dt = 0.1  # time step (s)
tf = 25   # time horizon (s)
N = Int(tf / dt) + 1

# Evaluate the continuous and discrete Jacobians
zeq = KnotPoint(xeq,ueq,dt)   # create a `KnotPoint` type that stores everything together
∇f = RobotDynamics.DynamicsJacobian(model)
jacobian!(∇f, model, zeq)
discrete_jacobian!(RK4, ∇f, model, zeq)

# Extract pieces of the Jacobian
A = ∇f.A
B = ∇f.B;

# Cost matrices
Q = Diagonal([
    1.0/5^2; 
    1.0/5^2; 
    1.0/(10*pi/180)^2; 
    10.0/5^2; 
    10.0/5^2; 
    10.0/(10*pi/180)^2;
    1/(.1*model.m*model.g)^2; 
    1/(10*pi/180)^2
])
R = Diagonal([
    1/(20*model.m*model.g)^2,
    1/(deg2rad(10))^2
])

## Visualizer
vis = initialize_visualizer(model)
open(vis)

## Nominal Trajectory
function nominal_trajectory(x0,N,dt)
    #Design a trajectory that linearly interpolates from x0 to the origin
    Xref = [zero(x0) for k = 1:N]
    for k = 1:N
        Xref[k][1:2] .= ((N-k)/(N-1))*x0[1:2]
        Xref[k][7:8] .= x0[7:8]
    end
    for k = 1:N-1
        Xref[k][4:5] .= (Xref[2][1:2]-Xref[1][1:2])/dt
    end
    
    return Xref
end
x0_ref = [-200, 500, 0, 0,0,0, 0,0.] + xeq
Xref = nominal_trajectory(x0_ref, N, dt)
Uref = [copy(ueq) for k = 1:N]
tref = range(0,tf, length=N)

## Infinite-Horizon LQR
struct LQRController
    K::Matrix{Float64}
    Xref::Vector{Vector{Float64}}
    Uref::Vector{Vector{Float64}}
    times::Vector{Float64}
end
get_k(controller, t) = searchsortedlast(controller.times, t)

function get_control(ctrl::LQRController, x, t)
    k = get_k(ctrl, t)
    return ctrl.Uref[k] - ctrl.K*(x - ctrl.Xref[k])
end

function lqr(A,B,Q,R; P=Matrix(Q), tol=1e-8, max_iters=400, verbose=false)
    # initialize the output
    n,m = size(B)
    K = zeros(m,n)
    Kprev = zero(K)
    
    # TODO: implement the Riccati recursion
    for k = 1:max_iters
        K .= (R + B'P*B)\(B'P*A)
        P .= Symmetric(Q + A'P*A - A'P*B*K)
        if norm(K-Kprev) < tol
            verbose && println("Converged in $k iters")
            break
        end
        Kprev .= K
        if k == max_iters
            @warn "Max iterations"
        end
    end
    
    # return the feedback gains and ctg matrices
    return K,P
end

K,Qf = lqr(A,B,Q,R)
ctrl = LQRController(K,Xref,Uref,tref);
Xlqr,Ulqr,tlqr = simulate(model, Xref[1], ctrl, tf=50)
visualize!(vis, model, tlqr[end] / 10, Xlqr)

## MPC Controller 
struct OSQPController
    P::SparseMatrixCSC{Float64,Int}
    q::Vector{Float64}
    A::SparseMatrixCSC{Float64,Int}
    lb::Vector{Float64}
    ub::Vector{Float64}
    Nmpc::Int
    model::OSQP.Model
    Xref::Vector{Vector{Float64}}
    Uref::Vector{Vector{Float64}}
    times::Vector{Float64}
end

function OSQPController(n::Integer, m::Integer, N::Integer, Nref::Integer=N, Nd::Integer=(N-1)*n)
    Np = (N-1)*(n+m)   # number of primals
    P = spzeros(Np,Np)
    q = zeros(Np)
    A = spzeros(Nd,Np)
    lb = zeros(Nd)
    ub = zeros(Nd)
    Xref = [zeros(n) for k = 1:Nref]
    Uref = [zeros(m) for k = 1:Nref]
    tref = zeros(Nref)
    model = OSQP.Model()
    OSQPController(P,q, A,lb,ub, N, model, Xref, Uref, tref)
end

function buildQP!(ctrl::OSQPController, A,B,Q,R,Qf; tol=1e-6, verbose=false)
    Nt = ctrl.Nmpc-1
    Nx = length(ctrl.Xref[1])    # number of states
    Nu = length(ctrl.Uref[1])    # number of controls

    max_roll = deg2rad(model.max_roll)
    xeq = Xref[end]
    aeq = xeq[7:8]  # actuators
    
    H = sparse([kron(Diagonal(I,Nt-1),[R zeros(Nu,Nx); zeros(Nx,Nu) Q]) zeros((Nx+Nu)*(Nt-1), Nx+Nu); zeros(Nx+Nu,(Nx+Nu)*(Nt-1)) [R zeros(Nu,Nx); zeros(Nx,Nu) Qf]])
    b = zeros(Nt*(Nx+Nu))
    C = sparse([
            [B -I zeros(Nx,(Nt-1)*(Nu+Nx))]; 
            zeros(Nx*(Nt-1),Nu) [kron(Diagonal(I,Nt-1), [A B]) zeros((Nt-1)*Nx,Nx)] + [zeros((Nt-1)*Nx,Nx) kron(Diagonal(I,Nt-1),[zeros(Nx,Nu) Diagonal(-I,Nx)])]
    ])
    Z = kron(Diagonal(I,Nt), [0 0 0 1 0 0 0 0 0 0]) #Matrix that picks out all x2 (height)
    Θ = kron(Diagonal(I,Nt), [0 0 0 0 1 0 0 0 0 0]) #Matrix that picks out all x3 (θ)
    U = kron(Diagonal(I,Nt), [zeros(Nu,Nx) I]) #Matrix that picks out all u
    D = [C; Z; Θ; U]
    
    umin = RobotZoo.umin(model)
    umax = RobotZoo.umax(model)
    lb = [zeros(Nx*Nt); zeros(Nt); -max_roll*ones(Nt); kron(ones(Nt),umin-aeq)]
    ub = [zeros(Nx*Nt); Inf*ones(Nt); max_roll*ones(Nt); kron(ones(Nt),umax-aeq)]
    
    Nd = length(ctrl.lb)
    if Nd == Nt*n
        D = C
        ub = zero(ctrl.ub)
        lb = zero(ctrl.lb)
    end
    lb = lb[1:Nd]
    ub = ub[1:Nd]
    
    ctrl.P .= H
    ctrl.A .= D
    ctrl.ub .= ub
    ctrl.lb .= lb
    OSQP.setup!(ctrl.model, P=ctrl.P, q=ctrl.q, A=ctrl.A, l=ctrl.lb, u=ctrl.ub, verbose=verbose, eps_rel=tol, eps_abs=tol, polish=1)
    return nothing
end

function update_QP!(ctrl::OSQPController, x, time)
    t = get_k(ctrl, time)
    
    Nt = ctrl.Nmpc-1     # horizon
    Nx = length(ctrl.Xref[1])    # number of states
    Nu = length(ctrl.Uref[1])    # number of controls
    
    #Update QP problem
    b = ctrl.q
    lb = ctrl.lb
    ub = ctrl.ub
    xref = ctrl.Xref
    xeq = Xref[end]
    N = length(ctrl.Xref)
    for t_h = 1:(Nt-1)
        if (t+t_h) <= N
            b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*(xref[t+t_h] - xeq)
        else
            b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*(xref[end] - xeq)
        end
    end
    if (t+Nt) <= N
        b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*(xref[t+Nt] - xeq)
    else
        b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*(xref[end] - xeq)
    end
    
    lb[1:Nx] .= -A*(x - xeq)
    ub[1:Nx] .= -A*(x - xeq)

    return nothing
end

function get_control(ctrl::OSQPController, x, time)
    k = get_k(ctrl, time)
    update_QP!(ctrl, x, time)
    OSQP.update!(ctrl.model, q=ctrl.q, l=ctrl.lb, u=ctrl.ub)

    #Solve QP
    results = OSQP.solve!(ctrl.model)
    Δu = results.x[1:2]

    return ctrl.Uref[k] + Δu 
end

## Build Unconstrained Controller
n,m = size(model)
Nmpc = 2
mpc1 = OSQPController(n, m, Nmpc, length(Xref))
mpc1.Xref .= Xref
mpc1.Uref .= Uref
mpc1.times .= tref
buildQP!(mpc1, A,B,Q,R,Qf, tol=1e-2)

##
Xmpc1,Umpc1,tmpc1 = simulate(model, Xref[1], mpc1, tf=50)
visualize!(vis, model, tmpc1[end] / 10, Xmpc1)
Xmpc1 - Xlqr

## Add Constraints
Nd = (Nmpc-1)*(n+4)
mpc2 = OSQPController(n, m, Nmpc, length(Xref), Nd)
mpc2.Xref .= Xref
mpc2.Uref .= Uref
mpc2.times .= tref
buildQP!(mpc2, A,B,Q,R,Qf, tol=1e-2, verbose=false)

##
Xmpc2,Umpc2,tmpc2 = simulate(model, Xref[1], mpc2, tf=40)
visualize!(vis, model, tmpc2[end] / 10, Xmpc2)

## Plots
p = plot(layout=(2,2), size=(1000,800))
plot!(p[1], tlqr,  Xlqr,  inds=2:2, xlabel="time (s)", label="LQR", ylabel="altitude (m)")
plot!(p[1], tmpc1, Xmpc1, inds=2:2, label="MPC")
plot!(p[1], tmpc2, Xmpc2, inds=2:2, label="MPC-Con", legend=:topright)

plot!(p[2], tlqr,  [rad2deg.(x) for x in Xlqr], inds=3:3, xlabel="time (s)", label="LQR", ylabel="roll angle (deg)")
plot!(p[2], tmpc1, [rad2deg.(x) for x in Xmpc1], inds=3:3, label="MPC")
plot!(p[2], tmpc2, [rad2deg.(x) for x in Xmpc2], inds=3:3, label="MPC-Con", legend=:bottomleft)
hline!(p[2], [model.max_roll, -model.max_roll], ls=:dash, c=:red, label="limits")

plot!(p[3], tlqr,  Xlqr,  inds=7:7, xlabel="time (s)", label="LQR", ylabel="Thrust (N)")
plot!(p[3], tmpc1, Xmpc1, inds=7:7, label="MPC")
plot!(p[3], tmpc2, Xmpc2, inds=7:7, label="MPC-Con", legend=:none)
hline!(p[3], [RobotZoo.umax(model)[1], RobotZoo.umin(model)[1]], ls=:dash, c=:red, label="limits")

plot!(p[4], tlqr,  [rad2deg.(x) for x in Xlqr],  inds=8:8, xlabel="time (s)", label="LQR", ylabel="Thrust angle (deg)")
plot!(p[4], tmpc1, [rad2deg.(x) for x in Xmpc1], inds=8:8, label="MPC")
plot!(p[4], tmpc2, [rad2deg.(x) for x in Xmpc2], inds=8:8, label="MPC-Con", legend=:none)
hline!(p[4], rad2deg.([RobotZoo.umax(model)[2], RobotZoo.umin(model)[2]]), ls=:dash, c=:red, label="limits")

plot(tlqr[1:end-1], Ulqr, inds=1:1, xlabel="time (s)", label="LQR", ylabel="thrust (N)")
plot!(tmpc1[1:end-1], Umpc1, inds=1:1, label="MPC")
plot!(tmpc2[1:end-1], Umpc2, inds=1:1, label="MPC-Con", legend=:none)
p
