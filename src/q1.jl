import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()
using SparseArrays
include("quadruped.jl")

struct LQRController
    K::Vector{Matrix{Float64}}
    times::Vector{Float64}
    xeq::Vector{Float64}
    ueq::Vector{Float64}
end
function LQRController(K,xeq,ueq,tf)
    LQRController(Matrix.(K), collect(range(0,tf,length=length(K)+1)), xeq, ueq)
end

function get_control!(controller::LQRController, x, t)
    k = get_k(controller, t)
    return controller.ueq - controller.K[k]*(x - controller.xeq)
end
get_k(controller::LQRController, t) = searchsortedlast(controller.times, t)

function visualize(mvis, tf, X)
    fps = Int(round((length(X)-1)/tf))
    anim = MeshCat.Animation(fps)
    for (k,x) in enumerate(X)
        atframe(anim, k) do
            set_configuration!(mvis, x[1:14])
        end
    end
    setanimation!(mvis, anim)
end

function lqr(Q,R,Qf,A,B,N)
    n,m = size(B)
    P = [zeros(n,n) for k = 1:N]
    K = [zeros(m,n) for k = 1:N-1]
    P[end] .= Qf
    for k = reverse(1:N-1) 
        K[k] .= (R + B'P[k+1]*B)\(B'P[k+1]*A)
        P[k] .= Q + A'P[k+1]*A - A'P[k+1]*B*K[k]
    end
    return K,P
end

function kkt_solve(Q,R,Qf,A,B,N,x0)
    n,m = size(B)
    ix,iu = 1:n, 1:m
    np = n*N + (N-1)*m
    nd = n*N
    P = cat(push!(repeat([Q,R],N-1),Qf)..., dims=[1,2])
    D = spzeros(nd,np)
    d = zeros(nd)
    D[ix,ix] .= I(n)
    d[ix] .= -x0
    for k = 1:N-1
        D[ix .+ k*n, (k-1)*(n+m) .+ (1:2n+m)] .= [A B -I(n)]
    end
    H = [P D'; D zeros(nd,nd)]
    g = [zeros(np); d]
    Y = -H\g
    Z = Y[1:np]
    xf = Z[end-n+1:end]
    z = reshape(Z[1:end-n], n+m, :)
    X = [col[ix] for col in eachcol(z)]
    U = [col[iu .+ n] for col in eachcol(z)]
    push!(X,xf)
    return X, U
end

function linear_sim(A,B,K,x0)
    n,m = size(B)
    N = length(K)+1
    X = fill(zeros(n), N)
    U = fill(zeros(m), N-1)
    X[1] .= x0
    for k = 1:N-1
        U[k] = -K[k]*X[k]
        X[k+1] = A*X[k] + B*U[k]
    end
    X,U
end

# Build Model
model = UnitreeA1()
n,m = state_dim(model), control_dim(model)

# Get equilibrium point (HW1)
x0 = initial_state(model)
u0 = zeros(m)
xeq,ueq = newton_solve(model, x0, verbose=false) 

## Visualizer
mvis = initialize_visualizer(model)
open(mvis)

## Continuous Time
using ControlSystems
function lqr_controller!(torques::AbstractVector, t, state::MechanismState)
    q = configuration(state)
    v = velocity(state)
    x = vcat(q,v)
    torques[3:end] .= ueq - K*(x - xeq)
    torques[1:2].= 0 #
    return nothing
end

∇f = jacobian(model, xeq, ueq)
A = Matrix(∇f.A)
B = Matrix(∇f.B)
Q = 130 * Diagonal(ones(size(x0)));
R = 10 * Diagonal(ones(size(u0)));
K = ControlSystems.lqr(A,B,Q,R)
K0 = copy(K)
maximum(real(eigvals(A-B*K)))

##
x_init = copy(xeq)
x_init[18] -= 0.1
x_init[19] -= 0.05
set_configuration!(mvis, xeq[1:14])
set_configuration!(mvis, x_init[1:14])
state = MechanismState(model.mech)
# state = a1.state
set_configuration!(state, x_init[1:14])
set_velocity!(state, x_init[15:end])
# zero_velocity!(state)

t_lqrs, q_lqrs, v_lqrs = simulate(state, 3, lqr_controller!; Δt=1e-3);
norm(v_lqrs[end])
norm(q_lqrs[end]-xeq[1:14])
animation = MeshCat.Animation(mvis, t_lqrs, q_lqrs);
setanimation!(mvis, animation);

## Solve LQR
dt = 0.01
tf = 2.0
times = range(0, tf, step=dt)
norm(dynamics(model, xeq, ueq)) < 1e-10   # test equilibrium
∇f = discrete_jacobian(RK4, model, xeq, ueq, dt)
A = ∇f.A
B = ∇f.B

Q = Diagonal(fill(130,n))
R = Diagonal(fill(10,m))
Qf = Diagonal(fill(130.0,n))
Qf = dare(A,B,Q,R)

set_configuration!(mvis, x_init[1:14])
K,P = lqr(Q,R,Qf,A,B, length(times))
maximum(abs.(eigvals(A+B*K[1])))
X,U = linear_sim(A, B, K, x_init-xeq)
visualize(mvis, [x + xeq for x in X], 100)
norm(X[end][15:end])

## Simulator
function quad_simulate(model, x0, controller; dt=0.05, tf=1.0, mvis=nothing)
    time = range(0, tf, step=dt)
    n,m = state_dim(model), control_dim(model) 
    N = Int(round(tf/dt)) + 1
    X = [@SVector zeros(n) for k = 1:N] 
    U = [@SVector zeros(m) for k = 1:N-1] 
    X[1] = x0

    for k = 1:length(time) - 1
        U[k] = get_control!(controller, X[k], time[k])
        X[k+1] = discrete_dynamics(RK4, model, X[k], U[k], time[k], dt)
        if !isnothing(mvis)
            set_configuration!(mvis, X[k+1][1:14])
            sleep(dt)
        end
    end
    return X,U
end
ctrl = LQRController(K, xeq, ueq, tf)
X,U = quad_simulate(model, x_init, ctrl, dt=1e-2, tf=2.0, mvis=mvis)
visualize(mvis, tf, X)