import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()
using RobotDynamics
using RobotZoo
# using TrajOptPlots
using StaticArrays
using LinearAlgebra
using SparseArrays
using Plots

## Generate Model and equilibrium point
include("quadruped.jl")
model = UnitreeA1() 
x0 = initial_state(model)
xeq,ueq = newton_solve(model, x0)
dt = 0.01
z = KnotPoint(xeq,ueq,dt)
norm(dynamics(model, xeq, ueq))

# Visualization
mvis = initialize_visualizer(model)
open(mvis)
set_configuration!(mvis, xeq[1:15])

## Setup LQR Problem
∇f = RobotDynamics.DynamicsJacobian(model)
discrete_jacobian!(RK4, ∇f, model, z)
A = RobotDynamics.get_static_A(∇f)
B = RobotDynamics.get_static_B(∇f)
Q = Diagonal([fill(1e-2,15); fill(1e-2, 15)]) 
R = Diagonal(@SVector fill(1e-5, 12))
Qf = Diagonal([fill(1e2,15); fill(1e2, 15)]) 
tf = 3.0
N = Int(round(tf/dt)) + 1
x_init = zeros(30) 
# x_init[4] = deg2rad(5)
# x_init[7] = deg2rad(5)    # FR hip
# x_init[8] = -0.02   # FL hip
x_init[9] = 0.3     # RL hip
x_init[10] = 0.1    # FR thigh
# x_init[15] = 0.05   # RL calf
set_configuration!(mvis, xeq[1:15])
set_configuration!(mvis, x_init[1:15] + xeq[1:15])
maximum(abs.(eigvals(Matrix(A+B*K[1]))))

## Solve LQR
K,P = lqr(Q,R,Qf,A,B,N)
Xlin,Ulin = linear_sim(A,B,K,x_init)
Xlin = SVector{30}.(Xlin)

jinds = sortperm([v[2] for v in values(STATEMAP)])
joint_names = reshape(collect(string.(keys(STATEMAP)[jinds])),1,:)
plot(Xlin, inds=1:15, title="Linear Solution", xlabel="time step", labels=joint_names)
visualize(mvis, Xlin, xeq, fps=Int(floor(N/tf)))
norm(K[1],Inf)

## Simulate the nonlinear dynamica
X = [@SVector zeros(30) for k = 1:N]
U = [@SVector zeros(12) for k = 1:N-1] 
X[1] = x_init + xeq
for k = 1:N-1
    U[k] = clamp.(K[k] * (X[k] - xeq) + ueq, -30, 30)
    X[k+1] = discrete_dynamics(RK4, model, X[k], U[k], 0.0, dt)
end
plot(X, inds=1:15, title="Linear Solution", xlabel="time step", labels=joint_names, ylims=[-2,2])
plot(U, ylims=[-30,30])
norm.(X)

visualize(mvis, X, fps=Int(floor(N/tf)))