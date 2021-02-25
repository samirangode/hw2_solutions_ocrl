using RobotDynamics
using RobotZoo
using TrajOptPlots
using StaticArrays
using LinearAlgebra
using SparseArrays

function visualize(mvis, X, xeq=zero(X[1]); fps=20)
    anim = MeshCat.Animation(fps)
    for (k,x) in enumerate(X)
        atframe(anim, k) do
            set_configuration!(mvis, x[1:15] + xeq[1:15])
        end
    end
    setanimation!(mvis, anim)
end

## Generate Model and equilibrium point
model = RobotZoo.Cartpole()
xeq = SA[0,pi,0,0]
ueq = SA[0.]
dt = 0.01
z = KnotPoint(xeq,ueq,dt)
norm(dynamics(model, xeq, ueq))

## Setup LQR Problem
∇f = RobotDynamics.DynamicsJacobian(model)
discrete_jacobian!(RK4, ∇f, model, z)
A = RobotDynamics.get_static_A(∇f)
B = RobotDynamics.get_static_B(∇f)
Q = Diagonal(@SVector [10,10,.1,.1])
R = Diagonal(@SVector fill(1e-2, 1))
Qf = Diagonal(@SVector fill(100.0, 4))
tf = 2.0
N = Int(round(tf/dt)) + 1
x_init = SA[0.1,0.1,0,0]

## Solve LQR
K,P = lqr(Q,R,Qf,A,B,N)
Xlin,Ulin = linear_sim(A,B,K,x_init)
Xlin = SVector{4}.(X)
plot(Xlin, inds=1:2, title="Linear Solution", labels=["Δx" "Δθ"], xlabel="time step")

## Simulate the nonlinear dynamica
X = [@SVector zeros(4) for k = 1:N] 
X[1] = [0,0.4,0.0,0] + xeq
for k = 1:N-1
    u = K[k] * (X[k] - xeq) + ueq
    X[k+1] = discrete_dynamics(RK4, model, X[k], SVector{1}(u), 0.0, dt)
end
plot(X, labels=["x" "θ" "v" "ω"], xlabel="time step", title="Nonlinear dynamics", legend=:bottomright)

Kmat = vcat(K...)
plot(Kmat, legend=:bottomleft, labels=["Kx" "Kθ" "Kv" "Kω"], title="Feedback gains", xlabel="time steps")

## Visualization
using MeshCat
vis = Visualizer()
open(vis)
TrajOptPlots.set_mesh!(vis, model)
visualize!(vis, model, tf, X)
