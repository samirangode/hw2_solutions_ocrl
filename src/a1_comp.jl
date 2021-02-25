import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()
include("../../quadruped-control/src/quadruped_control.jl")
include("quadruped.jl")
using .quadruped_control

function ev_dynamics(x::AbstractVector{T}, u::AbstractVector{F}) where {T,F}
    if typeof(x[1]) <: ForwardDiff.Dual
#         state = MechanismState{T}(a1.mechanism);
        a1.state = statecache[T]
        dr = dr_cache[T];
    else
#         state = MechanismState{F}(a1.mechanism);
        a1.state = statecache[F]
        dr = dr_cache[F]
        x = convert(AbstractVector{F},x);
    end
    ẋ = zeros(typeof(x[1]),size(x,1));
    u_total = vcat(SVector(0.,0.,0.),u)
    ẋ = dynamics!(ẋ,dr,a1.state,x,u_total)
    return ẋ
end

## Old Model
a1 = SingleLegPendulum(load_a1_urdf("../quadruped-control/src/a1/urdf/a1.urdf"));
build_rear_foot_constraints_revolute!(a1);
data = load("../quadruped-control/Experiments/single_leg_init_revolute.jld")
statecache = StateCache(a1.mechanism)
dr_cache = DynamicsResultCache(a1.mechanism)

## New Model
model = UnitreeA1()

## Compare Dynamics
x0 = initial_state(model)
u0 = zeros(12)
dynamics(model, x0, u0)
# (dynamics(model, x0, u0) - ev_dynamics(x0, u0))[16:end]
# q_guess = [0.7494121122172587, -0.30800825827034156, 0.0739469451153811, -0.8973, 1.106, -0.3587759235809219, -0.8007, 0.8020,0.158,0.1667,0.482,1.2031,-0.9088,-0.9091,-0.906814375366219]
# x_guess = [q_guess; zeros(15)]

x0 = zeros(14)
x0[2] += pi/4
x0 = initial_state(model)
set_configuration!(mvis, x0[1:14])
xeq, ueq = newton_solve(model, x0, R=1e-3, ρ=1e-3)
norm(dynamics(model, xeq, ueq))
# norm(dynamics(model, data["x0"], data["u0"]))
# norm(ev_dynamics(data["x0"], data["u0"]))
# xeq, ueq = data["x0"], data["u0"]

∇f = jacobian(model, xeq, ueq)
# norm(∇f.A - data["A"])
# norm(∇f.B - data["B"])
A = ∇f.A
B = ∇f.B


## Test LQR Controller
using ControlSystems
Q = 130 * Diagonal(ones(size(x0)));
R = 10 * Diagonal(ones(size(u0)));
K = ControlSystems.lqr(A,B,Q,R)
maximum(real(eigvals(A-B*K)))

function lqr_controller!(torques::AbstractVector, t, state::MechanismState)
    q = configuration(state)
    v = velocity(state)
    x = vcat(q,v)
    torques[3:end] .= ueq - K*(x - xeq)
    torques[1:2].= 0 #
end

##
x_init = copy(xeq)
x_init[18] -= 0.1
x_init[19] -= 0.05
# x_init[6] = x_init[6] - 0.1
# x_init[7] = x_init[7] + 0.1
# x_init[8] = x_init[8] + 0.1
# x_init[9] = x_init[9] - 0.1
# x_init[10] = x_init[10] - 0.1
# x_init[12] = x_init[12] + 0.3
# x_init[13] = x_init[13] + 0.3
# x_init[14] = x_init[14] + 0.3
# x_init[29] -= 2
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

##
animation = MeshCat.Animation(mvis, t_lqrs, q_lqrs);
setanimation!(mvis, animation);


## Visualizer
mvis = initialize_visualizer(model)

open(mvis)
