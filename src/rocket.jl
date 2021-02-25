using MeshCat
using CoordinateTransformations
using Rotations 
using Colors
using GeometryBasics
using RobotZoo: PlanarRocket

function set_mesh!(vis, model::PlanarRocket;
    L = model.ℓ,
    rad = L / 15,
    Ltip = L / 5,
    Lcone = L / 10,
    )
    fuselage = Cylinder(Point3(0,0,0.), Point3(0,0,L), rad)
    tip = Cone(Point3(0.,0,L), Point3(0.,0,L+Ltip), rad)
    gimbal = Cone(Point3(0.,0,-Lcone), Point3(0,0,0.), rad/2)
    fins = Pyramid(Point3(0,0,0.), L, 2rad*0.95)
    setobject!(vis["geom"]["fuselage"], fuselage, MeshPhongMaterial(color=colorant"gray"))
    setobject!(vis["geom"]["fins"], fins, MeshPhongMaterial(color=colorant"black"))
    setobject!(vis["geom"]["tip"], tip, MeshPhongMaterial(color=colorant"black"))
    setobject!(vis["geom"]["gimbal"]["cone"], gimbal, MeshPhongMaterial(color=colorant"red"))
    settransform!(vis["geom"]["gimbal"]["cone"], Translation(0,0,Lcone*0.4))
end

function visualize!(vis, model::PlanarRocket, x::StaticArray)
    settransform!(vis["robot"], compose(Translation(x[1] / 50,0,x[2] / 50), LinearMap(RotY(x[3]))))
    T = x[7] / (model.m * model.g)
    settransform!(vis["robot"]["geom"]["gimbal"], LinearMap(RotY(-x[8]*20)*T))
end

"""
    visualize!(vis, model, tf, X)

Visualize a trajectory for `model` give the total time `tf` and a trajectory of states `X`.
"""
function visualize!(vis, model::PlanarRocket, tf::Real, X)
    fps = Int(round((length(X)-1)/tf))
    anim = MeshCat.Animation(fps)
    for (k,x) in enumerate(X)
        atframe(anim, k) do
            x = X[k]
            visualize!(vis, model, SVector{8}(x)) 
        end
    end
    setanimation!(vis, anim)
end


"""
    initialize_visualizer(model)

Launch a MeshCat visualizer and insert the geometry for `model`.
"""
function initialize_visualizer(model::PlanarRocket)
    vis = Visualizer()
    set_mesh!(vis["robot"], model)
    return vis
end

"""
    simulate(model, x0, ctrl; [kwargs...])

Simulate `model` starting from `x0` using the `get_control(ctrl, x, t)` method to get the 
closed-loop feedback command.

# Keyword Arguments
* `tf`: final time
* `dt`: simulation time step
"""
function simulate(model::PlanarRocket, x0, ctrl; tf=ctrl.times[end], dt=1e-2)
    n,m = size(model)
    times = range(0, tf, step=dt)
    N = length(times)
    X = [@SVector zeros(n) for k = 1:N] 
    U = [@SVector zeros(m) for k = 1:N-1]
    X[1] = x0

    tstart = time_ns()
    for k = 1:N-1
        U[k] = get_control(ctrl, X[k], times[k])
#         u = clamp(U[k], umin, umax)
        X[k+1] = discrete_dynamics(RK4, model, X[k], U[k], times[k], dt)
    end
    tend = time_ns()
    rate = N / (tend - tstart) * 1e9
    println("Controller ran at $rate Hz")
    return X,U,times
end

function comparison_plot(model, Z...)
    p = plot(layout=(2,2), size=(1000,800))
    for z in Z
        plot!(p[1], z[3],  z[1],  inds=2:2, label=z[4], 
            xlabel="time (s)", ylabel="altitude (m)", legend=:topright)
        plot!(p[2], z[3],  [rad2deg.(x) for x in z[1]], inds=3:3, label=z[4],
            xlabel="time (s)", ylabel="roll angle (deg)", legend=:none)
        plot!(p[3], z[3],  z[1], inds=7:7, label=z[4],
            xlabel="time (s)", ylabel="Thrust (N)", legend=:none)
        plot!(p[4], z[3],  [rad2deg.(x) for x in z[1]],  inds=8:8, label=z[4],
            xlabel="time (s)", ylabel="Thrust angle (deg)", legend=:none)
    end
    hline!(p[2], [model.max_roll, -model.max_roll], ls=:dash, c=:red, label="limits")
    hline!(p[3], [RobotZoo.umax(model)[1], RobotZoo.umin(model)[1]], ls=:dash, c=:red, label="limits")
    hline!(p[4], rad2deg.([RobotZoo.umax(model)[2], RobotZoo.umin(model)[2]]), ls=:dash, c=:red, label="limits")
    p
end

function run_tests()
    include(joinpath(@__DIR__,"..","test","q3.jl"))
end