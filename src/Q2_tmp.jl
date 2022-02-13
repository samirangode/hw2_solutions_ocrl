import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()
using Plots
using Test
const resfile = joinpath(@__DIR__, "Q2.jld2")
include("car.jl");
const isautograder = @isdefined autograder

#   Question 2: TVLQR (25 pts)
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# 
#   In this problem we'll use time-varying LQR (TVLQR) to track a reference
#   trajectory for a simplified model of a car.
# 
#   The Model
#   ===========
# 
#   In this problem we'll be using the standard kinematic "bicycle" model for a
#   car. As a kinematic model, we don't consider the effects of things like tire
#   forces, friction, or aerodynamics forces. The bicycle model combines the
#   tires on each axle into a single tire, and uses simple trigonometric
#   relationships to describe the motion of car. While deriving these equations
#   is good practice, we'll just state them here:
# 
# :$
# 
#   x = \begin{bmatrix} px \ py \ \theta \ \delta \end{bmatrix}, \quad u =
#   \begin{bmatrix} v \ \phi \end{bmatrix}, \quad \dot{x} = \begin{bmatrix} v
#   \cos{(\theta + \beta)} \ v \sin{(\theta + \beta)} \ \frac{v \cos{\beta}
#   \tan{\delta}}{L} \ \phi \end{bmatrix} :$
# 
#   where \theta is the yaw angle, \delta is the steering angle, v is the
#   forward velocity, \phi is the steering angle rate, L is the distance between
#   the wheels, and \beta = \text{atan2}(\delta l_r, L) is the side-slip angle.
#   Here we have defined the x,y position (p_x,p_y) to be relative to the center
#   of mass of the vehicle, located a distance l_r from the rear wheel.
# 
#   In this problem, we use the RobotDynamics.jl package (developed by the REx
#   Lab) to define the model. This package allows some convenient methods to
#   evaluate the discrete dynamics and both the continuous or discrete-time
#   Jacobians automatically using either automatic differentiation or finite
#   differencing. See the code block below for some simple examples of using the
#   API.

# Define the model
model = BicycleModel()

# get the number of states and controls
n = state_dim(model)
m = control_dim(model)
n,m = size(model)  # alternate method

# Evaluate the continuous and discrete dynamics
x0 = SA[0,0,0,0]
u0 = SA[0,0]
t0 = 0.0
dt = 0.1
dynamics(model, x0, u0)
discrete_dynamics(RK4, model, x0, u0, t0, dt)  # use rk4 for integration

# Evaluate the continuous and discrete Jacobians
z0 = KnotPoint(x0,u0,dt,t0)   # create a `KnotPoint` type that stores everything together
∇f = RobotDynamics.DynamicsJacobian(model)
jacobian!(∇f, model, z0)
discrete_jacobian!(RK4, ∇f, model, z0)

# Extract pieces of the Jacobian
A = ∇f.A
B = ∇f.B;

#   The Reference Trajectory
#   ==========================
# 
#   In this problem we'll be tracking a reference trajectory generated via
#   trajectory optimization. The code below loads the reference trajectory,
#   plots it using Plots.jl (and some custom plotting recipes in RobotDynamics)
#   and visualizes it using MeshCat.

# Load trajectories
traj = load(TRAJFILE)
Xref,Uref,tref = traj["X"], traj["U"], traj["times"];

# Plot the states
isautograder || plot(tref, Xref, inds=1:4, labels=["x" "y" "θ" "δ"], legend=:bottomleft, title="states", xlabel="time (s)")

# TIP: Use the `inds` keyword to plot a subset of the states, e.g. `inds=1:2` to only plot the x,y positions

# Plots the controls
isautograder || plot(tref[1:end-1], Uref, labels=["v" "ϕ"], legend=:bottomleft, title="controls", xlabel="time (s)")

# Visualize in MeshCat
if !isautograder
    vis = initialize_visualizer(model)
    render(vis)
end

# Send the trajectory to the visualizer
isautograder || visualize!(vis, model, tref[end], Xref)

#   Part (a): Implement TVLQR (10 pts)
#   ====================================
# 
#   Using the types provided below, implement the method to calculate the
#   feedback gains K to track the provided trajectory

# TASK: Implement the following methods
#       calc_gains!  (5 pts)
#       get_control  (5 pts)
"""
    TVLQR{L,T}

A type that contains all the information needed to evaluate a time-varying LQR 
control policy tracking a trajectory specified by `X`, `U`, and `times`.

# Constructor
    TVLQR(model, X, U, times, Q, R, Qf)

where `model` is a `RobotDynamics.AbstractModel`, `X` and `U` are vectors of the reference
states and controls at times `times`. `Q`, `R` and `Qf` are the cost matrices for TVLQR.

# Methods
The following methods are defined on `TVLQR`:

    get_k(ctrl)
    linearize!(ctrl)
    calc_gains!(ctrl)
    get_control(ctrl, x, t)
"""
struct TVLQR{L,T}
    model::L                     # dynamics model
    X::Vector{Vector{T}}         # state reference trajectory (n,)
    U::Vector{Vector{T}}         # control reference trajectory (m,)
    times::Vector{T}             # times for each point in the trajectory (N,)
    Q::Diagonal{T,Vector{T}}     # state cost matrix for TVLQR (n,n)
    R::Diagonal{T,Vector{T}}     # control cost matrix for TVLQR (m,m)
    Qf::Matrix{T}                # terminal state cost matrix for TVLQR (n,n)
    A::Vector{Matrix{T}}         # discrete state Jacobian for each knot point (n,n)
    B::Vector{Matrix{T}}         # discrete control Jacobian for each knot point (n,m)
    K::Vector{Matrix{T}}         # feedback gain matrices (m,n)
    P::Vector{Matrix{T}}         # cost-to-go (n,n)
end
function TVLQR(model::L, X, U, times::AbstractVector, Q::AbstractMatrix, R::AbstractMatrix, Qf::AbstractMatrix) where L <: AbstractModel
    T = promote_type(eltype(X[1]), eltype(X[2]))
    n,m = size(model)
    N = length(times)
    @assert length(X[1]) == n
    @assert length(U[1]) == m
    @assert length(X) == N 
    @assert N-1 <= length(U) <= N
    A = [zeros(T,n,n)*NaN for k = 1:N-1]
    B = [zeros(T,n,m)*NaN for k = 1:N-1]
    K = [zeros(T,m,n)*NaN for k = 1:N-1]
    P = [zeros(T,n,n)*NaN for k = 1:N]
    TVLQR(model, Vector{T}.(X), Vector{T}.(U), Vector{T}(times), 
        Diagonal{T}(diag(Q)), Diagonal{T}(diag(R)), Matrix{T}(Qf),
        A, B, K, P
    )
end

"""
    get_k(ctrl, t)

Get the time index corresponding to time `t`. 
Useful for implementing zero-order hold control.
Uses binary search to find the time index.
"""
get_k(controller::TVLQR, t) = searchsortedlast(controller.times, t)

"""
    linearize!(ctrl::TVLQR)

Linearize the discretized model about the reference trajectory, storing the result in 
`A` and `B`.
"""
function linearize!(ctrl::TVLQR)
    model = ctrl.model
    N = length(ctrl.X)
    ∇f = RobotDynamics.DynamicsJacobian(model)
    X,U = ctrl.X, ctrl.U
    
    # loop over all the time steps in the reference trajectory
    for k = 1:N-1
        # some boilerplate code...
        dt = ctrl.times[k+1] - ctrl.times[k]
        z = KnotPoint(X[k], U[k], dt, ctrl.times[k])
        
        # evaluate the discrete jacobian at the current time step
        discrete_jacobian!(RK4, ∇f, model, z)
        
        # store the pieces in the controller
        ctrl.A[k] .= ∇f.A
        ctrl.B[k] .= ∇f.B
    end
end

"""
    calc_gains!(ctrl::TVLQR)

Calculate the locally-optimal feedback gains `K` about the current trajectory, 
using the linearized dynamics in `A` and `B`. Should use a Riccati recursion.

**NOTE**: `linearize!(ctrl)` must be called before calling this function!
"""
function calc_gains!(ctrl::TVLQR)
    # Extract some variables
    N = length(ctrl.X)
    A,B = ctrl.A, ctrl.B
    Q,R = ctrl.Q, ctrl.R
    P = ctrl.P
    K = ctrl.K
    
    # TODO: Implement Riccati recursion for TVLQR
    #       After this function, all the matrices in ctrl.K and ctrl.P should be updated
    
    # SOLUTION
    P[end] .= ctrl.Qf
    for k = reverse(1:N-1) 
        K[k] .= (R + B[k]'P[k+1]*B[k])\(B[k]'P[k+1]*A[k])
        P[k] .= Q + A[k]'P[k+1]*A[k] - A[k]'P[k+1]*B[k]*K[k]
    end
    
    # no need to return anything, since the result is stored in TVLQR type
    return nothing
end

"""
    get_control(ctrl::TVLQR, x, t)

Evaluate the TVLQR feedback policy at state `x` and time `t`, returning the control 
to be executed by the system.
"""
function get_control(ctrl::TVLQR, x, t)
    # TODO: implement this function
    #       should return a vector of size (m,), where m is the number of controls
    u = zeros(2)
    
    # SOLUTION
    k = get_k(ctrl, t)
    u = ctrl.U[k] - ctrl.K[k]*(x - ctrl.X[k])
    return u 
end

# LQR Cost weights
Q = Diagonal([1,1,1e-2,1e-2])
R = Diagonal([1e-1,1e-1])
Qf = Diagonal([1,1,1,1.])*10

# Build controller
ctrl = TVLQR(model, Xref, Uref, tref, Q, R, Qf)

# Linearize the model about the trajectory
linearize!(ctrl)

# Calculate the gains using Riccati recursion
calc_gains!(ctrl)

@testset "Q2a" begin   # POINTS = 10
    time = 1.0
    k = get_k(ctrl, time)
    xtest = Xref[k] + [0.2, 0.2, deg2rad(10), deg2rad(-5)]
    @test length(ctrl.K) == 100                   # POINTS = 1
    @test norm(ctrl.K[1] - ctrl.K[end]) > 1       # POINTS = 1
    @test norm(ctrl.P[1] - ctrl.P[end]) > 1       # POINTS = 1
    @test ctrl.K ≈ load(resfile, "K") rtol=1e-3   # POINTS = 3
    @test ctrl.P ≈ load(resfile, "P") rtol=1e-3   # POINTS = 2
    @test get_control(ctrl, xtest, time) ≈ load(resfile, "utest1") rtol=1e-4  # POINTS = 2
end;

#   Part (b): Simulate the system (0 pts)
#   =======================================
# 
#   We'll now simulate our system using our controller, and analyze how well it
#   tracks the reference trajectory under perturbations to the controls and the
#   initial conditions.

"""
    simulate(model, x0, ctrl; [kwargs...])

Simulate `model` starting from `x0` using the `get_control(ctrl, x, t)` method to get the 
closed-loop feedback command.

# Keyword Arguments
* `tf`: final time
* `dt`: simulation time step
* `ν`: standard deviation of the white noise on the controls
* `w` standard deviation of the white noise on the steering angle
"""
function simulate(model::BicycleModel, x0, ctrl; tf=ctrl.times[end], dt=1e-2, ν=0.0, w=0.00)
    n,m = size(model)
    times = range(0, tf, step=dt)
    N = length(times)
    X = [@SVector zeros(n) for k = 1:N] 
    U = [@SVector zeros(m) for k = 1:N-1]
    X[1] = x0

    for k = 1:N-1
        U[k] = get_control(ctrl, X[k], times[k]) + SA[randn(), randn()]*ν
        X[k+1] = discrete_dynamics(RK4, model, X[k], U[k], times[k], dt) + SA[0,0,0,randn()*w]
    end
    return X,U,times
end

"""
    sse(Xref, tref, X, t)

Evaluate the normalized sum-squared error (SSE) of two trajectories of different lengths.
This assumes `length(Xref)` < `length(X)` and that `tref[2]` is a multiple of `t[2]`, 
and that both time vectors have a constant step length.
"""
function sse(Xref,tref,X,t)
    @assert length(tref) < length(t)
    step = tref[2] / t[2]
    if abs(step - round(step)) < 1e-6
        step = Int(round(step))
    end
    inds = 1:step:length(t)
    @assert norm(tref[1:length(inds)] - t[inds]) < 1e-10
    sum(norm.(X[inds] - Xref[1:length(inds)]).^2) / length(inds)
end

# display a new visualizer pane
# TIP: you can also use `open(vis)` to open the visualizer in a tab in your browser (useful if you have multiple monitors)
isautograder || render(vis)

# Simulate the system and compute the errors

# try changing both of these inputs to the simulation!
ν = 0.0            # std of control noise
xinit = [0,0,0,0]  # initial condition
X,U,times = simulate(model, xinit, ctrl, ν=ν)

# compute the errors
err_term = norm(X[end] - Xref[end])
err = sse(Xref,tref,X,times)
@show err_term
@show err

# send trajectory to the visualizer
isautograder || visualize!(vis, model, times[end], X)

#   Let's look at how robust it is to white noise in the input channels,
#   velocity and the steering angle rate. Try changing the standard deviation of
#   the control noise and see how it affects the performance. Use the code below
#   to plot the commands we're sending.

isautograder || plot(times[1:end-1],U, labels=["v" "ϕ"], ylabel="time (s)", title="controls")

#   Part (c): Monte-Carlo Analysis (15 pts)
#   =========================================
# 
#   We'll dig a little deeper into finding out how robust our controller is to
#   the initial condition. Our goal is to find the largest symmetric rectangular
#   region around the origin such that at least 95% of the samples within the
#   region converge to the desired target state.
# 
#   To do this, we'll perform a simple Monte-Carlo analysis. In this case, a
#   deterministic sampling scheme will be more sample-efficient than random
#   sampling. Generate a grid of uniformly-distributed sample points in only x
#   and y, and then simulate the system from each of those initial conditions
#   and check if the terminal error, defined as $ ||X{ref,N} - XN ||_2 $, is
#   less than 0.2.
# 
#   Once you have a grid of boolean values, find a \Delta x and \Delta y such
#   that the area \Delta x \Delta y is maximized and that at least 95% of the
#   initial conditions within \begin{bmatrix} \pm \Delta x & \pm \Delta y & 0 &0
#   \end{bmatrix}^T have a terminal error less than 0.2.
# 
#   HINT: Generate a list of (x,y,success)::Tuple{Float64,Float64,Bool} tuples
#   and use filter and sort to get the information you need. The sort method
#   offers a by argument you may find useful.

# TASK: simulate the initial conditions for all the points in a grid of your choice (3 pts)
res = NTuple{3,Float64}[]   # a data structure you may find useful...

# TODO: pick a grid to sample (hint: start coarse and then refine it)
xmax = 0
ymax = 0
Nx = 0
Ny = 0

# SOLUTION
xmax = 11      # maximum x coordinate to sample
ymax = 2.5     # maximum y coordinate to sample
Nx = 101       # number of sample points in x
Ny = 51        # number of sample points in y

# Generate the ranges
xs = range(-xmax,xmax,length=Nx)
ys = range(-ymax,ymax,length=Ny)
for x in xs, y in ys
    # TODO: simulate the initial conditions and cache the terminal error 
    err_term = NaN
    
    # SOLUTION
    local xinit = [x,y,0,0]
    local X,U,time = simulate(model, xinit, ctrl, ν=0, w=0, dt=0.03)
    local err_term = norm(X[end] - Xref[end])

    # push to res structure
    push!(res, (x,y,err_term))
end
res = hcat(collect.(res)...)  # turn the list of tuples into a 2D matrix

# Plot the initial conditions in a scatter plot
# TODO: fill the following variables to generate the plot (2 pts)
xpoints = zeros(Nx*Ny)         # (Nx*Ny,) vector of x coordinates
ypoints = zeros(Nx*Ny)         # (Nx*Ny,) vector of y coordinates
valid = falses(Nx*Ny)          # (Nx*Ny,) vector of boolean values indicating success

# SOLUTION
xpoints = res[1,:]             # (Nx*Ny,) vector of x coordinates
ypoints = res[2,:]             # (Nx*Ny,) vector of y coordinates
valid = res[3,:] .< 0.2        # (Nx*Ny,) vector of boolean values indicating success

# Generate the plot
isautograder || scatter(xpoints[valid], ypoints[valid], color=:green, xlabel="x", ylabel="y", label="success")
isautograder || scatter!(xpoints[.!valid],ypoints[.!valid], color=:red, label="failure")

# TASK: find the maximum area such that 95% of the points converge to the goal. (10 pts)
#       save the result in Δx and Δy
Δx = NaN
Δy = NaN

# SOLUTION
xabs = unique(abs.(xs))
yabs = unique(abs.(ys))
acc = NTuple{3,Float64}[]
for x in xabs, y in yabs
    inregion = (abs.(res[1,:]) .<= x) .& (abs.(res[2,:]) .<= y)
    num_valid = sum(inregion .& valid)
    rate = num_valid / sum(inregion)
    push!(acc, (x,y,rate))
end
max_region = sort(filter(x->x[3] >= 0.95, acc), by=x->x[1]*x[2])[end]
Δx = max_region[1]   # maximum deviation in x (half of the width of the region)
Δy = max_region[2]   # maximum deviation in y (half of the height of the region)

# Calculate the area
area = Δx * Δy
println("region: Δx = $Δx, Δy = $Δy, area = $area")

isautograder || render(vis)

if !isautograder
    delete!(vis["ic"])
    plot_region!(vis, Δx, Δy)
end

v = (abs.(xpoints) .<= Δx ) .& (abs.(ypoints) .<= Δy)
(sum(valid[v]) / sum(v)) >= 0.95

@testset "Q2c" begin    # POINTS = 15
    @test 2 < Δx < 10   # POINTS = 5
    @test 0.5 < Δy < 2  # POINTS = 5
    v = (abs.(xpoints) .<= Δx) .& (abs.(ypoints) .<= Δy)
    acc = sum(valid[v]) / sum(v)
    @test acc >= 0.95  # POINTS = 3
    @test acc < 1.0    # POINTS = 2
end;