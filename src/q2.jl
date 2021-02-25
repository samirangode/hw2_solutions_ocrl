struct TVLQR{L,T}
    model::L
    X::Vector{Vector{T}}
    U::Vector{Vector{T}}
    times::Vector{T}
    Q::Diagonal{T,Vector{T}}
    R::Diagonal{T,Vector{T}}
    Qf::Matrix{T}
    A::Vector{Matrix{T}}
    B::Vector{Matrix{T}}
    K::Vector{Matrix{T}}
    P::Vector{Matrix{T}}
end
function TVLQR(model::L, X, U, times::AbstractVector, Q::AbstractMatrix, R::AbstractMatrix, Qf::AbstractMatrix) where L <: AbstractModel
    T = promote_type(eltype(X[1]), eltype(X[2]))
    n,m = size(model)
    N = length(times)
    @assert length(X[1]) == n
    @assert length(U[1]) == m
    @assert length(X) == N 
    @assert N-1 <= length(U) <= N
    A = [zeros(T,n,n) for k = 1:N-1]
    B = [zeros(T,n,m) for k = 1:N-1]
    K = [zeros(T,m,n) for k = 1:N-1]
    P = [zeros(T,n,n) for k = 1:N]
    TVLQR(model, Vector{T}.(X), Vector{T}.(U), Vector{T}(times), 
        Diagonal{T}(diag(Q)), Diagonal{T}(diag(R)), Matrix{T}(Qf),
        A, B, K, P
    )
end

function linearize!(ctrl::TVLQR)
    model = ctrl.model
    N = length(ctrl.X)
    ∇f = RobotDynamics.DynamicsJacobian(model)
    X,U = ctrl.X, ctrl.U
    for k = 1:N-1
        dt = ctrl.times[k+1] - ctrl.times[k]
        z = KnotPoint(X[k], U[k], dt)
        discrete_jacobian!(RK4, ∇f, model, z)
        ctrl.A[k] .= ∇f.A
        ctrl.B[k] .= ∇f.B
    end
end

function calc_gains!(ctrl::TVLQR)
    N = length(ctrl.X)
    A,B = ctrl.A, ctrl.B
    Q,R = ctrl.Q, ctrl.R
    P = ctrl.P
    K = ctrl.K
    P[end] .= ctrl.Qf
    for k = reverse(1:N-1) 
        K[k] .= (R + B[k]'P[k+1]*B[k])\(B[k]'P[k+1]*A[k])
        P[k] .= Q + A[k]'P[k+1]*A[k] - A[k]'P[k+1]*B[k]*K[k]
    end
    return nothing
end

function get_control(ctrl::TVLQR, x, t)
    k = get_k(ctrl, t)
    return ctrl.U[k] - ctrl.K[k]*(x - ctrl.X[k])
end
get_k(controller::TVLQR, t) = searchsortedlast(controller.times, t)

##
model = BicycleModel()
data = load(TRAJFILE)
X,U,times = data["X"], data["U"], data["times"]

Q = Diagonal([1,1,1e-2,1e-2])
R = Diagonal([1e0,1e0])
Qf = Diagonal([1,1,1,10.])

ctrl = TVLQR(model, X, U, times, Q, R, Qf)
linearize!(ctrl)
calc_gains!(ctrl)