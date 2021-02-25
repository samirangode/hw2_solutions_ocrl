using Test
using FileIO

@testset "Q1" begin
xtest = xeq + randn(28)*0.1
data = load(joinpath(@__DIR__, "q1.jld2"))

@testset "Part a" begin
# Part a
K1 = data["K1"] 
P1 = data["P1"] 
@test length(K) == 200
@test size(K[1]) == (12,28)
@test size(P[1]) == (28,28)
@test K[1] ≈ K1 rtol = 1e-4
@test P[1] ≈ P1 rtol = 1e-4
@test norm(K[1],Inf) >  norm(K[end],Inf)
end

@testset "Part b" begin
# Part b 
@test norm(get_control(ctrl, xeq, 0.0) - ueq) ≈ 0
@test norm(get_control(ctrl, xeq, 1.0) - ueq) ≈ 0
@test norm(get_control(ctrl, xtest, 0.0) - ueq)  > 0
@test norm(get_control(ctrl, xtest, 0.0) - get_control(ctrl, xtest, 1.0)) > 1
@test get_control(ctrl, xtest, 0.0) - ueq ≈ -K[1]*(xtest - xeq) rtol=1e-4
end

@testset "Part c" begin
# Part c
@test length(Kinf) >= 2*length(K)
@test Kinf[1] ≈ data["Kinf"] rtol = 1e-3
@test norm(Kinf[1] - Kinf[2]) < 1e-2
@test norm(Pinf[1] - Pinf[2]) < 1e2
@test norm(get_control(ctrl, xeq, 0.0) - ueq) ≈ 0
@test norm(get_control(ctrl_inf, xtest, 0.0) - get_control(ctrl_inf, xtest, 1.0)) ≈ 0
end

@testset "Part d" begin
# Part d
@test stability0 > 1
@test stability < stability0
@test stability < 1
end

@testset "Part e" begin
# Part e
@test norm(X[end] - xeq) < 0.2
@test norm(Xinf[end] - xeq) < 0.1
@test norm(Xinf[end] - xeq) < norm(X[end] - xeq)
end

end