using Test
using FileIO
using JLD2
data = load(joinpath(@__DIR__,"q2.jld2"))

@testset "Q2" begin
@testset "Part a" begin
# Part a
@test length(ctrl.K) == 100
@test norm(ctrl.K[1] - ctrl.K[end]) > 1
@test norm(ctrl.P[1] - ctrl.P[end]) > 1
@test ctrl.K[1] ≈ data["K1"] rtol=1e-3
@test ctrl.P[1] ≈ data["P1"] rtol=1e-3

for k = 1:length(Xref)-1
    @test norm(get_control(ctrl, Xref[k], tref[k]) - Uref[k]) ≈ 0
end
xtest = Xref[1] + randn(4)
@test norm(get_control(ctrl, xtest, tref[1]) - Uref[1]) > 1e-1 
@test norm(get_control(ctrl, xtest, tref[1]) - get_control(ctrl, xtest, tref[10])) > 1e-1 
end

# Part b 
# No points in this section 

@testset "Part c" begin
# Part c 
@test 2 < Δx < 10
@test 0.5 < Δy < 2
v = (abs.(xpoints) .<= Δx ) .& (abs.(ypoints) .<= Δy)
acc = sum(valid[v]) / sum(v)
@test acc >= 0.95
@test acc < 1.0
end
end