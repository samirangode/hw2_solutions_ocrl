using Test
using Statistics
using ControlSystems: dlqr, dare
using BlockArrays

@testset "Q3" begin

@testset "Part a" begin
# Part a
xdiff = diff(Xref)[1:end-1]
@test mean(xdiff)[1] ≈ 200 / (N-1) rtol = 1e-6
@test mean(xdiff)[2] ≈ -500 / (N-1) rtol = 1e-6
@test mean(xdiff)[3:end] ≈ zeros(6) rtol = 1e-6

Xref2 = nominal_trajectory([100,200,0, 0,0,0, 12.2, 0.], 101, 0.15)
xmean = mean(Xref2[1:end-1])
@test xmean[4] ≈ -100 / 100 / .15 rtol = 1e-6
@test xmean[5] ≈ -200 / 100 / .15 rtol = 1e-6
@test xmean[6] ≈ 0 
@test xmean[7] ≈ 12.2 
@test xmean[8] ≈ 0 
end

@testset "Part b" begin
# Part b
@test ctrl.K ≈ dlqr(A,B,Q,R) rtol = 1e-3
@test Qf ≈ dare(A,B,Q,R) rtol = 1e-3
xtest = xeq + randn(8)
for k = 1:N-1
    @test get_control(ctrl, Xref[k], tref[k]) ≈ zeros(2)
end
@test !(get_control(ctrl, xtest, tref[1]) ≈ get_control(ctrl, xtest, tref[end]))
@test get_control(ctrl, Xref[10] + xtest, tref[10]) ≈ 
    get_control(ctrl, Xref[end] + xtest, tref[end])

@test norm(Xlqr[end][4:6]) < 0.1
@test norm(Xlqr[end][1:3]) < 0.1
@test minimum([x[2] for x in Xlqr]) < -10
@test maximum([x[7] for x in Xlqr]) > RobotZoo.umax(model)[1] 
@test maximum([abs(x[3]) for x in Xlqr]) > deg2rad(model.max_roll)
end

@testset "Part c" begin
# Part c
Np = 50*(n+m)
Nd1 = 50*n
@test size(mpc1.P) == (Np,Np) 
size(mpc1.A) == (Nd1, Np)
parts1 = fill(n,50)
parts2 = repeat([m,n],50)
D = PseudoBlockArray(mpc1.A, parts1, parts2)
@test D[Block(1,1)] ≈ B
@test D[Block(2,1)] ≈ zero(B)
@test D[Block(2,2)] ≈ A
@test D[Block(2,4)] ≈ -I(n) 
@test D[Block(3,4)] ≈ A
@test D[Block(3,5)] ≈ B
@test D[Block(3,6)] ≈ -I(n) 

H = PseudoBlockArray(mpc1.P, parts2, parts2)
@test H[Block(1,1)] ≈ R
@test H[Block(2,2)] ≈ Q
@test H[Block(100,100)] ≈ Qf
get_control(mpc1, Xref[1], tref[1])
g = PseudoBlockArray(mpc1.q, parts2)
@test g[Block(1)] ≈ zeros(2)
@test g[Block(2)] ≈ -Q*(Xref[2] - xeq)
@test g[Block(3)] ≈ zeros(2)
@test g[Block(4)] ≈ -Q*(Xref[3] - xeq)
@test g[Block(100)] ≈ -Qf*(Xref[51] - xeq)

get_control(mpc1, Xref[1], tref[10])
@test g[Block(2)] ≈ -Q*(Xref[11] - xeq)
@test g[Block(4)] ≈ -Q*(Xref[12] - xeq)

get_control(mpc1, Xref[1], tref[end-10])
@test g[Block(2)] ≈ -Q*(Xref[242] - xeq)
@test g[Block(2*9)] ≈ -Q*(Xref[250] - xeq)
@test g[Block(2*11)] ≈ zeros(8)
@test mpc1.lb[1:n] ≈ -A*(Xref[1] - xeq)

@test norm(Xmpc1[end][4:6]) < 1e-2 
@test norm(Xmpc1[end][1:3]) < 1e-2 
@test minimum([x[2] for x in Xmpc1]) < 0 
@test minimum([x[2] for x in Xmpc1]) > -5
@test maximum([x[7] for x in Xmpc1]) - RobotZoo.umax(model)[1]  > 10
@test maximum([abs(x[3]) for x in Xmpc1]) - deg2rad(model.max_roll) > 1e-2
end

@testset "Part d" begin
# Part d
Np = 50*(n+m)
Nd2 = 50*(n+4)
@test size(mpc2.A) == (Nd2,Np)

@test norm(Xmpc2[end][4:6]) < 1e-2
@test norm(Xmpc2[end][1:3]) < 1e-2 
@test minimum([x[2] for x in Xmpc2]) < 0 
@test minimum([x[2] for x in Xmpc2]) > -0.1 
@test (maximum([x[7] for x in Xmpc2]) - RobotZoo.umax(model)[1]) < 1e-3
@test maximum([abs(x[3]) for x in Xmpc2]) - deg2rad(model.max_roll) < 1e-2
end

# Part e
# No Points
end