module DiscreteTests

using Test
using JitterTime

@testset "Discrete Tests" begin
    A1 = [-1 1; 0 -1]
    B1 = [0; 1]
    C1 = [1 0]
    D1 = 0

    A2 = 0
    B2 = 1
    C2 = -1
    D2 = 0

    R1 = [0 0 0; 0 1 0; 0 0 1]
    R2 = [0 0; 0 0.1]
    Q1 = [1 0 0; 0 2 3; 0 3 0] # (n+p)=⍴Qd
    Q2 = [0 0; 0 1]

    plant = DiscreteSystem(1,A1,B1,C1,D1,2,R1,Q1)
    controller = DiscreteSystem(2,A2,B2,C2,D2,1,R2,Q2)

    N = JitterTimeSystem([plant, controller])
    calcDynamics!(N)

    execSys!(N,1)
    execSys!(N,2)
    execSys!(N,2)
    execSys!(N,1)
    execSys!(N,2)

    @test isapprox(N.dJdt, 0.0)
    @test isapprox(N.J, 0.0)
    @test all(isapprox.(N.Rc, zeros(eltype(N.Rc),size(N.Rc))))
    @test all(isapprox.(N.Ac, zeros(eltype(N.Ac),size(N.Ac))))
    @test all(isapprox.(N.Qc, [1. 0 0 0; 0 2 3 0; 0 3 0 0; 0 0 0 1]))
    @test all(isapprox.(N.P, [1 -1 0 0 0; -1 3.1 0 0 1; 0 0 1 1 0; 0 0 1 1 0; 0 1 0 0 1.1]))

end # testset

end # module
