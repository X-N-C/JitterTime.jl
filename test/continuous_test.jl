module ContinuousTests

using Test
using JitterTime

@testset "Continuous Tests" begin

    A1 = [-1 1; 0 -1]
    B1 = [0; 1]
    C1 = [1. 0]
    D1 = 0

    A2 = 0.
    B2 = 1
    C2 = -1
    D2 = 0

    R1 = B1*B1'
    R2 = 1.
    Q1 = C1'*C1
    Q2 = 1.
    Q2 = zeros(Float64, 2, 2)

    plant = ContinuousSystem(1, A1, B1, C1, D1, 2, R1, [1. 0 0; 0 1 0; 0 0 1])
    controller = ContinuousSystem(2, A2, B2, C2, D2, 1, R2, Q2)

    N = JitterTimeSystem([plant, controller])
    calcDynamics!(N)
    passTime!(N, 1000)

    @test isapprox(N.dJdt, 7.47)
    @test isapprox(N.J, 7470.)
    @test all(isapprox.(N.Rc, [0 0 0; 0 1 0; 0 0 1]))
    @test all(isapprox.(N.Ac, [-1 1 0; 0 -1 -1; 1 0 0]))
    @test all(isapprox.(N.Qc, [1 0 0; 0 1 0; 0 0 1]))
    @test all(isapprox.(N.P, [1.5 1.5 -0.5; 1.5 2.5 -2.0; -0.5 -2 3.5]))

end # testset

end # module
