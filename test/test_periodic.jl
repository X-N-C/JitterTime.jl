module TestPeriodic

using JitterTime
using Test

@testset "Periodic Tests" begin

    p_id    = 10
    s_id    = 20
    c_id    = 30
    h       = 1

    plant   = ContinuousSystem(p_id, 0, 1, 1, 0, c_id, 1, [1 0; 0 0])
    sampler = DiscreteSystem(s_id, 1, p_id)
    ctrler  = DiscreteSystem(c_id, -(3+sqrt(3))/(2+sqrt(3))/h, s_id)
    N       = JitterTimeSystem([plant, sampler, ctrler])
    Np      = PeriodicJitterTimeSystem(N)

    @test_throws ErrorException PeriodicJitterTimeSystem(Np)

    calcDynamics!(Np)
    beginPeriodicAnalysis!(Np)
    for _ in 1:10
        passTime!(Np, 0.1)
    end # for
    execSys!(Np, s_id)
    execSys!(Np, c_id)
    endPeriodicAnalysis!(Np)

    @test isapprox(Np.J, 0.5)
    @test isapprox(Np.dJdt, 0.95)
    @test all(isapprox.(Np.Ac, [0 1; 0 0]))
    @test all(isapprox.(Np.Rc, [1 0; 0 0]))
    @test all(isapprox.(Np.Pper, [1.077350269189626 1.077350269189626 -1.366025403784438; 1.077350269189626 1.077350269189626 -1.366025403784438; -1.366025403784438 -1.366025403784438 1.732050807568876]))
    @test all(isapprox.(Np.mper, [0; 0; 0]))
    @test all(isapprox.(Np.Atot, [1.0 0 1.0; 1.0 0 1.0; -1.267949192431123 0 -1.267949192431122]))
    @test all(isapprox.(Np.Rtot, [1.0 1.0 -1.267949192431122; 1.0 1.0 -1.267949192431122; -1.267949192431122 -1.267949192431122 1.607695154586736]))
    @test all(isapprox.(Np.dtot, [0; 0; 0]))

    reset!(Np)
    beginPeriodicAnalysis!(Np)
    for _ in 1:10
        passTime!(Np, 0.1)
    end # for
    execSys!(Np, s_id)
    execSys!(Np, c_id)
    endPeriodicAnalysis!(Np)

    @test isapprox(Np.J, 0.5)
    @test isapprox(Np.dJdt, 0.95)
    @test all(isapprox.(Np.Ac, [0 1; 0 0]))
    @test all(isapprox.(Np.Rc, [1 0; 0 0]))
    @test all(isapprox.(Np.Pper, [1.077350269189626 1.077350269189626 -1.366025403784438; 1.077350269189626 1.077350269189626 -1.366025403784438; -1.366025403784438 -1.366025403784438 1.732050807568876]))
    @test all(isapprox.(Np.mper, [0; 0; 0]))
    @test all(isapprox.(Np.Atot, [1.0 0 1.0; 1.0 0 1.0; -1.267949192431123 0 -1.267949192431122]))
    @test all(isapprox.(Np.Rtot, [1.0 1.0 -1.267949192431122; 1.0 1.0 -1.267949192431122; -1.267949192431122 -1.267949192431122 1.607695154586736]))
    @test all(isapprox.(Np.dtot, [0; 0; 0]))
end # testset

end # module
