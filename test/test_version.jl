module TestVersion

using JitterTime
using Test

@testset "Test Versions" begin
    h = 1.0
    
    A1 = 0
    B1 = 1
    C1 = 1
    D1 = 0
    R1 = 1
    Q1 = [1 0; 0 1]

    D3 = -(3+sqrt(3))/(2+sqrt(3))/h


    plant = ContinuousSystem(10,A1,B1,C1,D1,30,R1,Q1)
    sampler = DiscreteSystem(20,1,10)
    sampler2 = DiscreteSystem(20,0,10)
    samplerver = VersionedSystem([sampler, sampler2])
    controller = DiscreteSystem(30,D3,20)

    N = JitterTimeSystem([plant, samplerver, controller])
    calcDynamics!(N)

    dt = h/10.0
    for k = 1:6
        for j = 1:10
            passTime!(N,dt)
        end
        if k == 3
            execSys!(N,20,1)
            execSys!(N,30)
        elseif k > 3
            execSys!(N,20,2)
            execSys!(N,30)
        end
    end

    @test all(isapprox.(N.J, 15.057713659400530))
    @test all(isapprox.(N.dJdt, 3.165390309173476))
    @test all(isapprox.(N.Ac, [0 1; 0 0]))
    @test all(isapprox.(N.P, [3.215390309173475 0 0; 0 0 0; 0 0 0]))
    @test all(isapprox.(N.Qc, [1 0; 0 1]))
    @test all(isapprox.(N.Rc, [1 0; 0 0]))
    samp_idx = N.idtoindex[20]
    @test all(isapprox.(N.Ad[samp_idx][1], [1 0 0; 1 0 0; 0 0 1]))
    @test all(isapprox.(N.Ad[samp_idx][2], [1 0 0; 0 0 0; 0 0 1]))

end # testset

end # module
