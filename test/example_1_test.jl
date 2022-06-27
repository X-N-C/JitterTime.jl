module Example1Tests

using Test
using JitterTime

@testset "Example 1" begin
    h = 1
    
    A1 = 0
    B1 = 1
    C1 = 1
    D1 = 0
    R1 = 1
    Q1 = [1 0; 0 0] #Punish only the state
    
    K = -(3+sqrt(3))/(2+sqrt(3))/h

    pvec = [0, .100000000000000, 0.200000000000000, 0.300000000000000, 0.400000000000000, 0.500000000000000, 0.600000000000000, 0.700000000000000, 0.800000000000000, 0.900000000000000, 1.000000000000000, 1.000000000000000, 1.100000000000000, 1.200000000000000, 1.300000000000000, 1.400000000000000, 1.500000000000000, 1.600000000000000, 1.700000000000000, 1.800000000000000, 1.900000000000001, 2.000000000000000, 2.000000000000000, 2.100000000000001, 2.200000000000001, 2.300000000000001, 2.400000000000001, 2.500000000000001, 2.600000000000001, 2.700000000000001, 2.800000000000001, 2.900000000000001, 3.000000000000001, 3.000000000000001, 2.387461339178930, 1.871384387633062, 1.451769145362399, 1.128615612366940, 0.901923788646685, 0.771693674201634, 0.737925269031787, 0.800618573137145, 0.959773586517707, 1.215390309173473, 1.215390309173473, 1.026719448082348, 0.877128129211021, 0.766616352559490, 0.695184118127758, 0.662831425915822, 0.669558275923684, 0.715364668151343, 0.800250602598800, 0.924216079266054, 1.087261098153105, 1.087261098153105, 0.929022575872818, 0.805743741577959, 0.717424595268528, 0.664065136944525, 0.645665366605950, 0.662225284252803, 0.713744889885084, 0.800224183502793, 0.921663165105930, 1.078061834694495]
    Jvec = [0, 0.005000000000000, 0.020000000000000, 0.045000000000000, 0.080000000000000, 0.125000000000000, 0.180000000000000, 0.245000000000000, 0.320000000000000, 0.405000000000000, 0.500000000000000, 0.500000000000000, 0.605000000000000, 0.720000000000000, 0.845000000000000, 0.980000000000000, 1.125000000000000, 1.280000000000000, 1.445000000000000, 1.620000000000000, 1.805000000000000, 2.000000000000000, 2.000000000000000, 2.205000000000000, 2.420000000000000, 2.645000000000000, 2.880000000000000, 3.125000000000000, 3.380000000000000, 3.645000000000000, 3.920000000000000, 4.205000000000000, 4.500000000000000, 4.500000000000000, 4.768569219381654, 4.980707658144960, 5.146061487217439, 5.274276877526613, 5.375000000000001, 5.457877025565123, 5.532554125149501, 5.608677469680654, 5.695893230086103, 5.803847577293369, 5.803847577293369, 5.915627402304328, 6.010494118317165, 6.092355679553859, 6.165120040236389, 6.232695154586737, 6.298988976826880, 6.367909461178800, 6.443364561864476, 6.529262233105887, 6.629510429125013, 6.629510429125013, 6.730033282093097, 6.816480267232424, 6.892347353341536, 6.961130509218977, 7.026325703663288, 7.091428905473014, 7.159936083446697, 7.235343206382879, 7.321146243080103, 7.420841162336913]
    pvecjl = Float64[]
    Jvecjl = Float64[]

    plant = ContinuousSystem(10,A1,B1,C1,D1,30,R1,Q1)
    sampler = DiscreteSystem(20,1,10)
    controller = DiscreteSystem(30,K,20)

    N = JitterTimeSystem([plant, sampler, controller])
    calcDynamics!(N)

    Nsteps = 6
    dt = h/10.0
    for k = 1:Nsteps
        for j = 1:10
            push!(pvecjl, N.P[1, 1])
            push!(Jvecjl, N.J)
            passTime!(N,dt)
        end
        push!(pvecjl, N.P[1, 1])
        push!(Jvecjl, N.J)
        if k >= 3
            execSys!(N,20)
            execSys!(N,30)
        end
    end
    
    @test isapprox(N.dJdt, 0.996949192568094)
    @test isapprox(N.J, 7.420841162336913)
    #@test all(isapprox.(N.Rc, [1. 0 0; 0 0 0; 0 0 0]))
    #@test all(isapprox.(N.Ac, [0. 0 1; 0 0 0; 0 0 0]))
    #@test all(isapprox.(N.Qc, [1. 0 0; 0 0 0; 0 0 0]))
    @test all(isapprox.(N.P, [1.078061834694495 1.078061834694495 -1.366927632691700; 1.078061834694495 1.078061834694495 -1.366927632691700; -1.366927632691700 -1.366927632691700 1.733194787983227]))
    @test all(isapprox.(pvecjl, pvec))
    @test all(isapprox.(Jvecjl, Jvec))

    reset!(N)
    for k = 1:Nsteps
        for j = 1:10
            passTime!(N,dt)
        end
        if k >= 3
            execSys!(N,20)
            execSys!(N,30)
        end
    end
    
    @test isapprox(N.dJdt, 0.996949192568094)
    @test isapprox(N.J, 7.420841162336913)
    #@test all(isapprox.(N.Rc, [1. 0 0; 0 0 0; 0 0 0]))
    #@test all(isapprox.(N.Ac, [0. 0 1; 0 0 0; 0 0 0]))
    #@test all(isapprox.(N.Qc, [1. 0 0; 0 0 0; 0 0 0]))
    @test all(isapprox.(N.P, [1.078061834694495 1.078061834694495 -1.366927632691700; 1.078061834694495 1.078061834694495 -1.366927632691700; -1.366927632691700 -1.366927632691700 1.733194787983227]))

end # testset

end # module