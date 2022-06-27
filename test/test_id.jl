module TestID

using Test
using JitterTime

function aux_cont(ids::Vector)
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
    Q2 = zeros(2, 2)

    plant = ContinuousSystem(ids[1], A1, B1, C1, D1, ids[2], R1, [1. 0 0; 0 1 0; 0 0 1])
    controller = ContinuousSystem(ids[2], A2, B2, C2, D2, ids[1], R2, Q2)

    N = JitterTimeSystem([plant, controller])
    calcDynamics!(N)
    passTime!(N, 1000.)
    M = JitterTimeSystem([controller, plant])
    calcDynamics!(M)
    passTime!(M, 1000.)

    Nindices    = Vector{UnitRange}(undef, length(N.systems))
    Mindices    = Vector{UnitRange}(undef, length(M.systems))
    Ntotstates  = 0
    Mtotstates  = 0

    for i in 1:length(N.systems)
        # TODO: Add "stateindex"?
        Nn          = N.systems[i].n
        Mn          = M.systems[i].n
        Nindices[i] = UnitRange(Ntotstates+1, Ntotstates+Nn)
        Mindices[i] = UnitRange(Mtotstates+1, Mtotstates+Mn)
        Ntotstates  += Nn
        Mtotstates  += Mn
    end

    return N, Nindices, M, Mindices
end # function

function aux_disc(ids::Vector)
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

    plant = DiscreteSystem(ids[1],A1,B1,C1,D1,ids[2],R1,Q1)
    controller = DiscreteSystem(ids[2],A2,B2,C2,D2,ids[1],R2,Q2)

    N = JitterTimeSystem([plant, controller])
    calcDynamics!(N)
    execSys!(N,ids[1])
    execSys!(N,ids[2])
    execSys!(N,ids[2])
    execSys!(N,ids[1])
    execSys!(N,ids[2])

    M = JitterTimeSystem([controller, plant])
    calcDynamics!(M)
    execSys!(M,ids[1])
    execSys!(M,ids[2])
    execSys!(M,ids[2])
    execSys!(M,ids[1])
    execSys!(M,ids[2])

    Nindices    = Vector{UnitRange}(undef, length(N.systems))
    Mindices    = Vector{UnitRange}(undef, length(M.systems))
    Ntotstates  = 0
    Mtotstates  = 0

    for i in 1:length(N.systems)
        # TODO: Add "stateindex"?
        Nn          = N.systems[i].n
        Mn          = M.systems[i].n
        Nindices[i] = UnitRange(Ntotstates+1, Ntotstates+Nn)
        Mindices[i] = UnitRange(Mtotstates+1, Mtotstates+Mn)
        Ntotstates  += Nn
        Mtotstates  += Mn
    end

    return N, Nindices, M, Mindices
end # function

@testset "ID test - Pure Continuous (Simple IDs)" begin

    N, Nindices, M, Mindices = aux_cont([1, 2])

    @test isapprox(N.dJdt, M.dJdt)
    @test isapprox(N.J, M.J)
    @test all(isapprox.(N.Ac[Nindices[1], Nindices[1]], M.Ac[Mindices[2], Mindices[2]]))
    @test all(isapprox.(N.Rc[Nindices[1], Nindices[1]], M.Rc[Mindices[2], Mindices[2]]))
    @test all(isapprox.(N.Qc[Nindices[1], Nindices[1]], M.Qc[Mindices[2], Mindices[2]]))
    @test all(isapprox.(N.P[Nindices[1], Nindices[1]], M.P[Mindices[2], Mindices[2]]))
end # testset

@testset "ID test - Pure Continuous (Complex IDs)" begin

    N, Nindices, M, Mindices = aux_cont([49_463, 12_289])

    @test isapprox(N.dJdt, M.dJdt)
    @test isapprox(N.J, M.J)
    @test all(isapprox.(N.Ac[Nindices[1], Nindices[1]], M.Ac[Mindices[2], Mindices[2]]))
    @test all(isapprox.(N.Rc[Nindices[1], Nindices[1]], M.Rc[Mindices[2], Mindices[2]]))
    @test all(isapprox.(N.Qc[Nindices[1], Nindices[1]], M.Qc[Mindices[2], Mindices[2]]))
    @test all(isapprox.(N.P[Nindices[1], Nindices[1]], M.P[Mindices[2], Mindices[2]]))
end # testset

@testset "ID test - Pure Discrete (Simple IDs)" begin
    N, Nindices, M, Mindices = aux_cont([1, 2])

    @test isapprox(N.dJdt, M.dJdt)
    @test isapprox(N.J, M.J)
    @test all(isapprox.(N.Ac[Nindices[1], Nindices[1]], M.Ac[Mindices[2], Mindices[2]]))
    @test all(isapprox.(N.Rc[Nindices[1], Nindices[1]], M.Rc[Mindices[2], Mindices[2]]))
    @test all(isapprox.(N.Qc[Nindices[1], Nindices[1]], M.Qc[Mindices[2], Mindices[2]]))
    @test all(isapprox.(N.P[Nindices[1], Nindices[1]], M.P[Mindices[2], Mindices[2]]))

end # testset

@testset "ID test - Pure Discrete (Complex IDs)" begin
    N, Nindices, M, Mindices = aux_cont([26_561, 64_591])

    @test isapprox(N.dJdt, M.dJdt)
    @test isapprox(N.J, M.J)
    @test all(isapprox.(N.Ac[Nindices[1], Nindices[1]], M.Ac[Mindices[2], Mindices[2]]))
    @test all(isapprox.(N.Rc[Nindices[1], Nindices[1]], M.Rc[Mindices[2], Mindices[2]]))
    @test all(isapprox.(N.Qc[Nindices[1], Nindices[1]], M.Qc[Mindices[2], Mindices[2]]))
    @test all(isapprox.(N.P[Nindices[1], Nindices[1]], M.P[Mindices[2], Mindices[2]]))

end # testset

end # module
