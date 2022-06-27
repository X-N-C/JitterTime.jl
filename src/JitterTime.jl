module JitterTime

using LinearAlgebra

include("types.jl")             # Containing all datastructures necessary for JitterTime
include("util.jl")              # Containing auxiliary functions
include("printing.jl")          # Containing all Pretty Printing code

export calcDynamics!, 
       addSys!, 
       execSys!, 
       passTime!, 
       passTimeUntil!, 
       reset!,
       beginPeriodicAnalysis!,
       endPeriodicAnalysis!

"""
    `calcDynamics!(N)`

Calculate the total system dynamics for the JitterTime system N. 

__Arguments__:
N       The JitterTime system.

__State dynamics__:
Let the state of the total system (including all subsystems) be x.

__Continous time evolution__:
dx(t) = N.Ac x(t)dt + dvc(t)  where vc is continuous white noise with
                              intensity N.Rc

__Discrete-event evolution when executing system S__:
x+(t_k) = S.Ad x(t_k) + vd(t_k) where vd is white Gaussian noise with
                                variance S.Rd
"""
function calcDynamics!(N::JTSystem{T}) where {T}
    # Step 1: Count states and inputs and build a system-to-state index mapping
    indices     = Vector{UnitRange}(undef, length(N.systems))
    totstates   = 0

    empty!(N.idtoindex)
    merge!(N.idtoindex,
           Dict(N.systems[x].id => x for x in eachindex(N.systems)))
    length(N.idtoindex) == length(N.systems) || error("Duplicate System ID's")

    for i in eachindex(N.systems)
        # TODO: Add "stateindex"? Used for stateDisturbance
        n           = N.systems[i].n
        indices[i]  = UnitRange(totstates+1, totstates+n)
        totstates   += n
    end # for

    # Step 2: Check that the number of outputs and inputs match for all system connections (including reset dynamics)
    # TODO Add reset dynamics
    # TODO Add continuous systems
    # FIXME Bug in original Matlab code? Can't use an arbitrary ID
    for s in N.systems
        if s isa VersionedSystem
            for v in s.versions
                sum(N.systems[N.idtoindex[x]].p for x in v.inputid) == v.r || error("Wrong number of inputs.")
            end # for
        else
            sum(N.systems[N.idtoindex[x]].p for x in s.inputid) == s.r || error("Wrong number of inputs.")
        end # if
    end # for

    # Initialice Ac, Rc, Qc
    Ac = zero(N.P)
    Rc = zero(N.P)
    Qc = zero(N.P)

    # Step 3: Formulate the total continuous-time dynamics, continuous noise, discrete-time dynamics (all versions), discrete noise and continuous cost
    # TODO Add reset dynamics?
    for (k, s) in enumerate(N.systems)
        xtou = zeros(T, s.n, totstates)
        xtou[:, indices[k]] = Matrix{T}(I, s.n, s.n) # = I
        if s isa ContinuousSystem
            # Construct Ac matrix
            Ac[indices[k], indices[k]] = s.A
            bix = 1

            for inputindex in map(x -> N.idtoindex[x], s.inputid)
                input = N.systems[inputindex]
                # Construct Ac matrix
                Ac[indices[k], indices[inputindex]] = s.B[:, bix:bix+input.p-1]*input.C
                bix += input.p
                xtoy = zeros(T, input.p, totstates)
                xtoy[:, indices[inputindex]] = input.C
                xtou = vcat(xtou, xtoy)
            end # for

            # Construct Rc matrix
            Rc[indices[k], indices[k]] = s.Rc
            
        elseif s isa DiscreteSystem # TODO Possible to combine this and VersionedSystem?
            sidx = N.idtoindex[s.id]
            push!(N.Ad[sidx], Matrix{T}(I, totstates, totstates)) # = I
            N.Ad[sidx][1][indices[k], indices[k]] = s.A
            push!(N.Rd[sidx], zeros(T, totstates, totstates))
            N.Rd[sidx][1][indices[k], indices[k]] = s.R
            
            bix = 1
            for inputindex in map(x -> N.idtoindex[x], s.inputid)
                input = N.systems[inputindex]
                N.Ad[sidx][1][indices[k], indices[inputindex]] = s.B[:, bix:bix+input.p-1]*input.C
                bix += input.p
            end # for

        elseif s isa VersionedSystem
            for (v, sver) in enumerate(s.versions)
                sidx = N.idtoindex[s.id]
                push!(N.Ad[sidx], Matrix{T}(I, totstates, totstates)) # = I
                N.Ad[sidx][v][indices[k], indices[k]] = sver.A
                push!(N.Rd[sidx], zeros(T, totstates, totstates))
                N.Rd[sidx][v][indices[k], indices[k]] = sver.R
                
                bix = 1
                for inputindex in map(x -> N.idtoindex[x], sver.inputid)
                    input = N.systems[inputindex]
                    N.Ad[sidx][v][indices[k], indices[inputindex]] = sver.B[:, bix:bix+input.p-1]*input.C
                    bix += input.p
                end # for
            end # for
        end # if
        
        # Construct Qc matrix
        Qc .+= xtou'*s.Qc*xtou # Combined cost of state, outputs and inputs
    end # for

    N.reduced = ReducedSystem(Ac, Rc, Qc)

    return nothing
end

"""
    `execSys!(N, sysid [, ver])`

Execute the discrete-time system 'sysid'

__Arguments__:
N           The JitterTime system.
sysid       ID of the discrete-time system to be updated.

__Optional Arguments__:
ver         What version to execute (default = 1)
"""
#function execSys!(N::T, sysid::SysID, ver::Integer = 1) where {T <: JitterTimeSystem}
#    #sidx = N.idtoindex[sysid]
#    sidx = getfield(N,:idtoindex)[sysid]
#    #N.systems[sidx] isa DiscreteSystem || N.systems[sidx] isa VersionedSystem || error("Can only execute (versioned) discrete-time systems.")
#    getfield(N,:systems)[sidx] isa Union{DiscreteSystem, VersionedSystem} || error("Can only execute (versioned) discrete-time systems.")
#
#    #Ad = N.Ad[sidx][ver]
#    Ad = getfield(N,:Ad)[sidx][ver]
#    #Rd = N.Rd[sidx][ver]
#    Rd = getfield(N,:Rd)[sidx][ver]
#
#    #N.m .= Ad * N.m
#    setfield!(N,:m, Ad*getfield(N,:m))
#
#    # N.P = Ad * N.P * Ad' .+ Rd
#    #_mat_mul_add!(Ad, N.P, Ad', Rd, N.bufl)
#    _mat_mul_add!(Ad, getfield(N,:P),  Ad', Rd, getfield(N,:reduced).bufl)
#
#    return nothing
#end

execSys!(N::FixedIntervalJitterTimeSystem, sysid::SysID, ver::Integer = 1) = execSys!(getfield(N, :jtsys), sysid, ver)

function execSys!(N::JTSystem, sysid::SysID, ver::Integer = 1)
    sidx = N.idtoindex[sysid]
    N.systems[sidx] isa Union{DiscreteSystem, VersionedSystem} || error("Can only execute (versioned) discrete-time systems.")

    Ad = N.Ad[sidx][ver]
    Rd = N.Rd[sidx][ver]

    N.m .= Ad * N.m

    # N.P = Ad * N.P * Ad' .+ Rd
    _mat_mul_add!(Ad, N.P, Ad', Rd, N.bufl)

    if N isa PeriodicJitterTimeSystem && N.periodicAnalysis
        # N.Atot = Ad * N.Atot
        mul!(N.bufl, Ad, N.Atot)
        N.Atot .= N.bufl

        # N.Rtot = Ad * N.Rtot * Ad' + Rd
        _mat_mul_add!(Ad, N.Rtot, Ad', Rd, N.bufl)
    end
    return nothing
end

"""
    `passTime!(N, T)`

Let time pass by T units, running all continuous-time systems and accumulating
cost.

__Arguments__:
N       The JitterTime system.
T       The time interval.
"""
function passTime!(N::JTSystem, time::Real)
    time >= 0 || error("Time must be positive!")
    time == 0 && return nothing # Special case when no time passed
    
    (Ad, Rd, Qd, Qconst) = calcC2D(N.reduced, time)

    #costm = N.m' * Qd * N.m
    costm = dot(N.m, Qd, N.m)
    #costP = tr(Qd * N.P) .+ Qconst
    mul!(N.bufl, Qd, N.P)
    costP = tr(N.bufl) .+ Qconst

    N.J += costm + costP
    N.dJdt = (costm + costP) / time

    N.m .= Ad * N.m

    # N.P = Ad * N.P * Ad' .+ Rd
    _mat_mul_add!(Ad, N.P, Ad', Rd, N.bufl)

    N.Tsim += time

    if N isa PeriodicJitterTimeSystem && N.periodicAnalysis
        # N.Atot = Ad * N.Atot
        mul!(N.bufl, Ad, N.Atot)
        N.Atot .= N.bufl

        # N.Rtot = Ad * N.Rtot * Ad' + Rd
        _mat_mul_add!(Ad, N.Rtot, Ad', Rd, N.bufl)

        N.dtot .= Ad * N.dtot
    end
    return nothing
end

"""
    `passTime!(N)`
Let time pass for a system with a fixed interval.

__Arguments__:
N       The fixed interval JitterTime system
"""
function passTime!(N::JTSystem)
    try N.FIQconst catch; error("JTSystem does not contain a FixedIntervalJitterTimeSystem!") end

    #costm = N.m' * Qd * N.m
    costm = dot(N.m, N.FIQd, N.m)
    #costP = tr(Qd * N.P) .+ Qconst
    mul!(N.bufl, N.FIQd, N.P)
    costP = tr(N.bufl) .+ N.FIQconst

    N.J += costm + costP
    N.dJdt = (costm + costP) / N.h

    N.m .= N.FIAd * N.m

    # N.P = Ad * N.P * Ad' .+ Rd
    _mat_mul_add!(N.FIAd, N.P, N.FIAd', N.FIRd, N.bufl)

    N.Tsim += N.h

    if N isa PeriodicJitterTimeSystem && N.periodicAnalysis
        # N.Atot = Ad * N.Atot
        mul!(N.bufl, N.FIAd, N.Atot)
        N.Atot .= N.bufl

        # N.Rtot = Ad * N.Rtot * Ad' + Rd
        _mat_mul_add!(N.FIAd, N.Rtot, N.FIAd', N.FIRd, N.bufl)

        N.dtot .= N.FIAd * N.dtot
    end
    return nothing
end

"""
    `passTimeUntil!(N, T)`

Let time pass until N.Tsim = T, running all continuous-time systems and
accumulating cost.

__Arguments__:
N       The JitterTime system.
T       The target time. Must be greater than or equal to N.Tsim.
"""
passTimeUntil!(N::JTSystem, endTime::Real) = passTime!(N, endTime - N.Tsim)

""" 
    `beginPeriodicAnalysis!(N)`

Start a periodic covariance analysis. Must be called after calcDynamics!.
"""
function beginPeriodicAnalysis!(N::PeriodicJitterTimeSystem)
    N.periodicAnalysis = true
    return nothing
end # function

""" 
    `endPeriodicAnalysis!(N)`

Start a periodic covariance analysis. Must be called after beginPeriodicAnalysis!
"""
function endPeriodicAnalysis!(N::PeriodicJitterTimeSystem)
    N.periodicAnalysis || error("Periodic analysis not started")
    N.Atot |> eigvalsCheck |> (x)->abs.(x) |> maximum < 1 || error("Periodic analysis failed due to unstable system mode(s)")

    N.Pper = dlyap(N.Atot, N.Rtot)
    N.mper = (I - N.Atot) \ N.dtot
    N.periodicAnalysis = false
    return nothing
end # function

"""
    `reset!(N)`

Reset dynamics of JitterTimeSystem

__Arguments__:
N       The JitterTime system to reset.
"""
function reset!(N::JTSystem)
    # Revert passTime!, execSys!
    N.J = 0.0
    N.dJdt = 0.0
    fill!(N.m, 0.0)
    fill!(N.P, 0.0)
    N.Tsim = 0.0
    return nothing
end
function reset!(N::PeriodicJitterTimeSystem{T}) where {T}
    N.periodicAnalysis = false
    N.Atot = Matrix{T}(I, size(N.P)...)
    fill!(N.Rtot, 0.0)
    fill!(N.dtot, 0.0)
    fill!(N.Pper, 0.0)
    fill!(N.mper, 0.0)
    reset!(N.jtsys)
end

end # module
