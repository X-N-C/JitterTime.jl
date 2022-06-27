export execSys!,
       passTime!,
       passTimeUntil!,
       beginPeriodicAnalysis!,
       endPeriodicAnalysis!


"""
    `execSys!(N, sysid [, ver])`

Execute the discrete-time system 'sysid'

__Arguments__:
N           The JitterTime system.
sysid       ID of the discrete-time system to be updated.

__Optional Arguments__:
ver         What version to execute (default = 1)
"""
execSys!(N::FixedIntervalJitterTimeSystem, sysid::Integer, ver::Integer = 1) = execSys!(getfield(N, :jtsys), sysid, ver)

function execSys!(N::JTSystem, sysid::Integer, ver::Integer = 1)
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

