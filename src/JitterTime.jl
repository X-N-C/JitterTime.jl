module JitterTime

using LinearAlgebra

include("types.jl")             # Containing all datastructures necessary for JitterTime
include("util.jl")              # Containing auxiliary functions
include("execution.jl")         # Containing all code for Executing the JitterTimeSystem
include("printing.jl")          # Containing all Pretty Printing code

export calcDynamics!, 
       reset!

"""
    calcDynamics!(N::JTSystem)

Calculate the total system dynamics for the JitterTime system `N`. 

Assuming the state of the total system (including all subsystems) to be x.

The continuous time evolution follows the equation
```
dx(t) = N.Ac x(t)dt + dvc(t)
```
where `vc` is continuous white noise with intensity `N.Rc`.

The discrete-event evolution (when executing system `S`) follows the equation
```
x+(t_k) = S.Ad x(t_k) + vd(t_k)
```
where vd is white Gaussian noise with variance `S.Rd`.
"""
function calcDynamics!(N::JTSystem{T}) where {T}
    # Step 1: Count states and inputs and build a system-to-state index mapping
    indices             = Vector{UnitRange{Int64}}(undef, length(N.systems))
    totstates::Int64    = 0
    n::Int64            = 0

    empty!(N.idtoindex)
    merge!(N.idtoindex,
           Dict(N.systems[x].id::Int64 => x for x in eachindex(N.systems)))
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
    if !_connections_correct(N) 
        error("Wrong number of inputs")
    end # if 

    # Initialice Ac, Rc, Qc
    Ac = zero(N.P)
    Rc = zero(N.P)
    Qc = zero(N.P)

    # Step 3: Formulate the total continuous-time dynamics, continuous noise, discrete-time dynamics (all versions), discrete noise and continuous cost
    # TODO Add reset dynamics?
    for (k, s) in enumerate(N.systems)
        xtou = zeros(T, s.n, totstates)
        xtou[:, indices[k]] = Matrix{T}(I, s.n, s.n) # = I

        _formulate_dynamics!(s, N, Ac, Rc, Qc, xtou, indices, totstates, k)
    end # for

    # Step 4: Create reduced system to improve performance
    N.reduced = ReducedSystem(Ac, Rc, Qc)

    return nothing
end

"""
    reset!(N::JTSystem)

Reset the dynamics of the JitterTime system.
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

###########################
### Auxiliary functions ###
###########################

#= Check that the number of outputs and inputs match for all system connections =#
function _connections_correct(N::JTSystem)
    for s in N.systems
        if !_connections_correct(s, N)
            return false
        end # if 
    end # for
    return true
end # function
function _connections_correct(s::VersionedSystem, N::JTSystem)
    for v in s.versions
        if sum(N.systems[N.idtoindex[x]].p for x in v.inputid) != v.r
            return false
        end # if
    end # for
    return true
end # function
function _connections_correct(s::LinearSystem, N::JTSystem)
    return sum(N.systems[N.idtoindex[x]].p for x in s.inputid) == s.r
end # function

#= Formulate dynamics depending on LinearSystem type =#
function _formulate_dynamics!(s::ContinuousSystem,
                              N::JTSystem{T},
                              Ac::Matrix{T},
                              Rc::Matrix{T},
                              Qc::Matrix{T},
                              xtou::Matrix,
                              indices::Vector{UnitRange{S}},
                              totstates::S,
                              k::S) where {T, S <: Integer}
    # Construct Ac matrix
    Ac[indices[k], indices[k]] = s.A
    bix::Int64 = 1

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

    # Construct Qc matrix
    Qc .+= xtou'*s.Qc*xtou # Combined cost of state, outputs and inputs
end # function
function _formulate_dynamics!(s::DiscreteSystem,
                              N::JTSystem{T},
                              Ac::Matrix{T},
                              Rc::Matrix{T},
                              Qc::Matrix{T},
                              xtou::Matrix,
                              indices::Vector{UnitRange{S}},
                              totstates::S,
                              k::S) where {T, S <: Integer}
    
    sidx::Int64 = N.idtoindex[s.id]
    push!(N.Ad[sidx], Matrix{T}(I, totstates, totstates))
    N.Ad[sidx][1][indices[k], indices[k]] = s.A
    push!(N.Rd[sidx], zeros(T, totstates, totstates))
    N.Rd[sidx][1][indices[k], indices[k]] = s.R
    
    bix::Int64 = 1
    for inputindex in map(x -> N.idtoindex[x], s.inputid)
        input = N.systems[inputindex]
        N.Ad[sidx][1][indices[k], indices[inputindex]] = s.B[:, bix:bix+input.p-1]*input.C
        bix += input.p
    end # for

    # Construct Qc matrix
    Qc .+= xtou'*s.Qc*xtou # Combined cost of state, outputs and inputs
end # function
function _formulate_dynamics!(s::VersionedSystem,
                              N::JTSystem{T},
                              Ac::Matrix{T},
                              Rc::Matrix{T},
                              Qc::Matrix{T},
                              xtou::Matrix,
                              indices::Vector{UnitRange{S}},
                              totstates::S,
                              k::S) where {T, S <: Integer}
    
    for (v, sver) in enumerate(s.versions)
        sidx::Int64 = N.idtoindex[s.id]
        push!(N.Ad[sidx], Matrix{T}(I, totstates, totstates)) # = I
        N.Ad[sidx][v][indices[k], indices[k]] = sver.A
        push!(N.Rd[sidx], zeros(T, totstates, totstates))
        N.Rd[sidx][v][indices[k], indices[k]] = sver.R
        
        bix::Int64 = 1
        for inputindex in map(x -> N.idtoindex[x], sver.inputid)
            input = N.systems[inputindex]
            N.Ad[sidx][v][indices[k], indices[inputindex]] = sver.B[:, bix:bix+input.p-1]*input.C
            bix += input.p
        end # for
    end # for

    # Construct Qc matrix
    Qc .+= xtou'*s.Qc*xtou # Combined cost of state, outputs and inputs
end # function

end # module
