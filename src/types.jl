export LinearSystem, 
       DiscreteSystem, 
       ContinuousSystem, 
       VersionedSystem, 
       JTSystem,
       JitterTimeSystem, 
       PeriodicJitterTimeSystem,
       FixedIntervalJitterTimeSystem

# Abstract types
""" 
    LinearSystem{T} 

Abstract supertype for all linear systems with elements of type `T`.
"""
abstract type LinearSystem{T} end

""" 
    JTSystem{T} 

Abstract supertype for all JitterTime systems with elements of type `T`.
"""
abstract type JTSystem{T} end

# Extracts the parametric type of an abstract type
_parametric_type(::LinearSystem{T}) where {T} = T
_parametric_type(::JTSystem{T}) where {T} = T

# Converts input to matrix form (and promotes it to type T)
_to_matrix(T, A::AbstractVector) = Matrix{T}(reshape(A, length(A), 1))
_to_matrix(T, A::AbstractMatrix) = (Ahat = similar(A, T); Ahat .= A) # Fallback
_to_matrix(T, A::Number)         = fill(T(A), 1, 1)
_to_matrix(T, A::Adjoint{AT, AM}) where {AT <: Number, AM <: AbstractMatrix} = _to_matrix(T, AM(A))

# Initial validation that the LinearSystem is feasible
function _init_validation(A, B, C, D, id, inputid)
    nx = size(A, 1)
    nu = size(B, 2)
    ny = size(C, 1)

    if size(A, 2) != nx && nx != 0 # FIXME size([],1)==0; size([],2)==1
        error("A must be square")
    elseif size(B, 1) != nx
        error("The number of rows of A ($(size(A,1))) and B ($(size(B,1))) are not equal")
    elseif size(C, 2) != nx # FIXME add nx != 0 check? This is wrong
        error("The number of columns of A ($(size(A,2))) and C ($(size(C,2))) are not equal")
    elseif nu != size(D, 2)
        error("The number of columns of B ($(size(B,2))) and D ($(size(D,2))) are not equal")
    elseif ny != size(D, 1)
        error("The number of rows of C ($(size(C,1))) and D ($(size(D,1))) are not equal")
    elseif id <= 0
        error("ID must be positive")
    elseif id in inputid
        error("System ($id) can't be connected to itself")
    end

    return nx, nu, ny
end # function

###############
### Structs ###
###############
"""
    ContinuousSystem(sysid, A, B, C, D, inputid [, Rc, Qc])

Create a continuous-time linear system
```
u(t) --> sys --> y(t)
         x(t)
```
The state initially has mean value `E(x) = 0` and covariance `V(x) = 0`.

## Arguments:
* `sysid`:          A unique positive `Int64` identifier for this system (pick any). Used when referred to from other systems. 
* `[A, B, C, D]`:   A strictly proper, delay-free continuous-time LTI system in state-space form. 
* `inputid`:        A vector of input system identifiers. 
                    The outputs of the corresponding systems will be used as inputs to this system.
                    The number of inputs in this system must equal the total number of outputs in the input systems.
                    An empty vector (or zero) indicates that the system inputs are unconnected.
* `Rc`:             (OPTIONAL) The continuous state or input noise intensity matrix. 
                    The noise vector is assumed to have the same size as x (for state-space systems).
                    Default assumed zero.
* `Qc`:             (OPTIONAL) The continuous cost function is `E(Integral [x(t); u(t)]' * Qc * [x(t); u(t)] dt)` (for state-space systems).
                    Default assumed zero.
"""
struct ContinuousSystem{T} <: LinearSystem{T}
    A::Matrix{T}
    B::Matrix{T}
    C::Matrix{T}
    Rc::Matrix{T} # Noise intensity matrix (continuous)
    Qc::Matrix{T} # Weight matrix (continuous)

    inputid::Vector{Int64}
    n::Integer # Number of states
    r::Integer # Number of inputs
    p::Integer # Number of outputs
    id::Int64
    resetDynamics::Bool
end

# Constructors
function ContinuousSystem(id::Int64, A, B, C, D, inputid::Vector{Int64}, 
                          Rc = zeros(size(A, 2), size(A, 2)), 
                          Qc = zeros(size(A, 2)+size(B, 2), size(A, 2)+size(B, 2)))
    n, r, p = _init_validation(A, B, C, D, id, inputid)
    (size(Rc,1),size(Rc,2)) == (n,n) || error("Rc ($(size(Rc,1))*$(size(Rc,2))) should be a n*n matrix; n = #states (= $n)")
    (size(Qc,1),size(Qc,2)) == (n+r,n+r) || error("Qc ($(size(Qc,1))*$(size(Qc,2))) should be an (n+r)*(n+r) matrix; n = #states (= $n), r = #inputs (= $r)")

    T = promote_type(eltype(A), eltype(B), eltype(C), eltype(D), eltype(Rc), eltype(Qc))
    return ContinuousSystem{T}(_to_matrix(T, A), _to_matrix(T, B), _to_matrix(T, C), _to_matrix(T, Rc), _to_matrix(T, Qc), inputid, n, r, p, id, false)
end
ContinuousSystem(id::Int64, A, B, C, D, inputid::Int64, 
                 Rc = zeros(size(A, 2), size(A, 2)), 
                 Qc = zeros(size(A, 2)+size(B, 2), size(A, 2)+size(B, 2))) = 
    ContinuousSystem(id, A, B, C, D, Int64[inputid], Rc, Qc)

"""
    DiscreteSystem(sysid, Phi, Gam, C, D, inputid [, Rd, Qc])

Create a discrete-time linear system
```
u(k) --> sys --> y(k)
         x(k)
```
The state/output initially has mean value `E([x; y]) = 0` and covariance `V([x; y]) = 0`.

## Arguments:
* `sysid`:              A unique positive `Int64` identifier for this system (pick any). Used when referred to from other systems. 
* `[Phi, Gam, C, D]`:   A discrete-time LTI system in state-space form (if `Phi`, `Gam`, and `C` are missing, interpreted as a static gain).
* `inputid`:            A vector of input system identifiers. 
                        The outputs of the corresponding systems will be used as inputs to this system.
                        The number of inputs in this system must equal the total number of outputs in the input systems.
                        An empty vector (or zero) indicates that the system inputs are unconnected.
* `Rd`:                 (OPTIONAL) The discrete state/output/input noise covariance matrix. 
                        The noise vector is assumed to have the same size as `[x; y]` (for state-space systems).
                        Noise is added each time the system is executed.
                        Default assumed zero.
* `Qc`:                 (OPTIONAL) The continuous cost function is `E(Integral [x(t); u(t)]' * Qc * [x(t); u(t)] dt)` (for state-space systems).
                        _Note_ that both `x(t)` and `y(t)` are held constant between system executions.
                        Default assumed zero.
"""
struct DiscreteSystem{T} <: LinearSystem{T}
    A::Matrix{T}
    B::Matrix{T}
    C::Matrix{T}
    R::Matrix{T}   # Noise
    Qc::Matrix{T}  # Cost

    inputid::Vector{Int64}
    n::Integer # Number of states
    r::Integer # Number of inputs
    p::Integer # Number of outputs
    id::Int64
end

# Constructor
function DiscreteSystem(id::Int64, Phi, Gam, C, D, inputid::Vector{Int64},
                        R = zeros(size(Phi, 2)+size(C, 1), size(Phi, 2)+size(C, 1)),
                        Qc = zeros(size(Phi, 2)+size(C, 1), size(Phi, 2)+size(C, 1)))
    n, r, p = _init_validation(Phi, Gam, C, D, id, inputid)
    (size(R,1),size(R,2)) == (n+p,n+p) || error("R ($(size(R,1))*$(size(R,2))) should be a (n+p)*(n+p) matrix; n = #states (= $n), p = #outputs (= $p).")
    (size(Qc,1),size(Qc,2)) == (n+p,n+p) || error("Qc ($(size(Qc,1))*$(size(Qc,2))) should be a (n+p)*(n+p) matrix; n = #states (= $n), p = #outputs (= $p).")

    Aarray = [Phi zeros(n, p); C zeros(p, p)]
    Barray = [Gam; D]

    T = promote_type(eltype(Phi), eltype(Gam), eltype(C), eltype(D), eltype(R), eltype(Qc))
    return DiscreteSystem{T}(_to_matrix(T, Aarray), _to_matrix(T, Barray), _to_matrix(T, [zeros(p, n) I]), _to_matrix(T, R), _to_matrix(T, Qc), inputid, n+p, r, p, id)
end
DiscreteSystem(id::Int64, Phi, Gam, C, D, inputid::Int64,
               R = zeros(size(Phi, 2)+size(C, 1), size(Phi, 2)+size(C, 1)),
               Qc = zeros(size(Phi, 2)+size(C, 1), size(Phi, 2)+size(C, 1))) = 
    DiscreteSystem(id, Phi, Gam, C, D, Int64[inputid], R, Qc)

function DiscreteSystem(id::Int64, D, inputid::Vector{Int64},
                        R = zeros(size(D, 1), size(D, 1)),
                        Qc = zeros(size(D, 1), size(D, 1)))
    n, r, p = 0, size(D, 2), size(D, 1)
    (size(R,1),size(R,2)) == (p,p) || error("R ($(size(R,1))*$(size(R,2))) should be a (n+p)*(n+p) matrix; n = #states (= $n), p = #outputs (= $p).")
    (size(Qc,1),size(Qc,2)) == (p,p) || error("Qc ($(size(Qc,1))*$(size(Qc,2))) should be a (n+p)*(n+p) matrix; n = #states (= $n), p = #outputs (= $p).")

    Aarray = zeros(p, p)
    Barray = D

    T = promote_type(eltype(D), eltype(R), eltype(Qc))
    return DiscreteSystem{T}(_to_matrix(T, Aarray), _to_matrix(T, Barray), Matrix{T}(I, p, p), _to_matrix(T, R), _to_matrix(T, Qc), inputid, n+p, r, p, id)
end
DiscreteSystem(id::Int64, D, inputid::Int64,
               R = zeros(size(D, 1), size(D, 1)),
               Qc = zeros(size(D, 1), size(D, 1))) = 
    DiscreteSystem(id, D, Int64[inputid], R, Qc)

"""
    VersionedSystem(versions::Vector{DiscreteSystem})

Create a versioned discrete-time linear system; the version identifier is given by its index in the input vector.
"""
struct VersionedSystem{T} <: LinearSystem{T}
    versions::Vector{DiscreteSystem}

    # Constructor
    function VersionedSystem(vers::Vector{DiscreteSystem{T}}) where {T}
        S = promote_type(_parametric_type.(vers)...)

        length(vers) > 1 || error("More than one DiscreteSystem necessary to use VersionedSystem")
        all(x -> x.id == vers[1].id, vers) || error("Versions do not share ID")
        all(x -> x.n == vers[1].n, vers) || error("Versions do not share state-space dimensions")
        all(x -> x.r == vers[1].r, vers) || error("Versions do not share input dimensions")
        all(x -> x.p == vers[1].p, vers) || error("Versions do not share output dimensions")
        all(x -> _to_matrix(S, x.Qc) == _to_matrix(S, vers[1].Qc), vers) || error("Qc is defined in continuous time and cannot have multiple versions")

        new{S}(vers)
    end
end

#= Override get function for structs of type VersionedSystem =#
function Base.getproperty(v::VersionedSystem{T}, s::Symbol) where {T}
    if s === :C || s === :Qc
        return _to_matrix(T, getproperty(v.versions[1], s))
    elseif s === :n || s === :r || s === :p || s === :id
        return getproperty(v.versions[1], s)
    end # if
    return getfield(v, s)
end # function

#= Internal ReducedSystem struct =#
mutable struct ReducedSystem{T}
    n::Integer
    Ac::Matrix{T}
    Rc::Matrix{T}
    Qc::Matrix{T}
    indices::Vector{Integer}

    # Submatrices used in calcC2D
    M112::Vector{Matrix{T}}
    M122::Vector{Matrix{T}}
    M213::Vector{Matrix{T}}
    M223::Vector{Matrix{T}}

    # Sanity checks
    maxAbsEig::Real
    Mnorm::T

    # Matrix buffers to reduce allocation time
    bufAd::Matrix{T}
    bufRd::Matrix{T}
    bufQd::Matrix{T}
    bufl::Matrix{T} # Large buffer

    # Commonly used matrices
    Id::Matrix{T}   # Large identity matrix
    Zl::Matrix{T}   # Large zero matrix

    function ReducedSystem(Ac::Matrix{T},
                           Rc::Matrix{T},
                           Qc::Matrix{T}) where {T}
        n = size(Ac, 1)
        # Find all closed-loop states which are continuous
        mask = iszero.(Ac) .& iszero.(Rc) .& iszero.(Qc)
        indices = filter(!iszero, (vec(any(!, mask, dims=1)) .| vec(any(!, mask, dims=2))) .* (1:n))

        # Reduce the continuous matrices
        reducedAc = Ac[indices,indices]
        reducedRc = Rc[indices,indices]
        reducedQc = Qc[indices,indices]

        # Create buffers for reduced discretised continuous system
        bufAd = Matrix{T}(undef, n, n)
        bufRd = Matrix{T}(undef, n, n)
        bufQd = Matrix{T}(undef, n, n)

        M122 = [reducedAc^k for k in 0:1]
        M112 = [zero(reducedQc), reducedQc]
        M223 = [zero(reducedRc), reducedRc']
        M213 = [zero(reducedRc), zero(reducedRc)]
        for k in 2:18
            c = 1.0/k

            push!(M213, c*(-reducedAc * last(M213) + last(M223))) # + zero matrix
            push!(M223, c*(-reducedAc * last(M223) + reducedRc'* last(M122)'))

            push!(M112, c*(-reducedAc'* last(M112) + reducedQc * last(M122)))
            push!(M122, c*( reducedAc * last(M122)))
        end

        # Norm used for scaling and squaring
        Meig = maximum(abs.(eigvalsCheck(reducedAc))) # TODO Remove eigvalsCheck and handle some other way...
        l = size(reducedAc, 1)
        M = zeros(T, 2 .* size(reducedAc))
        M[1:l, 1:l] = -reducedAc'
        M[1:l, l+1:2*l] = reducedQc
        M[l+1:2*l, l+1:2*l] = reducedAc
        M1norm = norm(M, 1)
        
        Id = Matrix{T}(I, size(reducedAc))
        M = zeros(T, 3 .* size(reducedAc))
        M[1:l, 1:l] = M[l+1:2*l, l+1:2*l] = -reducedAc
        M[2*l+1:3*l, 2*l+1:3*l] = reducedAc'
        M[1:l, l+1:2*l] = Id
        M[l+1:2*l, 2*l+1:3*l] = reducedRc'
        M2norm = norm(M, 1)
        return new{T}(n, reducedAc, reducedRc, reducedQc, indices, 
                      M112, M122, M213, M223, 
                      Meig, max(M1norm, M2norm),
                      bufAd, bufRd, bufQd, zeros(T, n, n),
                      Matrix{T}(I, n, n), zeros(T, n, n))
    end
end


"""
    JitterTimeSystem(systems::Vector{LinearSystem})

Initialize a new JitterTime system.
"""
mutable struct JitterTimeSystem{T} <: JTSystem{T}
    systems::Vector{LinearSystem}
    idtoindex::Dict{Int64, Int64}
    J::T
    dJdt::T
    Tsim::T

    Ad::Vector{Vector{Matrix{T}}} # TODO: Change this complex structure to individual structs instead
    Rd::Vector{Vector{Matrix{T}}} # TODO: Change this complex structure to individual structs instead
    m::Vector{T}
    P::Matrix{T}
    
    reduced::ReducedSystem{T}

    # Implicit Constructor
    function JitterTimeSystem(T, systems::Vector{S}) where {S <: LinearSystem}
        totstates::Int64 = sum(sys.n for sys in systems)

        Ad::Vector{Vector{Matrix{T}}}  = [Matrix{T}[] for _ in systems]
        Rd::Vector{Vector{Matrix{T}}}  = [Matrix{T}[] for _ in systems]
        P                              = zeros(T, totstates, totstates)
        m                              = zeros(T, totstates)

        return new{T}(systems, Dict{Integer, Integer}(), T(0), T(0), T(0), Ad, Rd, m, P)
    end # Constructor
end
JitterTimeSystem(systems::Vector{S}) where {S <: LinearSystem} = JitterTimeSystem(promote_type(_parametric_type.(systems)...), systems)

#= Override get function for structs of type JitterTimeSystem =#
function Base.getproperty(sys::JitterTimeSystem, s::Symbol)
    if hasproperty(sys, s)
        return getfield(sys, s)
    end # if
    return getproperty(sys.reduced, s)
end # function
#= Override set function for structs of type JitterTimeSystem =#
function Base.setproperty!(sys::JitterTimeSystem, s::Symbol, x)
    if hasproperty(sys, s)
        return setfield!(sys, s, x)
    end # if
    return setproperty!(sys.reduced, s, x)
end # function

"""
    PeriodicJitterTimeSystem(N::JTSystem)

Initialize a new Periodic JitterTime system.
"""
mutable struct PeriodicJitterTimeSystem{T} <: JTSystem{T}
    jtsys::JTSystem{T}
    periodicAnalysis::Bool

    Atot::Matrix{T}
    Rtot::Matrix{T}
    dtot::Matrix{T} # TODO: What is this?
    Pper::Matrix{T}
    mper::Matrix{T} # TODO: What is this?

    # Implicit Constructor
    function PeriodicJitterTimeSystem(jtsys::S) where {S <: JTSystem}
        !(jtsys isa PeriodicJitterTimeSystem) || error("Unecessary wrapping of PeriodicJitterTimeSystem")

        T = _parametric_type(jtsys)

        Atot = Matrix{T}(I, size(jtsys.P))
        Rtot = zeros(T, size(jtsys.P))
        dtot = zeros(T, size(jtsys.P, 1), 1)
        Pper = zero(Rtot)
        mper = zero(dtot)

        new{T}(jtsys, false, Atot, Rtot, dtot, Pper, mper)
    end
end

#= Override get function for structs of type PeriodicJitterTimeSystem =#
function Base.getproperty(sys::PeriodicJitterTimeSystem, s::Symbol)
    if hasproperty(sys, s)
        return getfield(sys, s)
    end # if
    return getproperty(sys.jtsys, s)
end # function
#= Override set function for structs of type PeriodicJitterTimeSystem =#
function Base.setproperty!(sys::PeriodicJitterTimeSystem, s::Symbol, x)
    if hasproperty(sys, s)
        return setfield!(sys, s, x)
    end # if
    return setproperty!(sys.jtsys, s, x)
end # function

"""
    FixedIntervalJitterTimeSystem(N::JTSystem, h::Real)

Initialize a new JitterTime system with fixed activation interval h.
"""
mutable struct FixedIntervalJitterTimeSystem{T} <: JTSystem{T}
    jtsys::JTSystem{T} # TODO: Should this be of type JTSystem{T}?
    h::Real

    FIAd::Matrix{T}
    FIRd::Matrix{T}
    FIQd::Matrix{T}
    FIQconst::T

    function FixedIntervalJitterTimeSystem(jtsys::S, h::Real) where {S <: JTSystem}
        h > 0 || error("h must be positive!")
        T = _parametric_type(jtsys)
        tmp_sys = deepcopy(jtsys) # create deepcopy to calculate dynamic on
        calcDynamics!(tmp_sys)
        (Ad, Rd, Qd, Qconst) = calcC2D(tmp_sys.reduced, h) # TODO: This is bad practice... We assume the jtsys has had its dynamics calculated before hand (even though we never do it elsewhere)

        new{T}(jtsys, h, Ad, Rd, Qd, Qconst)
    end
        
end

#= Override get function for structs of type FixedIntervalJitterTimeSystem =#
function Base.getproperty(sys::FixedIntervalJitterTimeSystem, s::Symbol)
    if hasproperty(sys, s)
        return getfield(sys, s)
    end # if
    return getproperty(sys.jtsys, s)
end # function
#= Override set function for structs of type FixedIntervalJitterTimeSystem =#
function Base.setproperty!(sys::FixedIntervalJitterTimeSystem, s::Symbol, x)
    if hasproperty(sys, s)
        return setfield!(sys, s, x)
    end # if
    return setproperty!(sys.jtsys, s, x)
end # function
