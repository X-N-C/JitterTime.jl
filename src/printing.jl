# Structs
function Base.show(io::IO, sys::DiscreteSystem{T}) where {T}
    n = sys.n - sys.p # Number of states in phi
    println(io, "$(n == 0 ? "Static " : "")DiscreteSystem{$T} (ID = $(sys.id))\n{")
    if n != 0
        println(io, "\tA:  $(sys.A[1:n, 1:n])") # TODO Which matrices should be printed?
        println(io, "\tB:  $(sys.B[1:n, :])")
        println(io, "\tC:  $(sys.A[n+1:end, 1:n])")
        println(io, "\tD:  $(sys.B[n+1:end, :])")
    else
        println(io, "\tD:  $(sys.B[n+1:end, :])")
    end
    println(io, "\tR:  $(sys.R)")
    println(io, "\tQc: $(sys.Qc)")
    print(io, "} Input from LinearSystem$((length(sys.inputid) == 1) ? ": $(sys.inputid[1])" : "s $(sys.inputid)")")
end # function

function Base.show(io::IO, sys::ContinuousSystem{T}) where {T}
    println(io, "ContinuousSystem{$T} (ID = $(sys.id))\n{")
    println(io, "\tA:  $(sys.A)")
    println(io, "\tB:  $(sys.B)")
    println(io, "\tC:  $(sys.C)")
    println(io, "\tRc: $(sys.Rc)")
    println(io, "\tQc: $(sys.Qc)")
    print(io, "} Input from LinearSystem$((length(sys.inputid) == 1) ? ": $(sys.inputid[1])" : "s: $(sys.inputid)")")
end # function

function Base.show(io::IO, sys::VersionedSystem{T}) where {T}
    println(io, "VersionedSystem{$T} (ID = $(sys.id)) containing $(length(sys.versions)) DiscreteSystems")
end # function

function Base.show(io::IO, sys::JitterTimeSystem{T}) where {T}
    print(io, "JitterTimeSystem{$T} containing LinearSystem$((length(sys.systems) == 1) ? ": $(sys.systems[1].id)" : "s: $(map(x -> x.id, sys.systems))")")
end # function

function Base.show(io::IO, sys::PeriodicJitterTimeSystem{T}) where {T}
    print(io, "PeriodicJitterTimeSystem{$T} containing LinearSystem$((length(sys.systems) == 1) ? ": $(sys.systems[1].id)" : "s: $(map(x -> x.id, sys.systems))")")
end # function

function Base.show(io::IO, sys::FixedIntervalJitterTimeSystem{T}) where {T}
    print(io, "FixedIntervalJitterTimeSystem{$T} (interval = $(sys.h)) containing LinearSystem$((length(sys.systems) == 1) ? ": $(sys.systems[1].id)" : "s: $(map(x -> x.id, sys.systems))")")
end # function
