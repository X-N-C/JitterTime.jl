
"""
    `(Phi, R, Q, Qconst) = calcC2D(a, r, q, h)`

Calculate the discrete-time (ZOH) version of the continuous
system 
                `xdot = a*x + w`
where the incremental variance of w is r. The cost of the system
is 
                `J = integral_0^h (x'*q*x) dt`.
The resulting discrete-time system is
                `x(n+1) = Phi*x(n) + e(n)`
where the variance of e(n) is R, and the cost is
                `J = x'*Q*x + Qconst`.
"""
function calcC2D(N::ReducedSystem, dt::Real)
    copy!(N.bufAd, N.Id)
    copy!(N.bufRd, N.Zl)
    copy!(N.bufQd, N.Zl)
    (N.bufAd[N.indices, N.indices], 
         N.bufRd[N.indices, N.indices], 
         N.bufQd[N.indices, N.indices], 
         Qconst) = _calcC2D(N, dt)
    return (N.bufAd, N.bufRd, N.bufQd, Qconst)
end
function _calcC2D(N::ReducedSystem, dt::Real)
    if N.Mnorm * dt > 1.0 || N.maxAbsEig * dt > 4.0
        (Phi, R, Q, Qconst) = _calcC2D(N, dt/2)

        Qconst = 2*Qconst + tr(Q * R)
        Q .+= Phi' * Q * Phi
        R .+= Phi * R * Phi'
        Phi .= Phi*Phi
        return Phi, R, Q, Qconst
    end

    phi12 = jitterExp(N.M112, dt)
    phi22 = jitterExp(N.M122, dt)
    phi13 = jitterExp(N.M213, dt)
    phi23 = jitterExp(N.M223, dt)

    # Q = (Q+Q')*0.5
    Q = phi22' * phi12
    lmul!(0.5, Q)
    Q .+= Q'

    # R = (R+R')*0.5
    R = phi22 * phi23 #phi33'*phi23
    lmul!(0.5, R)
    R .+= R'

    Qconst = tr(phi22 * phi13 * N.Qc)
    return (phi22, R, Q, Qconst)
end

# Theta from: Computing the Matrix Exponential with an Optimized Taylor Polynomial Approximation
function jitterExp(M::Vector{Matrix{T}}, dt::Real) where {T}
    return jitterTaylor(M, dt, 18)
end

function jitterTaylor(M::Vector{Matrix{T}}, dt::Real, m::Int64) where {T}
    h = one(dt)
    p = zero(M[1])
    for k in 1:m+1
        #p .+= h .* M[k]
        axpy!(h, M[k], p)
        h *= dt
    end
    return p
end

######################
### TRACE OVERLOAD ###
######################

""" 
    `tr(A)`

Return the trace of an arbitrary matrix A.
The trace sums all the elements on the diagonal of A. 
"""
function tr(A::Matrix{T}) where {T}
    return sum((A)[idx, idx] for idx in 1:size(A, 1))
end # function

######################
## EIGVALS OVERLOAD ##
######################

""" Temporary eigenvalue check for matrix types who don't support LinearAlgebra.eigvals """
function eigvalsCheck(A::Matrix{T}) where {T}
    if isdefined(first(A), :value)
        return eigvals(getproperty.(A, :value))
    else
        return eigvals(A)
    end
end # function

######################
#### FREXP CHECK #####
######################

""" Temporary function for matrix types who don't support frexp """
function frexpCheck(A)
    if isdefined(first(A), :value)
        return frexp(getproperty.(A, :value))
    else
        return frexp(A)
    end
end # function

######################
# DISCRETE LYAPUNOV ##
######################

"""
    `dlyap(A, Q)`

Compute the solution `X` to the discrete Lyapunov equation
`AXA' - X + Q = 0`.
"""
function dlyap(A::Matrix{T}, Q) where {T}
    lhs = kron(A, conj(A))
    lhs = I - lhs
    x = lhs\reshape(Q, prod(size(Q)), 1)
    return reshape(x, size(Q))
end

######################
### AUX MATRIX MUL ###
######################
"""
    `_mat_mul!(A, B, C, Z)`
Compute the solution to A*B*C and store it in B
"""
function _mat_mul!(A::Matrix{T}, B::Matrix{T}, C::Matrix{T}, Z::Matrix{T}) where {T}
    mul!(Z, B, C)
    mul!(B, A, Z)
end
function _mat_mul!(A::Adjoint{T, Matrix{T}}, B::Matrix{T}, C::Matrix{T}, Z::Matrix{T}) where {T}
    mul!(Z, B, C)
    mul!(B, A, Z)
end
function _mat_mul!(A::Matrix{T}, B::Matrix{T}, C::Adjoint{T, Matrix{T}}, Z::Matrix{T}) where {T}
    mul!(Z, B, C)
    mul!(B, A, Z)
end

"""
    `_mat_mul_add!(A, B, C, D, Z)`
Compute the solution to A*B*C + D and store it in B
"""
function _mat_mul_add!(A::Matrix{T}, B::Matrix{T}, C::Matrix{T}, D::Matrix{T}, Z::Matrix{T}) where {T}
    _mat_mul!(A, B, C, Z)
    axpy!(1, D, B)
end
function _mat_mul_add!(A::Adjoint{T, Matrix{T}}, B::Matrix{T}, C::Matrix{T}, D::Matrix{T}, Z::Matrix{T}) where {T}
    _mat_mul!(A, B, C, Z)
    axpy!(1, D, B)
end
function _mat_mul_add!(A::Matrix{T}, B::Matrix{T}, C::Adjoint{T, Matrix{T}}, D::Matrix{T}, Z::Matrix{T}) where {T}
    _mat_mul!(A, B, C, Z)
    axpy!(1, D, B)
end
