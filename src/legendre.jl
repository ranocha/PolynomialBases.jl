"""
    legendre(x, p::Integer)

Evaluate the Legendre polynomial of degree `p` at `x` using the three term
recursion [Kopriva, Implementing Spectral Methods for PDEs, Algorithm 20].
"""
function legendre(x, p::Integer)
    a::typeof(x) = one(x)
    b::typeof(x) = x

    if p <= 0
        return a
    elseif p == 1
        return b
    end

    for j in 2:p
        a, b = b, ( (2j-1)*x*b - (j-1)*a ) / j
    end

    b
end

"""
    legendre_and_derivative(x, p::Integer)

Evaluate the Legendre polynomial of degree `p` and its derivative at `x` using
the three term recursion [Kopriva, Implementing Spectral Methods for PDEs,
Algorithm 22].
"""
function legendre_and_derivative(x, p::Integer)
    # coefficients for the polynomial...
    a::typeof(x) = one(x)
    b::typeof(x) = x
    # ... and for the derivative
    aa::typeof(x) = zero(x)
    bb::typeof(x) = one(x)

    if p <= 0
        return a, aa
    elseif p == 1
        return b, bb
    end

    for j in 2:p
        a, b = b, ( (2j-1)*x*b - (j-1)*a ) / j
        aa, bb = bb, aa + (2j-1)*a
    end

    b, bb
end


"""
    gauss_legendre_nodes_and_weights(p, T=Float64::Type, tol=4*eps(T), maxit=100)

Compute the Gauss-Legendre nodes and weights for polynomials of degree `p`
using the scalar type `T`, tolerance `tol` and maximal number of Newton iterations
`maxit` [Kopriva, Implementing Spectral Methods for PDEs, Algorithm 23].
"""
function gauss_legendre_nodes_and_weights(p, T=Float64::Type, tol=4*eps(T), maxit=100)
    if p <= 0
        return T[0], T[2]
    elseif p == 1
        x0 = -1 / sqrt(T(3))
        return [x0, -x0], T[1, 1]
    end

    nodes = Vector{T}(undef, p+1)
    weights = Vector{T}(undef, p+1)
    for j in 0:((p+1)÷2-1)
        x = -cospi(T(2j+1)/(2p+2))
        for k in 1:maxit
            pol, der = legendre_and_derivative(x, p+1)
            Δ = -pol / der
            x = x + Δ
            abs(Δ) <= tol*abs(x) && break
        end
        pol, der = legendre_and_derivative(x, p+1)
        nodes[j+1] = x
        nodes[p-j+1] = -x
        weights[j+1] = weights[p-j+1] = 2 / ( (1-x^2)*der^2 )
    end

    if mod(p,2) == 0
        pol, der = legendre_and_derivative(zero(T), p+1)
        nodes[p÷2+1]   = 0
        weights[p÷2+1] = 2 / der^2
    end

    return nodes, weights
end


# helper function [Kopriva, Implementing Spectral Methods for PDEs, Algorithm 24]
function q_and_L_evaluation(x, p::Integer)
    a = one(x)
    b = x
    aa = zero(x)
    bb = one(x)

    for j in 2:p
        a, b = b, ( (2j-1)*x*b - (j-1)*a ) / j
        aa, bb = bb, aa + (2j-1)*a
    end
    pol = ( (2p+1)*x*b - p*a ) / (p+1)
    der = aa + (2p+1)*b

    pol-a, der-aa, b
end

"""
    lobatto_legendre_nodes_and_weights(p, T=Float64::Type, tol=4*eps(T), maxit=100)

Compute the Lobatto-Legendre nodes and weights for polynomials of degree `p`
using the scalar type `T`, tolerance `tol` and maximal number of Newton iterations
`maxit` [Kopriva, Implementing Spectral Methods for PDEs, Algorithm 25].
"""
function lobatto_legendre_nodes_and_weights(p, T=Float64::Type, tol=4*eps(T), maxit=100)
    if p <= 0
        return T[0], T[2]
    elseif p == 1
        return T[-1, 1], T[1, 1]
    end

    nodes = Vector{T}(undef, p+1)
    weights = Vector{T}(undef, p+1)

    nodes[1] = -1
    nodes[end] = 1
    weights[1] = weights[end] = T(2) / (p*(p+1))

    for j in 1:((p+1)÷2-1)
        x = -cos( T(j+1//4)*π/p - 3/(8p*T(j+1//4)*π) )
        for k in 1:maxit
            pol, der, leg = q_and_L_evaluation(x, p)
            Δ = -pol / der
            x = x + Δ
            abs(Δ) <= tol*abs(x) && break
        end
        pol, der, leg = q_and_L_evaluation(x, p)
        nodes[j+1] = x
        nodes[p-j+1] = -x
        weights[j+1] = weights[p-j+1] = 2 / ( p*(p+1)*leg^2 )
    end

    if mod(p,2) == 0
        pol, der, leg = q_and_L_evaluation(zero(T), p)
        nodes[p÷2+1]   = 0
        weights[p÷2+1] = 2 / ( p*(p+1)*leg^2 )
    end

    return nodes, weights
end



"""
    legendre_vandermonde(nodes)

Computes the Vandermonde matrix with respect to the Legendre polynomials and
the nodal basis on `nodes`.
The Vandermonde matrix `V` is the transformation matrix from the modal Legendre
basis to the nodal Lagrange basis associated with `nodes`.
"""
function legendre_vandermonde(nodes::AbstractVector)
    T = eltype(nodes)
    pp1 = length(nodes)
    V = Array{T}(undef, pp1, pp1)
    for j in 1:pp1, (i,x) in enumerate(nodes)
      V[i, j] = legendre(x, j-1)
    end
    V
end

legendre_vandermonde(basis::NodalBasis{Line}) = legendre_vandermonde(grid(basis))


"""
    legendre_M(p, T=Float64)

Computes the diagonal mass matrix in the modal Legendre basis up to degree `p`
using the scalar type `T`.
"""
function legendre_M(p, T=Float64)
  Diagonal( T[2//(2n+1) for n in 0:p] )
end


"""
    legendre_D(p, T=Float64)

Computes the derivative matrix in the modal Legendre basis up to degree `p`
using the scalar type `T`.
"""
function legendre_D(p, T=Float64)
  D = fill(zero(T), p+1, p+1)
  if p >= 1
    D[1, 2] = 1
  end
  if p >= 2
    D[2, 3] = 3
  end
  for col in 4:p+1
    D[:, col] = D[:, col-2]
    D[col-1, col] = 2*col-3
  end
  D
end
