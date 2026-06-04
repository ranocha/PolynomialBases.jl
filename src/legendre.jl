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
