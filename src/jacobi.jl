"""
    jacobi(x, p::Integer, α, β)

Evaluate the Legendre polynomial with parameters `α`, `β` of degree `p` at `x`
using the three term recursion [Karniadakis and Sherwin, Spectral/hp Element
Methods for CFD, Appendix A].
"""
function jacobi(x, p::Integer, α, β)
    T = typeof( (2+α+β)*x / 2 )
    a = one(T)
    b = ((2+α+β)*x + α - β) / 2

    if p <= 0
        return a
    elseif p == 1
        return b
    end

    for n in 2:p
        a1 = 2n*(n+α+β)*(2n-2+α+β)
        a2 = (2n-1+α+β)*(α+β)*(α-β)
        a3 = (2n-2+α+β)*(2n-1+α+β)*(2n+α+β)
        a4 = 2*(n-1+α)*(n-1+β)*(2n+α+β)
        a, b = b, ( (a2+a3*x)*b - a4*a ) / a1
    end

    b
end

"""
    jacobi_and_derivative(x, p::Integer, α, β)

Evaluate the Jacobi polynomial with parameters `α`, `β` of degree `p` and its
derivative at `x` using the three term recursion [Karniadakis and Sherwin,
Spectral/hp Element Methods for CFD, Appendix A].
"""
function jacobi_and_derivative(x, p::Integer, α, β)
    T = typeof( (2+α+β)*x / 2 )
    # Coefficients for the polynomial...
    a = one(T)
    b = ((2+α+β)*x + α - β) / 2
    # ... and for the derivative. The derivative is computed via eq. (A.1.8)
    # of Karniadakis and Sherwin as Jacobi polynomial with parameters α+1, β+1.
    aa = one(T)
    bb = ((4+α+β)*x + α - β) / 2

    if p <= 0
        return a, zero(T)
    elseif p == 1
        return b, (2+α+β)/2
    end

    for n in 2:p
        a1 = 2n*(n+α+β)*(2n-2+α+β)
        a2 = (2n-1+α+β)*(α+β)*(α-β)
        a3 = (2n-2+α+β)*(2n-1+α+β)*(2n+α+β)
        a4 = 2*(n-1+α)*(n-1+β)*(2n+α+β)
        a, b = b, ( (a2+a3*x)*b - a4*a ) / a1

        b1 = 2n*(n+2+α+β)*(2n+α+β)
        b2 = (2n+1+α+β)*(α+β+2)*(α-β)
        b3 = (2n+α+β)*(2n+1+α+β)*(2n+2+α+β)
        b4 = 2*(n+α)*(n+β)*(2n+2+α+β)
        aa, bb = bb, ( (b2+b3*x)*bb - b4*aa ) / b1
    end

    b, (α+β+p+1)*aa/2
end

"""
    jacobi_vandermonde(nodes, α, β)

Computes the Vandermonde matrix with respect to the Jacobi polynomials with
parameters `α`, `β` and the nodal basis on `nodes`.
The Vandermonde matrix `V` is the transformation matrix from the modal Jacobi
basis to the nodal Lagrange basis associated with `nodes`.
"""
function jacobi_vandermonde(nodes::AbstractVector, α, β)
    T = eltype(nodes)
    pp1 = length(nodes)
    V = Array{T}(undef, pp1, pp1)
    for j in 1:pp1, (i,x) in enumerate(nodes)
      V[i, j] = jacobi(x, j-1, α, β)
    end
    V
end

jacobi_vandermonde(basis::NodalBasis{Line}, α, β) = jacobi_vandermonde(grid(basis), α, β)
