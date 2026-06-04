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
