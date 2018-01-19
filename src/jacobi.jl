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
    gauss_jacobi_nodes_and_weights(p, α, β, T=Float64::Type, tol=4*eps(T), maxit=100)

Compute the Gauss-Jacobi nodes and weights for polynomials of degree `p` with
parameters `α`, `β` using the scalar type `T`, tolerance `tol` and maximal number
of Newton iterations `maxit` [Karniadakis and Sherwin, Spectral/hp Element
Methods for CFD, Appendix B].
"""
function gauss_jacobi_nodes_and_weights(p, α, β, T=Float64::Type, tol=4*eps(T), maxit=1000)
    T = promote_type(typeof(α), typeof(β), T)
    nodes = Vector{T}(p+1)
    weights = Vector{T}(p+1)
    for j in 0:p
        x = -cospi(T(2j+1)/(2p+2))
        if j > 0
            x = (x + nodes[j]) / 2
        end
        for k in 1:maxit
            s = zero(T)
            for l in 1:j
                s += 1 / (x - nodes[l])
            end
            pol, der = jacobi_and_derivative(x, p+1, α, β)
            Δ = -pol / (der - s*pol)
            x = x + Δ
            abs(Δ) <= tol*abs(x) && break
        end
        pol, der = jacobi_and_derivative(x, p+1, α, β)
        nodes[j+1] = x
        weights[j+1] = 2^(α+β+1) * gamma(α+p+2) * gamma(β+p+2) /
                        (gamma(p+2) * gamma(α+β+p+2) * (1-x^2) * der^2)
    end

    return nodes, weights
end



doc"
    jacobi_vandermonde(nodes, α, β)

Computes the Vandermonde matrix with respect to the Jacobi polynomials with
parameters `α`, `β` and the nodal basis on `nodes`.
The Vandermonde matrix $V$ is the transformation matrix from the modal Jacobi
basis to the nodal Lagrange basis associated with `nodes`.
"
function jacobi_vandermonde(nodes::AbstractVector, α, β)
    T = eltype(nodes)
    pp1 = length(nodes)
    V = Array{T}(pp1, pp1)
    for j in 1:pp1, (i,x) in enumerate(nodes)
      V[i, j] = jacobi(x, j-1, α, β)
    end
    V
end

jacobi_vandermonde(basis::NodalBasis{Line}, α, β) = jacobi_vandermonde(basis.nodes, α, β)
