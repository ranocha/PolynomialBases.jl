"""
    hahn(x, p::Integer, α, β, N::Integer)

Evaluate the Hahn polynomial with parameters `α`, `β`, `N` of degree `p` at `x`
using the three term recursion [Öffner, Zweidimensionale klassische und diskrete
orthogonale Polynome, Chapter 5].
"""
function hahn(x, p::Integer, α, β, N)
    T = typeof( (α+β+2)*x/(N*(α+1)) )
    a = one(T)
    b = 1 - (α+β+2)*x/(N*(α+1))

    if p <= 0
        return a
    elseif p == 1
        return b
    elseif typeof(N) <: Integer && p > N
        throw(ArgumentError("p==$p must not be larger than N==$N."))
    end

    for n in 2:p
        An = T( (n+α+β)*(n+α)*(N-n+1) ) / ( (2n-1+α+β)*(2n+α+β) )
        Cn = T( (n-1)*(n+α+β+N)*(n-1+β) ) / ( (2n-2+α+β)*(2n-1+α+β) )
        a, b = b, ( (An+Cn-x)*b - Cn*a ) / An
    end

    b
end
