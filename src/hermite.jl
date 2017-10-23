"""
    hermite(x, p::Integer)

Evaluate the Hermite polynomial of degree `p` at `x` using the three term recursion.
"""
function hermite(x, p::Integer)
    T = typeof( 2x )
    p₀ = one(T)
    p₁ = 2x

    if p <= 0
        return p₀
    elseif p == 1
        return p₁
    end

    for n in 2:p
        p₀, p₁ = p₁, 2x*p₁ - 2*(n-1)*p₀
    end

    p₁
end
