"""
    gegenbauer(x, p::Integer)

Evaluate the Gegenbauer polynomial with parameter `α` of degree `p` at `x` using
the three term recursion.
"""
function gegenbauer(x, p::Integer, α)
    T = typeof( 2α*x )
    p₀ = one(T)
    p₁ = 2α*x

    if p <= 0
        return p₀
    elseif p == 1
        return p₁
    end

    for n in 2:p
        p₀, p₁ = p₁, ( 2x*(n+α-1)*p₁ - (2α+n-2)*p₀ ) / n
    end

    p₁
end
