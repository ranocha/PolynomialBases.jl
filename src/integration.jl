
"""
    integrate(func, u, weights)

Map the function `func` to the coefficients `u` and integrate with respect to
the quadrature rule given by `weights`.
"""
function integrate(func, u, weights)
    Pp1 = length(weights)
    @boundscheck begin
        @assert Pp1 == length(u)
    end

    res = zero(func(first(u)))
    @inbounds for n in 1:Pp1
        res += func(u[n]) * weights[n]
    end
    res
end

function integrate(func, u, basis::NodalBasis)
    integrate(func, u, basis.weights)
end

integrate(u, weights) = integrate(identity, u, weights)
integrate(u, basis::NodalBasis) = integrate(identity, u, basis)
