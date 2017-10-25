"""
    map_to_canonical(x, xmin, xmax, basis::NodalBasis{Line})

Map the given coordinate `x` in the interval [`xmin`, `xmax`] to the corresponding
canonical coordinate `ξ` in the interval [-1, 1].
"""
function map_to_canonical(x, xmin, xmax, basis::NodalBasis{Line})
    ( 2x - (xmax + xmin) ) / (xmax - xmin)
end

"""
    map_to_canonical!(ξ, x, xmin, xmax, basis::NodalBasis{Line})

Map the given coordinate `x` in the interval [`xmin`, `xmax`] to the corresponding
canonical coordinate `ξ` in the interval [-1, 1], updating `ξ`.
"""
function map_to_canonical!(ξ, x, xmin, xmax, basis::NodalBasis{Line})
    @boundscheck begin
        @assert size(ξ) == size(x)
    end
    @inbounds for i in eachindex(ξ)
        ξ[i] = map_to_canonical(x[i], xmin, xmax, basis)
    end
    ξ
end


"""
    map_from_canonical(ξ, xmin, xmax, basis::NodalBasis{Line})

Map the given canonical coordinate `ξ` in the interval [-1, 1] to the corresponding
coordinate `x` in the interval [`xmin`, `xmax`].
"""
function map_from_canonical(ξ, xmin, xmax, basis::NodalBasis{Line})
    ( (xmax + xmin) + ξ * (xmax - xmin) ) / 2
end

"""
    map_from_canonical!(x, ξ, xmin, xmax, basis::NodalBasis{Line})

Map the given canonical coordinate `ξ` in the interval [-1, 1] to the corresponding
coordinate `x` in the interval [`xmin`, `xmax`], updating `x`.
"""
function map_from_canonical!(x, ξ, xmin, xmax, basis::NodalBasis{Line})
    @boundscheck begin
        @assert size(ξ) == size(x)
    end
    @inbounds for i in eachindex(x)
        x[i] = map_from_canonical(ξ[i], xmin, xmax, basis)
    end
    x
end
