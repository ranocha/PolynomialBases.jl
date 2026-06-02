module PolynomialBasesSymEngineExt

using PolynomialBases: PolynomialBases, barycentric_weights, derivative_matrix,
                       interpolation_matrix, LobattoLegendre, GaussLegendre,
                       GaussRadauLeft, GaussRadauRight
using SymEngine: SymEngine

include("general_symbolic_extension.jl")

function PolynomialBases.interpolation_matrix!(mat, dest, src::AbstractVector{SymEngine.Basic}, baryweights)
    symbolic_interpolation_matrix!(mat, dest, src, baryweights, SymEngine.expand)
end

PolynomialBases.LobattoLegendre(p::Int, T::Type{SymEngine.Basic}) = symbolic_lobatto_legendre(p, T, SymEngine.expand)
PolynomialBases.GaussLegendre(p::Int, T::Type{SymEngine.Basic}) = symbolic_gauss_legendre(p, T, SymEngine.expand)
PolynomialBases.GaussRadauLeft(p::Int, T::Type{SymEngine.Basic}) = symbolic_gauss_radau_left(p, T, SymEngine.expand)
PolynomialBases.GaussRadauRight(p::Int, T::Type{SymEngine.Basic}) = symbolic_gauss_radau_right(p, T, SymEngine.expand)

end # module
