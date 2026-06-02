module PolynomialBasesSymPyExt

using PolynomialBases: PolynomialBases
using SymPy: SymPy

include("general_symbolic_extension.jl")

function PolynomialBases.interpolation_matrix!(mat, dest, src::AbstractVector{SymPy.Sym}, baryweights)
    symbolic_interpolation_matrix!(mat, dest, src, baryweights, SymPy.simplify)
end

PolynomialBases.LobattoLegendre(p::Int, T::Type{SymPy.Sym}) = symbolic_lobatto_legendre(p, T, SymPy.simplify)
PolynomialBases.GaussLegendre(p::Int, T::Type{SymPy.Sym}) = symbolic_gauss_legendre(p, T, SymPy.simplify)
PolynomialBases.GaussRadauLeft(p::Int, T::Type{SymPy.Sym}) = symbolic_gauss_radau_left(p, T, SymPy.simplify)
PolynomialBases.GaussRadauRight(p::Int, T::Type{SymPy.Sym}) = symbolic_gauss_radau_right(p, T, SymPy.simplify)

end # module
