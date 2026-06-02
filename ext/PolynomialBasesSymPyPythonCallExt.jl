module PolynomialBasesSymPyPythonCallExt

using PolynomialBases: PolynomialBases
using SymPyPythonCall: SymPyPythonCall

# SymPyPythonCall.Sym is the UnionAll Sym{T} where T; the concrete type
# used at runtime is Sym{PythonCall.Py}.
# This mirrors the pattern used in BSeries.jl.
const _Sym = SymPyPythonCall.Sym{SymPyPythonCall.PythonCall.Core.Py}

include("general_symbolic_extension.jl")

function PolynomialBases.interpolation_matrix!(mat, dest, src::AbstractVector{<:SymPyPythonCall.Sym}, baryweights)
    symbolic_interpolation_matrix!(mat, dest, src, baryweights, SymPyPythonCall.simplify)
end

PolynomialBases.LobattoLegendre(p::Int, ::Type{<:SymPyPythonCall.Sym}) = symbolic_lobatto_legendre(p, _Sym, SymPyPythonCall.simplify)
PolynomialBases.GaussLegendre(p::Int, ::Type{<:SymPyPythonCall.Sym}) = symbolic_gauss_legendre(p, _Sym, SymPyPythonCall.simplify)
PolynomialBases.GaussRadauLeft(p::Int, ::Type{<:SymPyPythonCall.Sym}) = symbolic_gauss_radau_left(p, _Sym, SymPyPythonCall.simplify)
PolynomialBases.GaussRadauRight(p::Int, ::Type{<:SymPyPythonCall.Sym}) = symbolic_gauss_radau_right(p, _Sym, SymPyPythonCall.simplify)

end # module
