module PolynomialBasesSymPyExt

using PolynomialBases: PolynomialBases
using SymPy: SymPy

include("general_symbolic_extension.jl")

# SymPy v1 defines a concrete non-parametric Sym struct.
# SymPy v2 reexports Sym from SymPyCore, making it a UnionAll Sym{T} where T;
# the concrete type is Sym{PyCall.PyObject}, following the pattern in BSeries.jl.
# isconcretetype distinguishes the two at module-load time.
if isconcretetype(SymPy.Sym)
    const _Sym = SymPy.Sym
else
    const _Sym = SymPy.Sym{SymPy.PyCall.PyObject}
end

function PolynomialBases.interpolation_matrix!(mat, dest, src::AbstractVector{<:SymPy.Sym}, baryweights)
    symbolic_interpolation_matrix!(mat, dest, src, baryweights, SymPy.simplify)
end

PolynomialBases.LobattoLegendre(p::Int, ::Type{<:SymPy.Sym}) = symbolic_lobatto_legendre(p, _Sym, SymPy.simplify)
PolynomialBases.GaussLegendre(p::Int, ::Type{<:SymPy.Sym}) = symbolic_gauss_legendre(p, _Sym, SymPy.simplify)
PolynomialBases.GaussRadauLeft(p::Int, ::Type{<:SymPy.Sym}) = symbolic_gauss_radau_left(p, _Sym, SymPy.simplify)
PolynomialBases.GaussRadauRight(p::Int, ::Type{<:SymPy.Sym}) = symbolic_gauss_radau_right(p, _Sym, SymPy.simplify)

end # module
