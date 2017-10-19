module PolynomialBases


using Parameters
using FastGaussQuadrature


abstract type AbstractDomain{Dim} end
struct Line <: AbstractDomain{1} end

abstract type AbstractBasis{Domain} end
abstract type NodalBasis{Domain} <: AbstractBasis{Domain} end


include("interpolation.jl")
include("integration.jl")
include("derivative.jl")
include("nodal_bases.jl")
include("legendre.jl")

# types
export LobattoLegendre, GaussLegendre

# interpolation
export interpolate, interpolate!, interpolation_matrix, interpolation_matrix!,
        change_basis, change_basis!

# derivative
export derivative_at, derivative_at!, derivative_matrix, derivative_matrix!

# integration
export integrate

# Legender polynomials
export legendre

# other utilities
export utility_matrices


end # module
