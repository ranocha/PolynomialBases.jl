__precompile__()

module PolynomialBases

using ArgCheck
using Requires
using Parameters
using FastGaussQuadrature


# types
abstract type AbstractDomain{Dim} end
struct Line <: AbstractDomain{1} end

abstract type AbstractBasis{Domain} end
abstract type NodalBasis{Domain} <: AbstractBasis{Domain} end


# source files
include("canonical_mappings.jl")
include("interpolation.jl")
include("integration.jl")
include("derivative.jl")
include("nodal_bases.jl")
include("legendre.jl")
include("gegenbauer.jl")
include("jacobi.jl")
include("hermite.jl")
include("laguerre.jl")
include("hahn.jl")

# export
## types
export NodalBasis, LobattoLegendre, GaussLegendre, GaussJacobi

## mappings
export map_to_canonical, map_to_canonical!, map_from_canonical, map_from_canonical!

## interpolation
export interpolate, interpolate!, interpolation_matrix, interpolation_matrix!,
        change_basis, change_basis!

## evaluation of coefficients
export compute_coefficients, compute_coefficients!,
        evaluate_coefficients, evaluate_coefficients!

## derivative
export derivative_at, derivative_at!, derivative_matrix, derivative_matrix!

## integration
export integrate

## Continuous orthogonal polynomials: Legendre, Gegenbauer, Jacobi, Hermite
export legendre, legendre_vandermonde, legendre_D, legendre_M,
       gegenbauer,
       jacobi, jacobi_vandermonde, jacobi_M,
       hermite,
       laguerre

## Discrete orthogonal polynomials: Hahn
export hahn

## other utilities
export utility_matrices


end # module
