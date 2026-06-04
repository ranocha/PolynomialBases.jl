module PolynomialBases

using LinearAlgebra: Diagonal, rmul!

using ArgCheck: @argcheck
using AutoHashEquals: @auto_hash_equals
using FFTW: FFTW
using FastGaussQuadrature: FastGaussQuadrature, gausslegendre, gausslobatto, gaussjacobi, gaussradau
using SimpleUnPack: @unpack


# types
abstract type AbstractDomain{Dim} end
struct Line <: AbstractDomain{1} end

abstract type AbstractBasis{Domain, T} end
abstract type NodalBasis{Domain, T} <: AbstractBasis{Domain, T} end


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
export NodalBasis, LobattoLegendre, GaussLegendre, GaussRadauLeft, GaussRadauRight, GaussJacobi,
        ClosedNewtonCotes, ClenshawCurtis

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

## Continuous orthogonal polynomials
export legendre, legendre_vandermonde, legendre_D, legendre_M,
       gegenbauer,
       jacobi, jacobi_vandermonde,
       hermite,
       laguerre

## Discrete orthogonal polynomials
export hahn

## other utilities
export utility_matrices, grid, derivative_matrix, mass_matrix, mass_matrix_boundary,
       includes_boundaries, includes_left_boundary, includes_right_boundary,
       satisfies_sbp



end # module
