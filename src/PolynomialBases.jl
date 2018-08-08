__precompile__()

module PolynomialBases

using ArgCheck
using Requires
using Parameters
using LinearAlgebra
using SpecialFunctions
using FastGaussQuadrature
using FastTransforms: clenshawcurtisnodes, clenshawcurtisweights, chebyshevmoments1


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
export NodalBasis, LobattoLegendre, GaussLegendre, GaussJacobi, ClosedNewtonCotes,
        ClenshawCurtis

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
export utility_matrices, includes_boundaries, satisfies_sbp


function __init__()
    # interpolation_matrix
    @require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" begin
        function interpolation_matrix!(mat, dest, src::AbstractVector{SymPy.Sym}, baryweights)
            @boundscheck begin
                @assert length(src) == length(baryweights)
                @assert size(mat,1) == length(dest)
                @assert size(mat,2) == length(src)
            end
            fill!(mat, zero(eltype(mat)))
            @inbounds for k in 1:size(mat,1)
                row_has_match = false
                for j in 1:size(mat,2)
                    # no method matching rtoldefault(::Type{SymPy.Sym})
                    if dest[k] == src[j]
                        row_has_match = true
                        mat[k,j] = 1
                    end
                end
                if row_has_match == false
                    s = zero(eltype(mat))
                    for j in 1:size(mat,2)
                        t = baryweights[j] / (dest[k] - src[j])
                        mat[k,j] = t
                        s = s + t
                    end
                    for j in 1:size(mat,2)
                        mat[k,j] /= s
                    end
                end
            end
            @inbounds for idx in eachindex(mat)
                mat[idx] = SymPy.simplify(mat[idx])
            end
            nothing
        end
    end

    @require SymEngine="123dc426-2d89-5057-bbad-38513e3affd8" begin
        function interpolation_matrix!(mat, dest, src::AbstractVector{SymEngine.Basic}, baryweights)
            @boundscheck begin
                @assert length(src) == length(baryweights)
                @assert size(mat,1) == length(dest)
                @assert size(mat,2) == length(src)
            end
            fill!(mat, zero(eltype(mat)))
            @inbounds for k in 1:size(mat,1)
                row_has_match = false
                for j in 1:size(mat,2)
                    # no method matching rtoldefault(::Type{SymEngine.Basic})
                    if dest[k] == src[j]
                        row_has_match = true
                        mat[k,j] = 1
                    end
                end
                if row_has_match == false
                    s = zero(eltype(mat))
                    for j in 1:size(mat,2)
                        t = baryweights[j] / (dest[k] - src[j])
                        mat[k,j] = t
                        s = s + t
                    end
                    for j in 1:size(mat,2)
                        mat[k,j] /= s
                    end
                end
            end
            @inbounds for idx in eachindex(mat)
                mat[idx] = SymEngine.expand(mat[idx])
            end
            nothing
        end
    end

    # LobattoLegendre
    @require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" begin
        function LobattoLegendre(p::Int, T::Type{SymPy.Sym})
            if p == 0
                nodes   = T[0]
                weights = T[2]
            elseif p == 1
                nodes   = T[-1, 1]
                weights = T[1, 1]
            elseif p == 2
                nodes   = T[-1, 0, 1]
                weights = T[1//3, 4//3, 1//3]
            elseif p == 3
                sqrt_1_5 = sqrt(one(T) / 5)
                nodes   = T[-1, -sqrt_1_5, sqrt_1_5, 1]
                weights = T[1//6, 5//6, 5//6, 1//6]
            elseif p == 4
                sqrt_3_7 = sqrt(T(3) / 7)
                nodes   = T[-1, -sqrt_3_7, 0, sqrt_3_7, 1]
                weights = T[1//10, 49//90, 32//45, 49//90, 1//10]
            elseif p == 5
                sqrt_7 = sqrt(T(7))
                sqrt_m = sqrt(one(T)/3 - 2*sqrt_7/21)
                sqrt_p = sqrt(one(T)/3 + 2*sqrt_7/21)
                nodes   = T[-1, -sqrt_p, -sqrt_m, sqrt_m, sqrt_p, 1]
                weights = T[1//15, (14-sqrt_7)/30, (14+sqrt_7)/30, (14+sqrt_7)/30, (14-sqrt_7)/30, 1//15]
            elseif p == 6
                sqrt_5_3 = sqrt(T(5) / 3)
                sqrt_15 = sqrt(T(15))
                sqrt_m = sqrt((5 - 2*sqrt_5_3)/11)
                sqrt_p = sqrt((5 + 2*sqrt_5_3)/11)
                nodes   = T[-1, -sqrt_p, -sqrt_m, 0, sqrt_m, sqrt_p, 1]
                weights = T[1//21, (124-7*sqrt_15)/350, (124+7*sqrt_15)/350, 256//525, (124+7*sqrt_15)/350, (124-7*sqrt_15)/350, 1//21]
            else
                throw(ArgumentError("Polynomial degree p = $p not implemented yet."))
            end
            baryweights = SymPy.simplify.(barycentric_weights(nodes))
            D = SymPy.simplify.(derivative_matrix(nodes, baryweights))
            LobattoLegendre(nodes, weights, baryweights, D)
        end
    end

    @require SymEngine="123dc426-2d89-5057-bbad-38513e3affd8" begin
        function LobattoLegendre(p::Int, T::Type{SymEngine.Basic})
            if p == 0
                nodes   = T[0]
                weights = T[2]
            elseif p == 1
                nodes   = T[-1, 1]
                weights = T[1, 1]
            elseif p == 2
                nodes   = T[-1, 0, 1]
                weights = T[1//3, 4//3, 1//3]
            elseif p == 3
                sqrt_1_5 = sqrt(one(T) / 5)
                nodes   = T[-1, -sqrt_1_5, sqrt_1_5, 1]
                weights = T[1//6, 5//6, 5//6, 1//6]
            elseif p == 4
                sqrt_3_7 = sqrt(T(3) / 7)
                nodes   = T[-1, -sqrt_3_7, 0, sqrt_3_7, 1]
                weights = T[1//10, 49//90, 32//45, 49//90, 1//10]
            elseif p == 5
                sqrt_7 = sqrt(T(7))
                sqrt_m = sqrt(one(T)/3 - 2*sqrt_7/21)
                sqrt_p = sqrt(one(T)/3 + 2*sqrt_7/21)
                nodes   = T[-1, -sqrt_p, -sqrt_m, sqrt_m, sqrt_p, 1]
                weights = T[1//15, (14-sqrt_7)/30, (14+sqrt_7)/30, (14+sqrt_7)/30, (14-sqrt_7)/30, 1//15]
            elseif p == 6
                sqrt_5_3 = sqrt(T(5) / 3)
                sqrt_15 = sqrt(T(15))
                sqrt_m = sqrt((5 - 2*sqrt_5_3)/11)
                sqrt_p = sqrt((5 + 2*sqrt_5_3)/11)
                nodes   = T[-1, -sqrt_p, -sqrt_m, 0, sqrt_m, sqrt_p, 1]
                weights = T[1//21, (124-7*sqrt_15)/350, (124+7*sqrt_15)/350, 256//525, (124+7*sqrt_15)/350, (124-7*sqrt_15)/350, 1//21]
            else
                throw(ArgumentError("Polynomial degree p = $p not implemented yet."))
            end
            baryweights = SymEngine.expand.(barycentric_weights(nodes))
            D = SymEngine.expand.(derivative_matrix(nodes, baryweights))
            LobattoLegendre(nodes, weights, baryweights, D)
        end
    end

    # GaussLegendre
    @require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" begin
        function GaussLegendre(p::Int, T::Type{SymPy.Sym})
            if p == 0
                nodes   = T[0]
                weights = T[2]
            elseif p == 1
                sqrt_1_3 = sqrt(one(T) / 3)
                nodes   = T[-sqrt_1_3, sqrt_1_3]
                weights = T[1, 1]
            elseif p == 2
                sqrt_3_5 = sqrt(T(3) / 5)
                nodes   = T[-sqrt_3_5, 0, sqrt_3_5]
                weights = T[5//9, 8//9, 5//9]
            elseif p == 3
                sqrt_30 = sqrt(T(30))
                sqrt_6_5 = sqrt(T(6)/5)
                sqrt_m = sqrt((3-2*sqrt_6_5)/7)
                sqrt_p = sqrt((3+2*sqrt_6_5)/7)
                nodes   = T[-sqrt_p, -sqrt_m, sqrt_m, sqrt_p]
                weights = T[(18-sqrt_30)/36, (18+sqrt_30)/36, (18+sqrt_30)/36, (18-sqrt_30)/36]
            elseif p == 4
                sqrt_70 = sqrt(T(70))
                sqrt_10_7 = sqrt(T(10)/7)
                sqrt_m = sqrt(5 - 2*sqrt_10_7) / 3
                sqrt_p = sqrt(5 + 2*sqrt_10_7) / 3
                w_m = (322-13*sqrt_70) / 900
                w_p = (322+13*sqrt_70) / 900
                nodes   = T[-sqrt_p, -sqrt_m, 0, sqrt_m, sqrt_p]
                weights = T[w_m, w_p, 128//225, w_p, w_m]
            else
                throw(ArgumentError("Polynomial degree p = $p not implemented yet."))
            end
            baryweights = SymPy.simplify.(barycentric_weights(nodes))
            D = SymPy.simplify.(derivative_matrix(nodes, baryweights))
            R = interpolation_matrix([-1, 1], nodes, baryweights)
            GaussLegendre(nodes, weights, baryweights, D, R[1,:], R[2,:])
        end
    end

    @require SymEngine="123dc426-2d89-5057-bbad-38513e3affd8" begin
        function GaussLegendre(p::Int, T::Type{SymEngine.Basic})
            if p == 0
                nodes   = T[0]
                weights = T[2]
            elseif p == 1
                sqrt_1_3 = sqrt(one(T) / 3)
                nodes   = T[-sqrt_1_3, sqrt_1_3]
                weights = T[1, 1]
            elseif p == 2
                sqrt_3_5 = sqrt(T(3) / 5)
                nodes   = T[-sqrt_3_5, 0, sqrt_3_5]
                weights = T[5//9, 8//9, 5//9]
            elseif p == 3
                sqrt_30 = sqrt(T(30))
                sqrt_6_5 = sqrt(T(6)/5)
                sqrt_m = sqrt((3-2*sqrt_6_5)/7)
                sqrt_p = sqrt((3+2*sqrt_6_5)/7)
                nodes   = T[-sqrt_p, -sqrt_m, sqrt_m, sqrt_p]
                weights = T[(18-sqrt_30)/36, (18+sqrt_30)/36, (18+sqrt_30)/36, (18-sqrt_30)/36]
            elseif p == 4
                sqrt_70 = sqrt(T(70))
                sqrt_10_7 = sqrt(T(10)/7)
                sqrt_m = sqrt(5 - 2*sqrt_10_7) / 3
                sqrt_p = sqrt(5 + 2*sqrt_10_7) / 3
                w_m = (322-13*sqrt_70) / 900
                w_p = (322+13*sqrt_70) / 900
                nodes   = T[-sqrt_p, -sqrt_m, 0, sqrt_m, sqrt_p]
                weights = T[w_m, w_p, 128//225, w_p, w_m]
            else
                throw(ArgumentError("Polynomial degree p = $p not implemented yet."))
            end
            baryweights = SymEngine.expand.(barycentric_weights(nodes))
            D = SymEngine.expand.(derivative_matrix(nodes, baryweights))
            R = interpolation_matrix([-1, 1], nodes, baryweights)
            GaussLegendre(nodes, weights, baryweights, D, R[1,:], R[2,:])
        end
    end
end


end # module
