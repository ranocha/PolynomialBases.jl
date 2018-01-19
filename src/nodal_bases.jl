"""
    degree(basis::NodalBasis{Line})

Return the polynomial degree used by `basis`.
"""
Base.@pure degree(basis::NodalBasis{Line}) = length(basis.nodes)-1


"""
    change_basis{Domain<:AbstractDomain}(dest_basis::NodalBasis{Domain},
                                         values, src_basis::NodalBasis{Domain})

Perform the change of basis for the coefficients `values` from `src_basis` to
`dest_basis`.
"""
function change_basis{Domain<:AbstractDomain}(dest_basis::NodalBasis{Domain},
                                                values, src_basis::NodalBasis{Domain})
    @boundscheck begin
        @assert length(dest_basis.nodes) == length(src_basis.nodes) == length(values)
    end
    ret = similar(values)
    @inbounds change_basis!(ret, dest_basis, values, src_basis)
    ret
end

"""
    change_basis!{Domain<:AbstractDomain}(ret, dest_basis::NodalBasis{Domain},
                                          values, src_basis::NodalBasis{Domain})

Perform the change of basis for the coefficients `values` from `src_basis` to
`dest_basis` and store the resulting coefficients in `ret`.
"""
function change_basis!{Domain<:AbstractDomain}(ret, dest_basis::NodalBasis{Domain},
                                               values, src_basis::NodalBasis{Domain})
    @boundscheck begin
        @assert length(dest_basis.nodes) == length(src_basis.nodes) == length(values)
        @assert length(values) == length(ret)
    end
    interpolate!(ret, dest_basis.nodes, values, src_basis)
    nothing
end


"""
    compute_coefficients(u, basis::NodalBasis{Line})

Compute the nodal values of the function `u` at the nodes corresponding to the
nodal basis `basis`.
"""
function compute_coefficients(u, basis::NodalBasis{Line})
    xmin = first(basis.nodes)
    xmax = last(basis.nodes)
    uval = Array{typeof(u((xmin+xmax)/2))}(length(basis.nodes))
    compute_coefficients!(uval, u, basis)
    uval
end

"""
    compute_coefficients!(uval::AbstractVector, u, basis::NodalBasis{Line})

Compute the nodal values of the function `u` at the nodes corresponding to the
nodal basis `basis` and store the result in `uval`.
"""
function compute_coefficients!(uval::AbstractVector, u, basis::NodalBasis{Line})
    uval .= u.(basis.nodes)
    nothing
end


"""
    evaluate_coefficients(u, basis::NodalBasis{Line}, npoints=2*length(basis.nodes))

Evaluate the coefficients `u` of the nodal basis `basis` at `npoints` equally
spaced nodes. Returns `xplot, uplot`, where `xplot` contains the equally spaced
nodes and `uplot` the corresponding values of `u`.
"""
@noinline function evaluate_coefficients(u, basis::NodalBasis{Line}, npoints=2*length(basis.nodes))
    xplot = Array{eltype(basis.nodes)}(npoints)
    uplot = Array{eltype(u)}(npoints)

    evaluate_coefficients!(xplot, uplot, u, basis)
end

"""
    evaluate_coefficients!(xplot, uplot, u, basis::NodalBasis{Line})

Evaluate the coefficients `u` of the nodal basis `basis` at `npoints` equally
spaced nodes and store the result in `xplot, uplot`. Returns `xplot, uplot`,
where `xplot` contains the equally spaced nodes and `uplot` the corresponding
values of `u`.
"""
function evaluate_coefficients!(xplot, uplot, u, basis::NodalBasis{Line})
    @argcheck length(uplot) == length(xplot)
    npoints = length(xplot)
    T = eltype(xplot)

    xplot .= linspace(T(-1), T(1), npoints)
    interpolate!(uplot, xplot, u, basis)

    xplot, uplot
end


"""
    utility_matrices(basis::NodalBasis{Line})

Return the matrices
- `D`, derivative matrix
- `M`, mass matrix
- `R`, restriction matrix (interpolation to the boundaries)
- `B`, boundary normal matrix
- `MinvRtB = M \ R' * B`
used in the formulation of flux reconstruction / correction procedure via
reconstruction using summation-by-parts operators, cf. Ranocha, Öffner, Sonar
(2016) Summation-by-parts operators for correction procedure via reconstruction,
Journal of Computational Physics 311, 299-328.
"""
function utility_matrices(basis::NodalBasis{Line})
    D = basis.D
    M = Diagonal(basis.weights)
    R = interpolation_matrix([-1,1], basis)
    B = Diagonal([-1,1])
    MinvRtB = M \ R' * B

    D, M, R, B, MinvRtB
end



"""
    LobattoLegendre{T<:Real}

The nodal basis corresponding to Legendre Gauss Lobatto quadrature in [-1,1]
with scalar type `T`.
"""
struct LobattoLegendre{T} <: NodalBasis{Line}
    nodes::Vector{T}
    weights::Vector{T}
    baryweights::Vector{T}
    D::Matrix{T}

    function LobattoLegendre(nodes::Vector{T}, weights::Vector{T}, baryweights::Vector{T}, D::Matrix{T}) where T
        @assert length(nodes) == length(weights) == length(baryweights) == size(D,1) == size(D,2)
        new{T}(nodes, weights, baryweights, D)
    end
end

"""
    LobattoLegendre(p::Int, T=Float64)

Generate the `LobattoLegendre` basis of degree `p` with scalar type `T`.
"""
function LobattoLegendre(p::Int, T=Float64)
    if p == 0
        nodes = T[0]
        weights = T[2]
    elseif T == Float64
        nodes, weights = gausslobatto(p+1)
    else
        nodes, weights = lobatto_legendre_nodes_and_weights(p, T)
    end
    baryweights = barycentric_weights(nodes)
    D = derivative_matrix(nodes, baryweights)
    LobattoLegendre(nodes, weights, baryweights, D)
end

@require SymPy begin
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

@require SymEngine begin
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

function Base.show{T}(io::IO, basis::LobattoLegendre{T})
  print(io, "LobattoLegendre{", T, "}: Nodal Lobatto Legendre basis of degree ",
            degree(basis))
end


"""
    GaussLegendre{T<:Real}

The nodal basis corresponding to Legendre Gauss quadrature in [-1,1]
with scalar type `T`.
"""
struct GaussLegendre{T} <: NodalBasis{Line}
    nodes::Vector{T}
    weights::Vector{T}
    baryweights::Vector{T}
    D::Matrix{T}
    interp_left::Vector{T}
    interp_right::Vector{T}

    function GaussLegendre(nodes::Vector{T}, weights::Vector{T}, baryweights::Vector{T}, D::Matrix{T}, interp_left::Vector{T}, interp_right::Vector{T}) where T
        @assert length(nodes) == length(weights) == length(baryweights) == size(D,1) == size(D,2) == length(interp_left) == length(interp_right)
        new{T}(nodes, weights, baryweights, D, interp_left, interp_right)
    end
end

"""
    GaussLegendre(p::Int, T=Float64)

Generate the `GaussLegendre` basis of degree `p` with scalar type `T`.
"""
function GaussLegendre(p::Int, T=Float64)
    if T == Float64
        nodes, weights = gausslegendre(p+1)
    else
        nodes, weights = gauss_legendre_nodes_and_weights(p, T)
    end
    baryweights = barycentric_weights(nodes)
    D = derivative_matrix(nodes, baryweights)
    R = interpolation_matrix(T[-1, 1], nodes, baryweights)
    GaussLegendre(nodes, weights, baryweights, D, R[1,:], R[2,:])
end

@require SymPy begin
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

@require SymEngine begin
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

function Base.show{T}(io::IO, basis::GaussLegendre{T})
  print(io, "GaussLegendre{", T, "}: Nodal Gauss Legendre basis of degree ",
            degree(basis))
end


"""
    GaussJacobi{T<:Real}

The nodal basis corresponding to Jacobi Gauss quadrature in [-1,1]
with parameters `α`, `β` and scalar type `T`.
"""
struct GaussJacobi{T1<:Real, T<:Real} <: NodalBasis{Line}
    α::T1
    β::T1
    nodes::Vector{T}
    weights::Vector{T}
    baryweights::Vector{T}
    D::Matrix{T}

    function GaussJacobi(α::T1, β::T1, nodes::Vector{T}, weights::Vector{T}, baryweights::Vector{T}, D::Matrix{T}) where {T1,T}
        @assert length(nodes) == length(weights) == length(baryweights) == size(D,1) == size(D,2)
        new{T1, T}(α, β, nodes, weights, baryweights, D)
    end
end

"""
    GaussJacobi(p::Int, α, β, T=Float64)

Generate the `JacobiLegendre` basis of degree `p` with parameters `α`, `β` and
scalar type `T`.
"""
function GaussJacobi(p::Int, α, β, T=Float64)
    if T == Float64
        nodes, weights = gaussjacobi(p+1, α, β)
    else
        nodes, weights = gauss_jacobi_nodes_and_weights(p, α, β, T)
    end
    baryweights = barycentric_weights(nodes)
    D = derivative_matrix(nodes, baryweights)
    GaussJacobi(promote(α, β)..., nodes, weights, baryweights, D)
end

GaussJacobi(p::Int, T=Float64) = GaussJacobi(p, 0, 0, T)

function Base.show{T}(io::IO, basis::GaussJacobi{T})
  print(io, "GaussJacobi{", T, "}: Nodal Gauss Jacobi basis of degree ",
            degree(basis), " with parameters α=", basis.α, " and β = ", basis.β)
end
