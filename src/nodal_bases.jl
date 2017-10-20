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
struct LobattoLegendre{T<:Real} <: NodalBasis{Line}
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

function Base.show{T}(io::IO, basis::LobattoLegendre{T})
  print(io, "LobattoLegendre{", T, "}: Nodal Lobatto Legendre basis of degree ",
            degree(basis))
end


"""
    GaussLegendre{T<:Real}

The nodal basis corresponding to Legendre Gauss quadrature in [-1,1]
with scalar type `T`.
"""
struct GaussLegendre{T<:Real} <: NodalBasis{Line}
    nodes::Vector{T}
    weights::Vector{T}
    baryweights::Vector{T}
    D::Matrix{T}

    function GaussLegendre(nodes::Vector{T}, weights::Vector{T}, baryweights::Vector{T}, D::Matrix{T}) where T
        @assert length(nodes) == length(weights) == length(baryweights) == size(D,1) == size(D,2)
        new{T}(nodes, weights, baryweights, D)
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
    GaussLegendre(nodes, weights, baryweights, D)
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
