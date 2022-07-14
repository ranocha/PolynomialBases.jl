"""
    degree(basis::NodalBasis{Line})

Return the polynomial degree used by `basis`.
"""
Base.@pure degree(basis::NodalBasis{Line}) = length(basis.nodes)-1


"""
    function change_basis(dest_basis::NodalBasis{Domain},
                          values, src_basis::NodalBasis{Domain}) where {Domain<:AbstractDomain}

Perform the change of basis for the coefficients `values` from `src_basis` to
`dest_basis`.
"""
function change_basis(dest_basis::NodalBasis{Domain},
                      values, src_basis::NodalBasis{Domain}) where {Domain<:AbstractDomain}
    @boundscheck begin
        @assert length(dest_basis.nodes) == length(src_basis.nodes) == length(values)
    end
    ret = similar(values)
    @inbounds change_basis!(ret, dest_basis, values, src_basis)
    ret
end

"""
    change_basis!(ret, dest_basis::NodalBasis{Domain},
                  values, src_basis::NodalBasis{Domain}) where {Domain<:AbstractDomain}

Perform the change of basis for the coefficients `values` from `src_basis` to
`dest_basis` and store the resulting coefficients in `ret`.
"""
function change_basis!(ret, dest_basis::NodalBasis{Domain},
                       values, src_basis::NodalBasis{Domain}) where {Domain<:AbstractDomain}
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
    uval = Array{typeof(u((xmin+xmax)/2))}(undef, length(basis.nodes))
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
function evaluate_coefficients(u, basis::NodalBasis{Line}, npoints=2*length(basis.nodes))
    xplot = Array{eltype(basis.nodes)}(undef, npoints)
    uplot = Array{eltype(u)}(undef, npoints)

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

    xplot .= range(T(-1), stop=T(1), length=npoints)
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
- `MinvRtB = M \\ R' * B`
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


# for broadcasting; treat bases as scalars
Base.broadcastable(basis::NodalBasis) = Ref(basis)



"""
    LobattoLegendre{T}

The nodal basis corresponding to Legendre Gauss Lobatto quadrature in [-1,1]
with scalar type `T`.
"""
struct LobattoLegendre{T} <: NodalBasis{Line}
    nodes::Vector{T}
    weights::Vector{T}
    baryweights::Vector{T}
    D::Matrix{T}

    function LobattoLegendre(nodes::Vector{T}, weights::Vector{T}, baryweights::Vector{T}, D::Matrix{T}) where {T}
        @argcheck length(nodes) == length(weights) == length(baryweights) == size(D,1) == size(D,2)
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
    else
        nodes, weights = lobatto_legendre_nodes_and_weights_impl(p, T)
    end
    baryweights = barycentric_weights(nodes)
    D = derivative_matrix(nodes, baryweights)
    LobattoLegendre(nodes, weights, baryweights, D)
end

function lobatto_legendre_nodes_and_weights_impl(p, ::Type{Float64})
    nodes::Vector{Float64}, weights::Vector{Float64} = gausslobatto(p+1)
    nodes, weights
end

function lobatto_legendre_nodes_and_weights_impl(p, T::DataType)
    nodes, weights = lobatto_legendre_nodes_and_weights(p, T)
    nodes, weights
end

# special methods of LobattoLegendre for SymPy and SymEngine are in __init__

function Base.show(io::IO, basis::LobattoLegendre{T}) where {T}
  print(io, "LobattoLegendre{", T, "}: Nodal Lobatto Legendre basis of degree ",
            degree(basis))
end

@inline includes_boundaries(basis::LobattoLegendre) = Val{true}()

@inline satisfies_sbp(basis::LobattoLegendre) = Val{true}()


"""
    GaussLegendre{T}

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

    function GaussLegendre(nodes::Vector{T}, weights::Vector{T}, baryweights::Vector{T}, D::Matrix{T}, interp_left::Vector{T}, interp_right::Vector{T}) where {T}
        @argcheck length(nodes) == length(weights) == length(baryweights) == size(D,1) == size(D,2) == length(interp_left) == length(interp_right)
        new{T}(nodes, weights, baryweights, D, interp_left, interp_right)
    end
end

"""
    GaussLegendre(p::Int, T=Float64)

Generate the `GaussLegendre` basis of degree `p` with scalar type `T`.
"""
function GaussLegendre(p::Int, T=Float64)
    nodes, weights = gauss_legendre_nodes_and_weights_impl(p, T)
    baryweights = barycentric_weights(nodes)
    D = derivative_matrix(nodes, baryweights)
    R = interpolation_matrix(T[-1, 1], nodes, baryweights)
    GaussLegendre(nodes, weights, baryweights, D, R[1,:], R[2,:])
end

function gauss_legendre_nodes_and_weights_impl(p, ::Type{Float64})
    nodes::Vector{Float64}, weights::Vector{Float64} = gausslegendre(p+1)
    nodes, weights
end

function gauss_legendre_nodes_and_weights_impl(p, T::DataType)
    nodes, weights = gauss_legendre_nodes_and_weights(p, T)
    nodes, weights
end

# special methods of GaussLegendre for SymPy and SymEngine are in __init__

function Base.show(io::IO, basis::GaussLegendre{T}) where {T}
  print(io, "GaussLegendre{", T, "}: Nodal Gauss Legendre basis of degree ",
            degree(basis))
end

@inline includes_boundaries(basis::GaussLegendre) = Val{false}()

@inline satisfies_sbp(basis::GaussLegendre) = Val{true}()


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
        @argcheck length(nodes) == length(weights) == length(baryweights) == size(D,1) == size(D,2)
        new{T1, T}(α, β, nodes, weights, baryweights, D)
    end
end

"""
    GaussJacobi(p::Int, α, β, T=Float64)

Generate the `JacobiLegendre` basis of degree `p` with parameters `α`, `β` and
scalar type `T`.
"""
function GaussJacobi(p::Int, α, β, T=Float64)
    nodes, weights = gauss_jacobi_nodes_and_weights_impl(p, α, β, T)
    baryweights = barycentric_weights(nodes)
    D = derivative_matrix(nodes, baryweights)
    GaussJacobi(promote(α, β)..., nodes, weights, baryweights, D)
end

function gauss_jacobi_nodes_and_weights_impl(p, α, β, ::Type{Float64})
    nodes::Vector{Float64}, weights::Vector{Float64} = gaussjacobi(p+1, α, β)
    nodes, weights
end

function gauss_jacobi_nodes_and_weights_impl(p, α, β, T::DataType)
    nodes, weights = gauss_jacobi_nodes_and_weights(p, α, β, T)
    nodes, weights
end

GaussJacobi(p::Int, T=Float64) = GaussJacobi(p, 0, 0, T)

function Base.show(io::IO, basis::GaussJacobi{T1, T}) where {T1, T}
  print(io, "GaussJacobi{", T1, ", ", T, "}: Nodal Gauss Jacobi basis of degree ",
            degree(basis), " with parameters α=", basis.α, " and β = ", basis.β)
end

@inline includes_boundaries(basis::GaussJacobi) = Val{false}()

@inline satisfies_sbp(basis::GaussJacobi) = Val{false}()


"""
    ClosedNewtonCotes{T}

The nodal basis corresponding to the closed Newton Cotes quadrature in [-1,1]
with scalar type `T`.
"""
struct ClosedNewtonCotes{T} <: NodalBasis{Line}
    nodes::Vector{T}
    weights::Vector{T}
    baryweights::Vector{T}
    D::Matrix{T}

    function ClosedNewtonCotes(nodes::Vector{T}, weights::Vector{T}, baryweights::Vector{T}, D::Matrix{T}) where {T}
        @argcheck length(nodes) == length(weights) == length(baryweights) == size(D,1) == size(D,2)
        new{T}(nodes, weights, baryweights, D)
    end
end

"""
    ClosedNewtonCotes(p::Int, T=Float64)

Generate the `ClosedNewtonCotes` basis of degree `p` with scalar type `T`.
"""
function ClosedNewtonCotes(p::Int, T=Float64)
    if p == 0
        nodes = T[0]
        weights = T[2]
    elseif 1 <= p <= 7
        nodes = T.(-1:2//p:1)
        if p == 1
            weights = T[1, 1]
        elseif p == 2
            weights = T[1//3, 4//3, 1//3]
        elseif p == 3
            weights = T[2//8, 6//8, 6//8, 6//8]
        elseif p == 4
            weights = T[14//90, 32//45, 4//15, 32//45, 14//90]
        elseif p == 5
            weights = T[38//288, 50//96, 50//144, 50//144, 50//96, 38//288]
        elseif p == 6
            weights = T[82//840, 18//35, 18//280, 68//105, 18//280, 18//35, 82//840]
        elseif p == 7
            weights = T[2*751//17280, 2*3577//17280, 98//640, 2*2989//17280, 2*2989//17280, 98//640, 2*3577//17280, 2*751//17280]
        end
    else
        throw(DomainError("Only implemented for p <= 6 due to nonnegative weights."))
    end
    baryweights = barycentric_weights(nodes)
    D = derivative_matrix(nodes, baryweights)
    ClosedNewtonCotes(nodes, weights, baryweights, D)
end

@inline includes_boundaries(basis::ClosedNewtonCotes) = Val{true}()

@inline satisfies_sbp(basis::ClosedNewtonCotes) = Val{false}()


"""
    ClenshawCurtis{T}

The nodal basis corresponding to the Clenshaw Curtis quadrature in [-1,1] with
scalar type `T`.
"""
struct ClenshawCurtis{T} <: NodalBasis{Line}
    nodes::Vector{T}
    weights::Vector{T}
    baryweights::Vector{T}
    D::Matrix{T}

    function ClenshawCurtis(nodes::Vector{T}, weights::Vector{T}, baryweights::Vector{T}, D::Matrix{T}) where {T}
        @argcheck length(nodes) == length(weights) == length(baryweights) == size(D,1) == size(D,2)
        new{T}(nodes, weights, baryweights, D)
    end
end

"""
    ClenshawCurtis(p::Int, T=Float64)

Generate the `ClenshawCurtis` basis of degree `p` with scalar type `T`.
"""
function ClenshawCurtis(p::Int, T=Float64)
    if p == 0
        nodes = T[0]
        weights = T[2]
    else
        nodes = Vector{T}(undef, p + 1)
        for k in eachindex(nodes)
            nodes[k] = sinpi(convert(T, p + 2 - 2 * k) / (2 * p))
        end

        weights = zeros(T, p + 1)
        for i in 0:2:p
            @inbounds weights[i+1] = convert(T, 2) / convert(T, (1 - i^2))
        end
        rmul!(weights, inv(convert(T, p)))
        plan = FFTW.plan_r2r!(weights, FFTW.REDFT00)
        plan * weights
        half = inv(convert(T, 2))
        weights[1] *= half
        weights[end] *= half
    end
    baryweights = barycentric_weights(nodes)
    D = derivative_matrix(nodes, baryweights)
    ClenshawCurtis(nodes, weights, baryweights, D)
end

@inline includes_boundaries(basis::ClenshawCurtis) = Val{true}()

@inline satisfies_sbp(basis::ClenshawCurtis) = Val{false}()
