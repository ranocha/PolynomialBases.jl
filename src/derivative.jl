
"""
    derivative_at(x::Real, values, nodes, baryweights)

Compute the derivative of the function represented by `values` on the `nodes`
at `x` using the corresponding barycentric weights `baryweights`.
[Kopriva, Implementing Spectral Methods for PDEs, Algorithm 36].
"""
function derivative_at(x::Real, values, nodes, baryweights)
    @boundscheck begin
        @assert size(values) == size(nodes) == size(baryweights)
    end
    T = promote_type(typeof(x), eltype(values), eltype(nodes), eltype(baryweights))
    at_node = false
    num = zero(T)
    den = zero(T)
    i = typemin(Int)
    @inbounds for j in eachindex(nodes)
        if x ≈ nodes[j]
            at_node = true
            p = values[j]
            den = -baryweights[j]
            i = j
        end
    end
    if at_node
        @inbounds for j in eachindex(nodes)
            if j != i
                num += baryweights[j] * (p - values[j]) / (x - nodes[j])
            end
        end
    else
        @inbounds p = interpolate(x, values, nodes, baryweights)
        @inbounds for j in eachindex(nodes)
            t = baryweights[j] / (x - nodes[j])
            num += t * (p - values[j]) / (x - nodes[j])
            den += t
        end
    end
    num / den
end

function derivative_at(x::AbstractArray, values, nodes, baryweights)
    @boundscheck begin
        @assert size(values) == size(nodes) == size(baryweights)
    end
    T = promote_type(eltype(x), eltype(values), eltype(nodes), eltype(baryweights))
    res = Array{T}(undef, size(x))
    @inbounds derivative_at!(res, x, values, nodes, baryweights)
    res
end

function derivative_at!(res::AbstractArray, x::AbstractArray, values, nodes, baryweights)
    @boundscheck begin
        @assert size(res) == size(x)
        @assert size(values) == size(nodes) == size(baryweights)
    end
    @inbounds for idx in eachindex(res)
        res[idx] = derivative_at(x[idx], values, nodes, baryweights)
    end
    nothing
end

function derivative_at(x, values, basis::NodalBasis)
    @unpack nodes, baryweights = basis
    derivative_at(x, values, nodes, baryweights)
end

function derivative_at!(res, x, values, basis::NodalBasis)
    @unpack nodes, baryweights = basis
    derivative_at!(res, x, values, nodes, baryweights)
end


"""
    derivative_matrix(nodes, baryweights)

Compute the derivative matrix corresponding to `nodes` and `baryweights`.
[Kopriva, Implementing Spectral Methods for PDEs, Algorithm 37].
"""
function derivative_matrix(nodes, baryweights)
    @boundscheck begin
        @assert length(nodes) == length(baryweights)
    end
    T = promote_type(eltype(nodes), eltype(baryweights))
    mat = Array{T}(undef, length(nodes), length(nodes))
    @inbounds derivative_matrix!(mat, nodes, baryweights)
    mat
end

derivative_matrix(nodes) = derivative_matrix(nodes, barycentric_weights(nodes))

function derivative_matrix!(D::AbstractMatrix, nodes, baryweights)
    @boundscheck begin
        @assert size(D,1) == size(D,2) == length(nodes) == length(baryweights)
    end
    @inbounds for i in eachindex(nodes)
        D[i,i] = 0
        for j in eachindex(nodes)
            if i != j
                D[i,j] = baryweights[j] / baryweights[i] / (nodes[i] - nodes[j])
                D[i,i] -= D[i,j]
            end
        end
    end
    nothing
end
