
"""
    barycentric_weights{T<:Real}(x::AbstractVector{T})

Compute the barycentric weights corresponding to the nodes `x`.
[Kopriva, Implementing Spectral Methods for PDEs, Algorithm 30].
"""
function barycentric_weights{T<:Real}(x::AbstractVector{T})
    w = ones(x)
    @inbounds barycentric_weights!(w, x)
    w
end

function barycentric_weights!(w, x)
    @boundscheck begin
        length(w) == length(x)
    end
    fill!(w, 1)
    @inbounds for j in 2:length(w), k in 1:j-1
        diff = x[k] - x[j]
        w[k] *= diff
        w[j] *= -diff
    end
    @inbounds for j in 1:length(w)
        w[j] = 1 / w[j]
    end
    nothing
end


"""
    interpolate(x::Real, nodes, values, baryweights)

Interpolate the function represented by `values` on the `nodes` using the
corresponding barycentric weights `baryweights`.
[Kopriva, Implementing Spectral Methods for PDEs, Algorithm 31].
"""
function interpolate(x::Real, values, nodes, baryweights)
    @boundscheck begin
        @assert size(values) == size(nodes) == size(baryweights)
    end
    T = promote_type(typeof(x), eltype(values), eltype(nodes), eltype(baryweights))
    num = zero(T)
    den = zero(T)
    @inbounds for idx in eachindex(nodes)
        xval = nodes[idx]
        if xval ≈ x
            return values[idx]
        end
        t = baryweights[idx] / (x - xval)
        num += t*values[idx]
        den += t
    end
    num / den
end

function interpolate(x::AbstractVector, values, nodes, baryweights)
    ret = similar(x)
    interpolate!(ret, x, values, nodes, baryweights)
    ret
end

function interpolate!(ret::AbstractVector, x::AbstractVector, values, nodes, baryweights)
    @boundscheck begin
        @assert size(ret) == size(x)
        @assert size(values) == size(nodes) == size(baryweights)
    end
    @inbounds for idx in eachindex(ret)
        ret[idx] = interpolate(x[idx], values, nodes, baryweights)
    end
    nothing
end

function interpolate(x, values, basis::NodalBasis)
    @unpack nodes, baryweights = basis
    interpolate(x, values, nodes, baryweights)
end

function interpolate!(ret, x, values, basis::NodalBasis)
    @unpack nodes, baryweights = basis
    interpolate!(ret, x, values, nodes, baryweights)
end


"""
    interpolation_matrix(destination, source, baryweights)

Compute the matrix performing interpolation from `src` to `dest`, where
`baryweights` are the barycentric weights corresponding to `src`.
[Kopriva, Implementing Spectral Methods for PDEs, Algorithm 32].
"""
function interpolation_matrix(dest, src, baryweights)
    @boundscheck begin
        @assert length(src) == length(baryweights)
    end
    T = promote_type(eltype(dest), eltype(src), eltype(baryweights))
    mat = zeros(T, length(dest), length(src))
    @inbounds interpolation_matrix!(mat, dest, src, baryweights)
    mat
end

function interpolation_matrix!(mat, dest, src, baryweights)
    @boundscheck begin
        @assert length(src) == length(baryweights)
        @assert size(mat,1) == length(dest)
        @assert size(mat,2) == length(src)
    end
    @inbounds for k in 1:size(mat,1)
        row_has_match = false
        for j in 1:size(mat,2)
            if dest[k] ≈ src[j]
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
    nothing
end

function interpolation_matrix(dest, basis::NodalBasis)
    @unpack nodes, baryweights = basis
    interpolation_matrix(dest, nodes, baryweights)
end

function interpolation_matrix!(mat, dest, basis::NodalBasis)
    @unpack nodes, baryweights = basis
    interpolation_matrix!(mat, dest, nodes, baryweights)
end
