
doc"
    jacobi_nodes_and_weights{T<:Real}(degree::Int, α::T, β::T)

Compute the `degree-1` nodes and weights of the Jacobi Gauss quadrature, i.e.
a polynomial of degree `degree` can be represented exactly on the nodes. The
weights of the Jacobi polynomial are `α,β`.
"
function jacobi_gauss_nodes_and_weights{T<:Real}(degree::Int, α::T, β::T)
    if degree < 0
        DomainError("Polynomial degree must be non-negative (`degree` = $degree).")
    end

    if degree == 0
        nodes = [(α-β)/(α+β+2)]
        weights = T[2]
    else
        h = 2 .* (0:degree) .+ (α+β)
        hh = h[1:degree]
        dv = @. -(α^2-β^2) / ((h+2) * h)
        if α+β < 10*eps(T)
            dv[1] = 0
        end
        ev = @. 2/(hh+2)*sqrt((1:degree)*((1:degree)+α+β)*((1:degree)+α)*
                                ((1:degree)+β)/((hh+1)*(hh+3)))
        mat = SymTridiagonal(dv, ev)
        nodes, V = eig(mat)
        weights = V[1,:].^2 * 2^(α+β+1)/(α+β+1)*gamma(α+1)*gamma(β+1)/gamma(α+β+1)
    end

    nodes, weights
end

function jacobi_gauss_nodes_and_weights(degree::Int, α::Real, β::Real)
    jacobi_gauss_nodes_and_weights(degree, promote(α,β)...)
end

function jacobi_gauss_nodes_and_weights(degree::Int, α::Real)
    jacobi_gauss_nodes_and_weights(degree, α, zero(α))
end

function jacobi_gauss_nodes_and_weights(degree::Int, T=Float64)
    jacobi_gauss_nodes_and_weights(degree, zero(T), zero(T))
end


doc"
    jacobi_lobatto_nodes_and_weights{T<:Real}(degree::Int, α::T, β::T)

Compute the `degree-1` nodes and weights of the Jacobi Lobatto quadrature, i.e.
a polynomial of degree `degree` can be represented exactly on the nodes. The
weights of the Jacobi polynomial are `α,β`.
"
function jacobi_lobatto_nodes_and_weights{T<:Real}(degree::Int, α::T, β::T)
    if degree < 0
        DomainError("Polynomial degree must be non-negative (`degree` = $degree).")
    end

    if degree == 0
        nodes = [(α-β)/(α+β+2)]
        weights = T[2]
    elseif degree == 1
        nodes = T[-1, 1]
        weights = [T(1)/2, T(1)/2]
    else
        interior_nodes, interior_weights = jacobi_gauss_nodes_and_weights(degree-2, α+1, β+1)
        nodes = vcat(-one(T), interior_nodes, one(T))

    end

    nodes, weights
end
