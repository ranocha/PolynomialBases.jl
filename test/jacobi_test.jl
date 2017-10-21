using Base.Test, PolynomialBases, FastGaussQuadrature

# Gauss Jacobi nodes and weights
for p in 0:10, α in linspace(-0.9, 4, 50), β in linspace(-0.9, 4, 50)
    x1, w1 = gaussjacobi(p+1, α, β)
    @inferred PolynomialBases.gauss_jacobi_nodes_and_weights(p, α, β)
    x2, w2 = PolynomialBases.gauss_jacobi_nodes_and_weights(p, α, β)
    @test norm(x1 - x2) < 1.e-11
    @test norm(w1 - w2) < 1.e-11
end

# Vandermonde matrices for Gauss Jacobi nodes and weights
# NOTE: Only some tests for low polynomial degrees since the Vandermonde matrix
#       becomes very ill-conditioned for higher values of p.
for p in 0:8, α in linspace(-0.9, 4, 50), β in linspace(-0.9, 4, 50)
    basis1 = GaussJacobi(p, α, β)
    V1 = jacobi_vandermonde(basis1, α, β)

    basis2 = GaussLegendre(p+1)
    V2 = jacobi_vandermonde(basis2, α, β)

    R = interpolation_matrix(basis1.nodes, basis2)

    P = eye(p+2, p+1)

    @test norm(R * V2 * P / V1 - I) < 2.e-13
end
