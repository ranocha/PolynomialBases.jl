using Base.Test, PolynomialBases, FastGaussQuadrature

# Gauss Legendre nodes and weights
for p in 0:20
    x1, w1 = gausslegendre(p+1)
    @inferred PolynomialBases.gauss_legendre_nodes_and_weights(p)
    x2, w2 = PolynomialBases.gauss_legendre_nodes_and_weights(p)
    @test norm(x1 - x2) < 1.e-14
    @test norm(w1 - w2) < 1.e-14
end

# Lobatto Legendre nodes and weights
for p in 1:20
    x1, w1 = gausslobatto(p+1)
    @inferred PolynomialBases.lobatto_legendre_nodes_and_weights(p)
    x2, w2 = PolynomialBases.lobatto_legendre_nodes_and_weights(p)
    @test norm(x1 - x2) < 1.e-14
    @test norm(w1 - w2) < 1.e-14
end

# derivative and mass matrices for Gauss Legendre nodes and weights
for p in 0:11
    basis = GaussLegendre(p)
    V = legendre_vandermonde(basis)

    Dhat = legendre_D(p)
    @test norm(basis.D - V * Dhat / V) < 1.e-13

    Mhat = legendre_M(p)
    @test norm(Diagonal(basis.weights) - V' \ Mhat / V) < 1.e-13
end
