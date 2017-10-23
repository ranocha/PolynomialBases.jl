using Base.Test, PolynomialBases, FastGaussQuadrature
import SymPy

x = SymPy.symbols("x")

@test 0 == SymPy.simplify( legendre(x, 0) - 1 )
@test 0 == SymPy.simplify( legendre(x, 1) - x )
@test 0 == SymPy.simplify( legendre(x, 2) - ( 3x^2 - 1 ) / 2 )
@test 0 == SymPy.simplify( legendre(x, 3) - ( 5x^3 - 3x) / 2 )
@test 0 == SymPy.simplify( legendre(x, 4) - ( 35x^4 - 30x^2 + 3 ) / 8 )
@test 0 == SymPy.simplify( legendre(x, 5) - ( 63x^5 - 70x^3 + 15x ) / 8 )
@test 0 == SymPy.simplify( legendre(x, 6) - ( 231x^6 - 315x^4 + 105x^2 - 5 ) / 16 )

@inferred legendre(x, 4)
@inferred legendre(10., 4)

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

# Derivative and mass matrices for Gauss Legendre nodes and weights
# NOTE: Only some tests for low polynomial degrees since the Vandermonde matrix
#       becomes very ill-conditioned for higher values of p.
for p in 0:8
    basis = GaussLegendre(p)
    V = legendre_vandermonde(basis)

    Dhat = legendre_D(p)
    @test norm(basis.D - V * Dhat / V) < 1.e-13

    Mhat = legendre_M(p)
    @test norm(Diagonal(basis.weights) - V' \ Mhat / V) < 1.e-13
end
