using Test, PolynomialBases, PolynomialBases.FastGaussQuadrature

if !haskey(ENV, "JULIA_PKGEVAL") # sympy is not installed on https://github.com/JuliaComputing/NewPkgEval.jl
  import SymPy
  x = SymPy.symbols("x")

  @test 0 == SymPy.simplify( legendre(x, 0) - 1 )
  @test 0 == SymPy.simplify( legendre(x, 1) - x )
  @test 0 == SymPy.simplify( legendre(x, 2) - ( 3x^2 - 1 ) / 2 )
  @test 0 == SymPy.simplify( legendre(x, 3) - ( 5x^3 - 3x) / 2 )
  @test 0 == SymPy.simplify( legendre(x, 4) - ( 35x^4 - 30x^2 + 3 ) / 8 )
  @test 0 == SymPy.simplify( legendre(x, 5) - ( 63x^5 - 70x^3 + 15x ) / 8 )
  @test 0 == SymPy.simplify( legendre(x, 6) - ( 231x^6 - 315x^4 + 105x^2 - 5 ) / 16 )

  @test_skip legendre(x, 4)
end
@inferred legendre(10., 4)

# Gauss Legendre nodes and weights
for T in (Float32, Float64, BigFloat)
    for p in 0:20
        x1, w1 = gausslegendre(T, p+1)
        x2, w2 = @inferred PolynomialBases.gauss_legendre_nodes_and_weights(p, T)
        @test x1 ≈ x2
        @test w1 ≈ w2
    end
end

# Lobatto Legendre nodes and weights
for T in (Float32, Float64, BigFloat)
    for p in 1:20
        x1, w1 = gausslobatto(T, p+1)
        x2, w2 = @inferred PolynomialBases.lobatto_legendre_nodes_and_weights(p, T)
        @test x1 ≈ x2
        @test w1 ≈ w2
    end
end

# Derivative and mass matrices for Gauss Legendre nodes and weights
# NOTE: Only some tests for low polynomial degrees since the Vandermonde matrix
#       becomes very ill-conditioned for higher values of p.
for p in 0:8
    basis = GaussLegendre(p)
    V = legendre_vandermonde(basis)

    Dhat = legendre_D(p)
    @test derivative_matrix(basis) ≈ V*Dhat/V atol=1.e-13

    Mhat = legendre_M(p)
    @test mass_matrix(basis) ≈ V'\Mhat/V atol=1.e-13
end
