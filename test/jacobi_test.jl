using Test, PolynomialBases, FastGaussQuadrature

if !haskey(ENV, "JULIA_PKGEVAL") # sympy is not install on https://github.com/JuliaComputing/NewPkgEval.jl
    import SymPy
    x, α, β = SymPy.symbols("x, alpha, beta")

    @test 0 == SymPy.simplify( jacobi(x, 0, α, β) - 1 )
    # https://www.wolframalpha.com/input/?i=JacobiP%5B1,+alpha,+beta,+x%5D
    @test 0 == SymPy.simplify( jacobi(x, 1, α, β) - ( α-β+x*(2+α+β) ) / 2 )
    # https://www.wolframalpha.com/input/?i=JacobiP%5B2,+alpha,+beta,+x%5D
    @test 0 == SymPy.simplify( jacobi(x, 2, α, β) - (
                (α+1)*(α+2)/2 + (x-1)^2*(α+β+3)*(α+β+4)/8
                + (α+2)*(x-1)*(α+β+3)/2
            ) )
    # https://www.wolframalpha.com/input/?i=JacobiP%5B3,+alpha,+beta,+x%5D
    @test 0 == SymPy.simplify( jacobi(x, 3, α, β) - (
                (α+1)*(α+2)*(α+3)/6 + (x-1)^3*(α+β+4)*(α+β+5)*(α+β+6)/48
                + (α+3)*(x-1)^2*(α+β+4)*(α+β+5)/8 + (α+2)*(α+3)*(x-1)*(α+β+4)/4
            ) )

    @test_skip jacobi(x, 4, α, β)
end
@inferred jacobi(10., 4, 3, 2)

# Gauss Jacobi nodes and weights
for p in 0:10, α in range(-0.9, stop=4, length=50), β in range(-0.9, stop=4, length=50)
    x1, w1 = gaussjacobi(p+1, α, β)
    @inferred PolynomialBases.gauss_jacobi_nodes_and_weights(p, α, β)
    x2, w2 = PolynomialBases.gauss_jacobi_nodes_and_weights(p, α, β)
    @test x1 ≈ x2 atol=1.e-11
    @test w1 ≈ w2 atol=1.e-11
end

# Vandermonde matrices for Gauss Jacobi nodes and weights
# NOTE: Only some tests for low polynomial degrees since the Vandermonde matrix
#       becomes ill-conditioned for higher values of p.
for p in 0:8, α in range(-0.9, stop=4, length=50), β in range(-0.9, stop=4, length=50)
    basis1 = GaussJacobi(p, α, β)
    V1 = jacobi_vandermonde(basis1, α, β)

    basis2 = GaussLegendre(p+1)
    V2 = jacobi_vandermonde(basis2, α, β)

    R = interpolation_matrix(basis1.nodes, basis2)

    P = Matrix(I, p+2, p+1)

    @test R * V2 * P / V1 ≈ I
end
