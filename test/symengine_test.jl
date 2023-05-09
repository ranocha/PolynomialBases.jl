using Test, PolynomialBases
import SymEngine

tol = 5.e-14

for p in 0:6
    basis_symbolic = LobattoLegendre(p, SymEngine.Basic)
    basis_float    = LobattoLegendre(p, Float64)
    @test maximum(abs.( float.(basis_symbolic.nodes) - basis_float.nodes )) < tol
    @test maximum(abs.( float.(basis_symbolic.weights) - basis_float.weights )) < tol
    @test maximum(abs.( float.(basis_symbolic.baryweights) - basis_float.baryweights )) < tol
    @test maximum(abs.( float.(basis_symbolic.D) - basis_float.D, )) < tol
end
@test_throws ArgumentError LobattoLegendre(7, SymEngine.Basic)

for p in 0:4
    basis_symbolic = GaussLegendre(p, SymEngine.Basic)
    basis_float = GaussLegendre(p, Float64)
    @test maximum(abs.( float.(basis_symbolic.nodes) - basis_float.nodes )) < tol
    @test maximum(abs.( float.(basis_symbolic.weights) - basis_float.weights )) < tol
    @test maximum(abs.( float.(basis_symbolic.baryweights) - basis_float.baryweights )) < tol
    @test maximum(abs.( float.(basis_symbolic.D) - basis_float.D, )) < tol
end
@test_throws ArgumentError GaussLegendre(5, SymEngine.Basic)

interpolation_matrix([-1, 1], LobattoLegendre(4, SymEngine.Basic))
interpolation_matrix([-1, 1], GaussLegendre(4, SymEngine.Basic))
