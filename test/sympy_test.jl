using Test, PolynomialBases

if !haskey(ENV, "JULIA_PKGEVAL") # sympy is not installed on https://github.com/JuliaComputing/NewPkgEval.jl
  import SymPy

  tol = 5.e-15

  for p in 0:6
      basis_sympy = LobattoLegendre(p, SymPy.Sym)
      basis_float = LobattoLegendre(p, Float64)
      @test maximum(abs.( float.(basis_sympy.nodes) - basis_float.nodes )) < tol
      @test maximum(abs.( float.(basis_sympy.weights) - basis_float.weights )) < tol
      @test maximum(abs.( float.(basis_sympy.baryweights) - basis_float.baryweights )) < tol
      @test maximum(abs.( float.(basis_sympy.D) - basis_float.D, )) < 2tol
  end
  @test_throws ArgumentError LobattoLegendre(7, SymPy.Sym)

  for p in 0:4
      basis_sympy = GaussLegendre(p, SymPy.Sym)
      basis_float = GaussLegendre(p, Float64)
      @test maximum(abs.( float.(basis_sympy.nodes) - basis_float.nodes )) < tol
      @test maximum(abs.( float.(basis_sympy.weights) - basis_float.weights )) < tol
      @test maximum(abs.( float.(basis_sympy.baryweights) - basis_float.baryweights )) < tol
      @test maximum(abs.( float.(basis_sympy.D) - basis_float.D, )) < 2tol
  end
  @test_throws ArgumentError GaussLegendre(5, SymPy.Sym)

  interpolation_matrix([-1, 1], LobattoLegendre(4, SymPy.Sym))
  interpolation_matrix([-1, 1], GaussLegendre(4, SymPy.Sym))
end
