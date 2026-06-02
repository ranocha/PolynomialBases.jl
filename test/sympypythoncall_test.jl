using Test, PolynomialBases

# PyCall and PythonCall cannot share the same Python interpreter if PyCall
# initialized it first; skip these tests when PyCall is already loaded.
pycall_loaded = any(k -> k.name == "PyCall", keys(Base.loaded_modules))

if !haskey(ENV, "JULIA_PKGEVAL") && !pycall_loaded
  import SymPyPythonCall

  tol = 5.e-15

  for p in 0:6
      basis_sympypythoncall = LobattoLegendre(p, SymPyPythonCall.Sym)
      basis_float = LobattoLegendre(p, Float64)
      @test maximum(abs.( float.(basis_sympypythoncall.nodes) - basis_float.nodes )) < tol
      @test maximum(abs.( float.(basis_sympypythoncall.weights) - basis_float.weights )) < tol
      @test maximum(abs.( float.(basis_sympypythoncall.baryweights) - basis_float.baryweights )) < tol
      @test maximum(abs.( float.(basis_sympypythoncall.D) - basis_float.D, )) < 2tol
  end
  @test_throws ArgumentError LobattoLegendre(7, SymPyPythonCall.Sym)

  for p in 0:4
      basis_sympypythoncall = GaussLegendre(p, SymPyPythonCall.Sym)
      basis_float = GaussLegendre(p, Float64)
      @test maximum(abs.( float.(basis_sympypythoncall.nodes) - basis_float.nodes )) < tol
      @test maximum(abs.( float.(basis_sympypythoncall.weights) - basis_float.weights )) < tol
      @test maximum(abs.( float.(basis_sympypythoncall.baryweights) - basis_float.baryweights )) < tol
      @test maximum(abs.( float.(basis_sympypythoncall.D) - basis_float.D, )) < 2tol
  end
  @test_throws ArgumentError GaussLegendre(5, SymPyPythonCall.Sym)

  interpolation_matrix([-1, 1], LobattoLegendre(4, SymPyPythonCall.Sym))
  interpolation_matrix([-1, 1], GaussLegendre(4, SymPyPythonCall.Sym))
end
