using Test, PolynomialBases

if !haskey(ENV, "JULIA_PKGEVAL") # sympy is not install on https://github.com/JuliaComputing/NewPkgEval.jl
  import SymPy
  x = SymPy.symbols("x")

  @test 0 == SymPy.simplify( hermite(x, 0) - 1 )
  @test 0 == SymPy.simplify( hermite(x, 1) - ( 2x ) )
  @test 0 == SymPy.simplify( hermite(x, 2) - ( 4x^2 - 2 ) )
  @test 0 == SymPy.simplify( hermite(x, 3) - ( 8x^3 - 12x ) )
  @test 0 == SymPy.simplify( hermite(x, 4) - ( 16x^4 - 48x^2 + 12 ) )

  @test_skip hermite(x, 4)
end
@inferred hermite(10., 4)
