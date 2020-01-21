using Test, PolynomialBases

if !haskey(ENV, "JULIA_PKGEVAL") # sympy is not install on https://github.com/JuliaComputing/NewPkgEval.jl
  import SymPy
  x, α = SymPy.symbols("x, alpha")

  @test 0 == SymPy.simplify( gegenbauer(x, 0, α) - 1 )
  @test 0 == SymPy.simplify( gegenbauer(x, 1, α) - ( 2α*x ) )
  @test 0 == SymPy.simplify( gegenbauer(x, 2, α) - ( -α + 2α*(1+α)*x^2 ) )
  @test 0 == SymPy.simplify( gegenbauer(x, 3, α) - ( -2α*(1+α)*x + 4*α*(1+α)*(2+α)*x^3/3 ) )

  @test_skip gegenbauer(x, 4, α)
end
@inferred gegenbauer(10., 4, 3)
