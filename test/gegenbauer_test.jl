using Test, PolynomialBases
import SymPy

x, α = SymPy.symbols("x, alpha")

@test 0 == SymPy.simplify( gegenbauer(x, 0, α) - 1 )
@test 0 == SymPy.simplify( gegenbauer(x, 1, α) - ( 2α*x ) )
@test 0 == SymPy.simplify( gegenbauer(x, 2, α) - ( -α + 2α*(1+α)*x^2 ) )
@test 0 == SymPy.simplify( gegenbauer(x, 3, α) - ( -2α*(1+α)*x + 4*α*(1+α)*(2+α)*x^3/3 ) )

@inferred gegenbauer(x, 4, α)
@inferred gegenbauer(10., 4, 3)
