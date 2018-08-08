using Test, PolynomialBases
import SymPy

x = SymPy.symbols("x")

@test 0 == SymPy.simplify( hermite(x, 0) - 1 )
@test 0 == SymPy.simplify( hermite(x, 1) - ( 2x ) )
@test 0 == SymPy.simplify( hermite(x, 2) - ( 4x^2 - 2 ) )
@test 0 == SymPy.simplify( hermite(x, 3) - ( 8x^3 - 12x ) )
@test 0 == SymPy.simplify( hermite(x, 4) - ( 16x^4 - 48x^2 + 12 ) )

@inferred hermite(x, 4)
@inferred hermite(10., 4)
