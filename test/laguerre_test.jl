using Base.Test, PolynomialBases
import SymPy

x = SymPy.symbols("x")

@test 0 == SymPy.simplify( laguerre(x, 0) - 1 )
@test 0 == SymPy.simplify( laguerre(x, 1) - ( 1-x ) )
@test 0 == SymPy.simplify( laguerre(x, 2) - ( x^2 - 4x + 2 ) / 2 )
@test 0 == SymPy.simplify( laguerre(x, 3) - ( -x^3 + 9x^2 - 18x + 6 ) / 6 )
@test 0 == SymPy.simplify( laguerre(x, 4) - ( x^4 - 16x^3 + 72x^2 - 96x + 24 ) / 24 )
@test 0 == SymPy.simplify( laguerre(x, 5) - ( -x^5 + 25x^4 - 200x^3 + 600x^2 - 600x + 120 ) / 120 )
@test 0 == SymPy.simplify( laguerre(x, 6) - ( x^6 - 36x^5 + 450x^4 - 2400x^3 + 5400x^2 - 4320x + 720 ) / 720 )

@inferred laguerre(x, 4)
@inferred laguerre(10., 4)
