using Base.Test, PolynomialBases
import SymPy

x, α = SymPy.symbols("x, alpha")

@test 0 == SymPy.simplify( laguerre(x, 0) - 1 )
@test 0 == SymPy.simplify( laguerre(x, 1) - ( 1-x ) )
@test 0 == SymPy.simplify( laguerre(x, 2) - ( x^2 - 4x + 2 ) / 2 )
@test 0 == SymPy.simplify( laguerre(x, 3) - ( -x^3 + 9x^2 - 18x + 6 ) / 6 )
@test 0 == SymPy.simplify( laguerre(x, 4) - ( x^4 - 16x^3 + 72x^2 - 96x + 24 ) / 24 )
@test 0 == SymPy.simplify( laguerre(x, 5) - ( -x^5 + 25x^4 - 200x^3 + 600x^2 - 600x + 120 ) / 120 )
@test 0 == SymPy.simplify( laguerre(x, 6) - ( x^6 - 36x^5 + 450x^4 - 2400x^3 + 5400x^2 - 4320x + 720 ) / 720 )

@test 0 == SymPy.simplify( laguerre(x, 0, α) - 1 )
# https://www.wolframalpha.com/input/?i=LaguerreL%5B1,+alpha,+x%5D
@test 0 == SymPy.simplify( laguerre(x, 1, α) - ( α-x+1 ) )
# https://www.wolframalpha.com/input/?i=LaguerreL%5B2,+alpha,+x%5D
@test 0 == SymPy.simplify( laguerre(x, 2, α) - ( α^2 + 3α + x^2 - 2α*x - 4x + 2 ) / 2 )
# https://www.wolframalpha.com/input/?i=LaguerreL%5B3,+alpha,+x%5D
@test 0 == SymPy.simplify( laguerre(x, 3, α) - ( α^3+6α^2+11α-x^3+3α*x^2+9x^2-3α^2*x-15α*x-18x+6 ) / 6 )
# https://www.wolframalpha.com/input/?i=LaguerreL%5B4,+alpha,+x%5D
@test 0 == SymPy.simplify( laguerre(x, 4, α) - ( α^4+10α^3+35α^2+50α+x^4-4α*x^3-16x^3+6α^2*x^2+42α*x^2+72x^2-4α^3*x-36α^2*x-104α*x-96x+24 ) / 24 )

@inferred laguerre(x, 4)
@inferred laguerre(x, 4, α)
@inferred laguerre(x, 4, 2)
@inferred laguerre(10., 4)
@inferred laguerre(10., 4, 2.)
@inferred laguerre(10., 4, 2)
