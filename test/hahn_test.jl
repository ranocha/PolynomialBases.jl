using Test, PolynomialBases

if !haskey(ENV, "JULIA_PKGEVAL") # sympy is not installed on https://github.com/JuliaComputing/NewPkgEval.jl
    import SymPy
    x, α, β, N = SymPy.symbols("x, alpha, beta, N")

    @test 0 == SymPy.simplify( hahn(x, 0, α, β, N) - 1 )
    # https://www.wolframalpha.com/input/?i=HypergeometricPFQ%5B+%7B-1,+1%2Balpha%2Bbeta%2B1,+-x%7D,+%7Balpha%2B1,+-N%7D,+1%5D
    @test 0 == SymPy.simplify( hahn(x, 1, α, β, N) - ( 1 - (α+β+2)*x/(N*(α+1)) ) )
    # https://www.wolframalpha.com/input/?i=HypergeometricPFQ%5B+%7B-2,+2%2Balpha%2Bbeta%2B1,+-x%7D,+%7Balpha%2B1,+-N%7D,+1%5D
    @test 0 == SymPy.simplify( hahn(x, 2, α, β, N) - (
            (1-x)*x*(α+β+4)*(α+β+3) / ((α+1)*(α+2)*(1-N)*N)
            - 2x*(α+β+3)/(N*(α+1)) + 1
        ) )
    # https://www.wolframalpha.com/input/?i=HypergeometricPFQ%5B+%7B-3,+3%2Balpha%2Bbeta%2B1,+-x%7D,+%7Balpha%2B1,+-N%7D,+1%5D
    @test 0 == SymPy.simplify( hahn(x, 3, α, β, N) - (
            3*(1-x)*x*(α+β+5)*(α+β+4) / ( (α+1)*(α+2)*N*(1-N) )
            - (1-x)*(2-x)*x*(α+β+5)*(α+β+6)*(α+β+4) / ( (α+1)*(α+2)*(α+3)*N*(1-N)*(2-N) )
            - 3x*(α+β+4) / ( N*(α+1) ) + 1
        ) )
    # https://www.wolframalpha.com/input/?i=HypergeometricPFQ%5B+%7B-4,+4%2Balpha%2Bbeta%2B1,+-x%7D,+%7Balpha%2B1,+-N%7D,+1%5D
    @test 0 == SymPy.simplify( hahn(x, 4, α, β, N) - (
            6*(1-x)*x*(α+β+6)*(α+β+5) / ( (α+1)*(α+2)*N*(1-N) )
            - 4*(1-x)*(2-x)*x*(α+β+6)*(α+β+7)*(α+β+5) / ( (α+1)*(α+2)*(α+3)*N*(1-N)*(2-N) )
            + (1-x)*(2-x)*(3-x)*x*(α+β+6)*(α+β+7)*(α+β+8)*(α+β+5) / ( (α+1)*(α+2)*(α+3)*(α+4)*N*(1-N)*(2-N)*(3-N) )
            - 4x*(α+β+5) / ( N*(α+1) )
            + 1
        ) )

    @test_throws ArgumentError hahn(x, 6, α, β, 5)
    @test_skip hahn(x, 4, α, β, 5)
end
@inferred hahn(10., 4, -0.5, -0.9, 5)
