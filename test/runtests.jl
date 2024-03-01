using Test
using PolynomialBases

@elapsed begin
    @time @testset "Explicit Imports" begin include("explicit_imports.jl") end
    @time @testset "Canonical Mappings" begin include("canonical_mappings_test.jl") end
    @time @testset "Interpolation" begin include("interpolation_test.jl") end
    @time @testset "Integration" begin include("integration_test.jl") end
    @time @testset "Derivatives" begin include("derivative_test.jl") end
    @time @testset "Legendre" begin include("legendre_test.jl") end
    @time @testset "Gegenbauer" begin include("gegenbauer_test.jl") end
    @time @testset "Jacobi" begin include("jacobi_test.jl") end
    @time @testset "Hermite" begin include("hermite_test.jl") end
    @time @testset "Laguerre" begin include("laguerre_test.jl") end
    @time @testset "Hahn" begin include("hahn_test.jl") end
    @time @testset "Utilities" begin include("utilities_test.jl") end
    @time @testset "Symbolic Bases (SymPy)" begin include("sympy_test.jl") end
    @time @testset "Symbolic Bases (SymEngine)" begin include("symengine_test.jl") end
end
