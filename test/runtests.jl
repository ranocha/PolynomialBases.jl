using Base.Test
using PolynomialBases

tic()

@time @testset "Interpolation" begin include("interpolation_test.jl") end
@time @testset "Integration" begin include("integration_test.jl") end
@time @testset "Derivatives" begin include("derivative_test.jl") end
@time @testset "Utilities" begin include("utilities_test.jl") end

toc()
