using Test, PolynomialBases

ufunc(x) = sinpi(x)

for basis_type in subtypes(PolynomialBases.NodalBasis{PolynomialBases.Line})
    for p in 1:5, T in (Float32, Float64)
        basis = basis_type(p, T)
        @inferred basis_type(p, T)
        @test eltype(basis.nodes) == T
        @test eltype(basis) == T
        @test basis == basis_type(p, T)
        u = ufunc.(basis.nodes)
        D, M, R, B, MinvRtB = utility_matrices(basis)
        @test Matrix(D) == D

        @test integrate(u, basis) ≈ sum(M * u)
        Ru = R*u
        @test Ru[1] ≈ interpolate(-1, u, basis)
        @test Ru[2] ≈ interpolate(+1, u, basis)

        @inferred includes_boundaries(basis)
        @inferred includes_left_boundary(basis)
        @inferred includes_right_boundary(basis)
        if includes_left_boundary(basis) == Val{true}()
            @test Ru[1] ≈ u[1]
        else
            @test !(Ru[1] ≈ u[1])
        end
        if includes_right_boundary(basis) == Val{true}()
            @test Ru[end] ≈ u[end]
        else
            @test !(Ru[end] ≈ u[end])
        end

        # SBP property
        if satisfies_sbp(basis) == Val{true}()
            @test M*D + D'*M ≈ mass_matrix_boundary(basis)
        end
    end
end
