using Base.Test, PolynomialBases

ufunc(x) = sinpi(x)

for basis_type in subtypes(PolynomialBases.NodalBasis{PolynomialBases.Line})
    for p in 1:5
        basis = basis_type(p)
        u = ufunc.(basis.nodes)
        D, M, R, B, MinvRtB = utility_matrices(basis)

        @test integrate(u, basis) ≈ sum(M * u)
        Ru = R*u
        @test Ru[1] ≈ interpolate(-1, u, basis)
        @test Ru[2] ≈ interpolate(+1, u, basis)

        @inferred includes_boundaries(basis)
        if includes_boundaries(basis) == Val{true}()
            @test Ru[1] ≈ u[1]
            @test Ru[end] ≈ u[end]
        else
            @test !(Ru[1] ≈ u[1])
            @test !(Ru[end] ≈ u[end])
        end

        basis_type <: ClosedNewtonCotes && continue

        # SBP property
        @test norm(M*D + D'*M - R'*B*R) < 1.e-14
    end
end
