using Base.Test, PolynomialBases

ufunc(x) = sinpi(x)
uprim(x) = π*cospi(x)

function tolerance(p)
    if p <= 2
        5.
    elseif p <= 6
        0.5
    elseif p <= 8
        0.04
    elseif p <= 10
        0.002
    elseif p <= 12
        3.e-5
    elseif p <= 14
        6.e-7
    elseif p <= 16
        7.e-9
    else
        7.e-11
    end
end

# regression test
for basis_type in subtypes(PolynomialBases.NodalBasis{PolynomialBases.Line})
    println("  ", basis_type(1))
    for p in 5:20
        basis = basis_type(p)
        u = ufunc.(basis.nodes)
        @test norm(basis.D * u - uprim.(basis.nodes)) < tolerance(p)
    end
end

# compare direct derivative evaluation and matrix
for basis_type in subtypes(PolynomialBases.NodalBasis{PolynomialBases.Line})
    for p in 0:15
        basis = basis_type(p)
        u = ufunc.(basis.nodes)
        xplot = linspace(-1, 1, 100)
        u1 = derivative_at(xplot, u, basis)
        u2 = interpolate(xplot, basis.D*u, basis)
        @test norm(u1 - u2, Inf) < 5.e-12
    end
end
