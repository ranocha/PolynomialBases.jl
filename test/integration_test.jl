using Base.Test, PolynomialBases

ufunc(x) = sinpi(x)^2

function tolerance(p)
    if p <= 2
        1.
    elseif p <= 4
        0.65
    elseif p <= 6
        0.02
    elseif p <= 8
        1.e-4
    elseif p <= 10
        1.5e-7
    elseif p <= 12
        9.e-11
    elseif p <= 14
        3.e-14
    else
        5.e-16
    end
end

for basis_type in subtypes(PolynomialBases.NodalBasis{PolynomialBases.Line})
    for p in 3:15
        basis = basis_type(p)
        u = ufunc.(basis.nodes)
        @test abs(integrate(u, basis) - 1) < tolerance(p)
    end
end
