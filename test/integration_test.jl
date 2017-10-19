using Base.Test, PolynomialBases

ufunc(x) = sinpi(x)^2

function tolerance(p, T=Float64)
    if p <= 2
        1.
    elseif p <= 4
        0.65
    elseif p <= 6
        0.02
    elseif p <= 8
        1.e-4
    elseif p <= 10
        if T==Float32
            3.e-7
        else
            1.5e-7
        end
    elseif p <= 12
        if T==Float32
            3.e-7
        else
            9.e-11
        end
    elseif p <= 14
        if T==Float32
            3.e-7
        else
            3.e-14
        end
    else
        if T==Float32
            5.e-7
        else
            5.e-16
        end
    end
end

for basis_type in subtypes(PolynomialBases.NodalBasis{PolynomialBases.Line})
    for p in 3:15, T in (Float32, Float64)
        basis = basis_type(p, T)
        u = ufunc.(basis.nodes)
        @test abs(integrate(u, basis) - 1) < tolerance(p, T)
        @test abs(integrate(u, basis.weights) - 1) < tolerance(p, T)
    end
end
