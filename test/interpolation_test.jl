using Base.Test, PolynomialBases

ufunc₁(x) = sinpi(x)
ufunc₂(x) = exp(sinpi(x))
ufunc₃(x) = cospi(x^2)^5

function tolerance(p, ufunc::typeof(ufunc₁))
    if p <= 4
        1.
    elseif p <= 6
        0.03
    elseif p <= 8
        2.e-3
    elseif p <= 10
        5.e-5
    else
        2.e-7
    end
end
function tolerance(p, ufunc::typeof(ufunc₂))
    if p <= 4
        2.
    elseif p <= 6
        0.7
    elseif p <= 8
        0.2
    elseif p <= 10
        0.06
    else
        6.e-3
    end
end
function tolerance(p, ufunc::typeof(ufunc₃))
    if p <= 2
        2.
    elseif p <= 6
        0.75
    elseif p <= 8
        0.7
    elseif p <= 10
        0.4
    else
        0.06
    end
end

# some regression tests
for basis_type in subtypes(PolynomialBases.NodalBasis{PolynomialBases.Line})
    for ufc in (ufunc₁, ufunc₂, ufunc₃), p in 3:10
        basis = basis_type(p)
        u = ufc.(basis.nodes)
        xplot = linspace(-1, 1, 100)
        uplot = interpolate(xplot, u, basis)
        @test norm(ufc.(xplot) - uplot, Inf) < tolerance(p, ufc)
    end
end

# compare direct interpolation and matrix
for basis_type in subtypes(PolynomialBases.NodalBasis{PolynomialBases.Line})
    for ufc in (ufunc₁, ufunc₂, ufunc₃), p in 3:10
        basis = basis_type(p)
        u = ufc.(basis.nodes)
        xplot = linspace(-1, 1, 100)
        u1 = interpolate(xplot, u, basis)
        u2 = interpolation_matrix(xplot, basis) * u
        @test norm(u1 - u2) < 1.e-14
    end
end

# change of bases
for ufc in (ufunc₁, ufunc₂, ufunc₃), p in 3:10
    basis1 = LobattoLegendre(p)
    basis2 = GaussLegendre(p)

    u1 = ufc.(basis1.nodes)
    u2 = ufc.(basis2.nodes)

    u12 = change_basis(basis1, u2, basis2)
    u21 = change_basis(basis2, u1, basis1)

    xplot = linspace(-1, 1, 100)

    @test norm(interpolate(xplot,u1,basis1) - interpolate(xplot,u12,basis1), Inf) < tolerance(p, ufc)
    @test norm(interpolate(xplot,u2,basis2) - interpolate(xplot,u21,basis2), Inf) < tolerance(p, ufc)
end
