using Base.Test, PolynomialBases

# regression test
for basis_type in subtypes(PolynomialBases.NodalBasis{PolynomialBases.Line})
    println(DevNull, "  ", basis_type(1))
    for p in 0:20, T in (Float32, Float64)
        basis = basis_type(p, T)
        xmin, xmax = 0.25, 0.9
        x = similar(basis.nodes)
        ξ = similar(x)
        map_from_canonical!(x, basis.nodes, xmin, xmax, basis)
        map_to_canonical!(ξ, x, xmin, xmax, basis)
        @test norm(ξ - basis.nodes) < 10*eps(T)
    end
end
