using Test
using PolynomialBases
using ExplicitImports

@test check_no_implicit_imports(PolynomialBases) === nothing

@test check_no_stale_explicit_imports(PolynomialBases) === nothing
