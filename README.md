# PolynomialBases

[![Build Status](https://travis-ci.org/ranocha/PolynomialBases.jl.svg?branch=master)](https://travis-ci.org/ranocha/PolynomialBases.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/i1saoodeqrepiodl?svg=true)](https://ci.appveyor.com/project/ranocha/PolynomialBases-jl)
[![Coverage Status](https://coveralls.io/repos/ranocha/PolynomialBases.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/ranocha/PolynomialBases.jl?branch=master)
[![codecov.io](http://codecov.io/github/ranocha/PolynomialBases.jl/coverage.svg?branch=master)](http://codecov.io/github/ranocha/PolynomialBases.jl?branch=master)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/PolynomialBases.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)

A library of functions for polynomial bases used in spectral element methods using the quadrature rules from
[FastGaussQuadrature.jl](https://github.com/ajt60gaibb/FastGaussQuadrature.jl) for `Float64` and root finding
via the Newton algorithm for other scalar types (such as `BigFloat`). The algorithms for interpolation and
differentiation use barycentric weights as described in the book "Implementing Spectral Methods for PDEs"
by David Kopriva. If [SymPy.jl](https://github.com/JuliaPy/SymPy.jl)/[SymEngine.jl](https://github.com/symengine/symengine)
is loaded, symbolic computations using `SymPy.Sym`/`SymEngine.Basic` are supported.

A brief tutorial is given as
[notebook](http://nbviewer.ipython.org/github/ranocha/PolynomialBases.jl/blob/master/notebooks/Tutorial.ipynb).
