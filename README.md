# PolynomialBases

[![Build Status](https://github.com/ranocha/PolynomialBases.jl/workflows/CI/badge.svg)](https://github.com/ranocha/PolynomialBases.jl/actions)
[![Codecov](http://codecov.io/github/ranocha/PolynomialBases.jl/coverage.svg?branch=master)](http://codecov.io/github/ranocha/PolynomialBases.jl?branch=master)
[![Coveralls](https://coveralls.io/repos/ranocha/PolynomialBases.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/ranocha/PolynomialBases.jl?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![GitHub commits since tagged version](https://img.shields.io/github/commits-since/ranocha/PolynomialBases.jl/v0.4.8.svg?style=social&logo=github)](https://github.com/ranocha/PolynomialBases.jl)
<!-- [![Build status](https://ci.appveyor.com/api/projects/status/i1saoodeqrepiodl?svg=true)](https://ci.appveyor.com/project/ranocha/PolynomialBases-jl)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/PolynomialBases.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html) -->

A library of functions for polynomial bases used in spectral element methods using the quadrature rules from
[FastGaussQuadrature.jl](https://github.com/ajt60gaibb/FastGaussQuadrature.jl) for `Float64` and root finding
via the Newton algorithm for other scalar types (such as `BigFloat`). The algorithms for interpolation and
differentiation use barycentric weights as described in the book "Implementing Spectral Methods for PDEs"
by David Kopriva. If [SymPy.jl](https://github.com/JuliaPy/SymPy.jl)/[SymEngine.jl](https://github.com/symengine/symengine)
is loaded, symbolic computations using `SymPy.Sym`/`SymEngine.Basic` are supported.

A brief tutorial is given as
[notebook](http://nbviewer.ipython.org/github/ranocha/PolynomialBases.jl/blob/master/notebooks/Tutorial.ipynb).
