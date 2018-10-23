# SMM

[![Build Status](https://travis-ci.com/JulienPascal/SMM.jl.svg?branch=master)](https://travis-ci.com/JulienPascal/SMM.jl)

[![Coverage Status](https://coveralls.io/repos/JulienPascal/SMM.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JulienPascal/SMM.jl?branch=master)

[![codecov.io](http://codecov.io/github/JulienPascal/SMM.jl/coverage.svg?branch=master)](http://codecov.io/github/JulienPascal/SMM.jl?branch=master)

`SMM.jl` is a package designed to facilitate the estimation of economic models
via the [Simulated Method of Moments](https://en.wikipedia.org/wiki/Method_of_simulated_moments). It relies on [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) to perform the minimization.


## Philosophy

`SMM.jl` is being developed with the following constraints in mind:
* the minimizing algorithm should be able to run in **parallel**, as the computational cost of simulating moments, for a given parameter value, is potentially high.
* Parallelization within the function generating simulated moments is difficult
to achieve. This is generally the case when working with the simulated method of moments,
 as the time series generated are serially correlated. This is why parallelization is done at the level of the minimization
algorithm itself.
* The minimizing algorithm should search for a **global** minimum, as the
objective function may have multiple local minima.


## Installation

This package is still in its development phase. Yet, If you feel brave enough:
```
Pkg.clone("https://github.com/JulienPascal/SMM.jl.git")
```

## Usage

See the following notebooks:
* [`examples/example2dNormalParallel.ipynb`](examples/example2dNormalParallel.ipynb)
* [`examples/example2dNormalSerial.ipynb`](examples/example2dNormalSerial.ipynb)
* [`examples/exampleAR1Serial`](examples/exampleAR1Serial)

## Related Packages

* [MomentOpt](https://github.com/floswald/MomentOpt.jl) : a package to do SMM using MCMC algorithms in parallel
