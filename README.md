# SMM

[![Build Status](https://travis-ci.com/JulienPascal/SMM.jl.svg?branch=master)](https://travis-ci.com/JulienPascal/SMM.jl)

[![Coverage Status](https://coveralls.io/repos/JulienPascal/SMM.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JulienPascal/SMM.jl?branch=master)

[![codecov.io](https://codecov.io/gh/JulienPascal/MSM.jl/branch/julia_1.5/graphs/badge.svg]
(https://codecov.io/github/<Your Organization/Acc.>/<YourRepo>?branch=master)

[![codecov.io](http://codecov.io/github/JulienPascal/SMM.jl/coverage.svg?branch=master)](http://codecov.io/github/JulienPascal/SMM.jl?branch=master)

## WARNING

This package is now deprecated and not maintained. Its (renamed) successor can be
found [here](https://github.com/JulienPascal).

---


`SMM.jl` is a package designed to facilitate the estimation of economic models
via the [Simulated Method of Moments](https://en.wikipedia.org/wiki/Method_of_simulated_moments). It relies on [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) and [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) to perform the minimization.

## Why

An economic theory can be written as a system of equations that depends on primitive
parameters. The aim of the econometrician is to recover the unknown parameters
using empirical data. As Thomas Sargent said, ["a rational expectations equilibrium is a likelihood function"](http://www.econ.ucl.ac.uk/downloads/denardi/Sargent_Interview.pdf), so a natural route to estimate the unknown parameters of a given model
is to maximize the [likelihood function](https://en.wikipedia.org/wiki/Likelihood_function), which captures how "likely" it is to obtain
the current observations for any given parameter value. However, in many
instances, this likelihood function is untractable.

An alternative approach to estimate the unknown parameters is to minimize a (weighted) distance between
the empirical [moments](https://en.wikipedia.org/wiki/Moment_(mathematics)) (means, variances, etc.) and their theoretical counterparts (which depend on paramaters). If the model is correctly specified, it is natural
to expect that the "true" parameter values are the ones generating theoretical moments close to the empirical ones. When the function mapping
the set of parameter values to the theoretical moments (the "expected response function") is known, this method is called
the [Generalized Method of Moments](https://en.wikipedia.org/wiki/Generalized_method_of_moments).
Yet, in many cases the expected response function is unknown. This issue may be circumvented by
simulating the expected response function, which is often an easy task. In this case, the method is called the [Method of Simulated Moments](https://en.wikipedia.org/wiki/Method_of_simulated_moments) or equivalently the **Simulated Method of Moments (SMM)**. This package is designed to facilitate the estimation of economic models via the SMM.

## Philosophy

`SMM.jl` is being developed with the following constraints in mind:
* the minimizing algorithm should be able to run in **parallel**, as the computational cost of simulating moments is potentially high.
* Parallelization within the function generating simulated moments is difficult
to achieve. This is generally the case when working with the simulated method of moments,
 as the time series generated are often serially correlated. This is why parallelization is done at the level of the minimization
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
* [`examples/exampleAR1Serial`](examples/exampleAR1Serial.ipynb)
* [`examples/exampleLocalToGlobal.ipynb`](examples/exampleLocalToGlobal.ipynb)
* [`examples/exampleLinearModel.ipynb`](examples/exampleLinearModel.ipynb)

## Related Packages

* [MomentOpt](https://github.com/floswald/MomentOpt.jl) : a package to do SMM using MCMC algorithms in parallel

## Acknowledgement

This work is supported by a public grant overseen by the French National Research Agency (ANR) as part of the “Investissements d’Avenir” program LIEPP (reference: ANR-11-LABX0091, ANR-11-IDEX-0005-02).
