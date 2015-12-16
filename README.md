[![Coverage Status](https://coveralls.io/repos/alancrawford/Smolyak/badge.svg?branch=master&service=github)](https://coveralls.io/github/alancrawford/Smolyak?branch=master)
## Description

Module contains code to define Smolyak Grids and construct corresponding Interpolating Smolyak Polynomials. Both Anisotrophic and Isotrophic Grids/Polynomials are supported and are constructed efficiently following the methodology outlined in Judd, Maliar, Maliar, Valero (2014). 

The code is designed for Julia version: 0.4.

The module has 3 main components with corresponding types:

- `SmolyakGrid` : Smolyak Grid
- `SmolyakBasis` : Smolyak Basis
- `SmolyakPoly` : Smolyak Polynomial

See help functions in REPL for details.

## Example

Below is an example of how to use module to:

1. Create a Smolyak Grid
2. Define Smolyak Basis functions of an Interpolating Smolyak Polynomial on the Smolyak Grid
3. Create and evaluate the Smolyak Polynomial

```
using Smolyak
mu = [2,2,2]
lb = -2.*ones(length(mu))
ub = 3.*ones(length(mu))
sg = SmolyakGrid(mu,lb,ub)
sb = SmolyakBasis(sg)
makeBasis!(sb)
sp = SmolyakPoly(sb)
getValue!(sp,sb)
getGrad!(sp,sb)
getHess!(sp,sb)
```

Further examples can be found in help files in REPL and at [Interpolation Example](./test/Interpolation_Example.jl).