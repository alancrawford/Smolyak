[![Coverage Status](https://coveralls.io/repos/alancrawford/Smolyak/badge.svg?branch=master&service=github)](https://coveralls.io/github/alancrawford/Smolyak?branch=master)
## Description

This module contains code to define Smolyak polynomials, basis functions and grids. 

The Smolyak algorithm provides and efficient way to construct multivariate function approximations. The accuracy of the approximation in each dimension of the state vector is linked to the 'level'. In an isotrophic Smolyak polynomial the level is the same for all dimensions of the state vector. Anisotrophic Smolyak polynomials allow the user to vary the accuracy of the interpolating polynomial in each dimension by specifiying a dimension specific level. 

Crudely put, a higher level of accuracy corresponds to interpolating higher order basis functions in each dimension. The cost of this is computational expense of computing more tensor products in the approximating polynomial. 

However, using the Smolyak algorithm ensures that the number of tensor products increase at a polynomial - rather than exponential - rate linked the level of accuracy chosen by the user. Therefore they are very effective at mitigating the curse of dimensionality.

Both Anisotrophic and Isotrophic Grids/Polynomials are supported and are constructed efficiently following the methodology outlined in Judd, Maliar, Maliar, Valero (2014). 

Moreover, the package allows the user to choose the univariate basis function type: ordinary, chebyshev and an augmented version of spread polynomials are available. 

Chebyshev polynomials have many benefits over ordinary polynomials. However, unlike ordinary polynomial they are not - unlike ordinary polynomials - inherently sparse. This is unlikely to be costly in many applications. However, when sparsity is important the user might like to use spread polynomials. 

The augmented spread polynomials are a translated version of chebyshev polynomials. Therefore they share some key benefits of Chebyshev polynomials (i.e. low correlation of basis functions), and they are - like ordinary polynomials - inherently sparse. 

The code is designed for Julia version: 0.7/ 1.0.

The module is designed around a `SmolyakKernel` type. It contains all information necessary to construct the 3 main Smolyak outputs:

- `SmolyakGrid` : Smolyak Grid
- `SmolyakBasis` : Smolyak Basis
- `SmolyakPoly` : Smolyak Polynomial

## Example

Below is an example of how to use module to:

1. Create a Smolyak Kernel
2. Define Smolyak Basis functions of an Interpolating Smolyak Polynomial on the Smolyak Grid
3. Create and evaluate the Smolyak Polynomial

```
using Smolyak

# Accuracy levels and bounds on x variables
mu = [2,2]
xbnds = [[-2. , 3.] for d in 1:length(mu)]

# Create a Smolyak Kernel 
sk = SmolyakKernel(mu, xbnds);
```

Using the SmolyakKernel to create Smolyak Grids and Polynomials.

# Creating a Smolyak grid

```
# Create Smolyak grid corresponding to Smolyak Kernel above
sg = SmolyakGrid(sk);

# View the vector of grid points
sg.grid

# Convert grid from vector of vectors to a matrix (i.e. for plot)
grid = VVtoMatrix(sg.grid)
```

See [SmolyakGridExample.jl](./Examples/SmolyakGridExample.) for plots of 2-dimensional Smolyak Grids. 

![](./Examples/IsotrophicSmolyakGridExample.png)
![](./Examples/AnisotrophicSmolyakGridExample.png)

# Creating a Smolyak Polynomial

Smolyak polynomials require a SmolyakKernel and choice of basis function to use when building up tensor products. 

- SmolyakKernel
- SmolyakBasis defined using sk at state vector x 
- Coefficients

The Smolyak basis can be automatically called in the SmolyakPoly() call. However, the user must specify type of basis function to use: ∈{:ordinary, :chebyshev, :spread} in the SmolyakPoly() function.
 
Initiate Basis function - with memory allcoated for up to 2 derivatives (optional)

```
# Choose type of polynomial basis function
basis_fun_type = :chebyshev;

# Option 1: Setup Smolyak polynomial using Smolyak basis 
sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=2);
sp = SmolyakPoly(sb; NumDeriv=2);

# Option 2: Setup Smolyak polynomial directly
sp = SmolyakPoly(basis_fun_type, sk; NumDeriv=2);

# Initialise coefficients: update Smolyak polynomial at stored x
θ = rand(length(sp.sb.BF));
makeSmolyakPoly!(θ, sp; NumDeriv=2);

# Print updated fields
sp.value
sp.gradient
sp.hessian

# Add a new x: update Smolyak polynomial at stored θ
x = [-1.0, 2.3];
makeSmolyakPoly!(sp, x; NumDeriv=2);

# Print updated fields
sp.value
sp.gradient
sp.hessian

# Add a new θ  & x: update Smolyak polynomial
θ = rand(length(sp.sb.BF));
x = [-1.0, 2.3];
makeSmolyakPoly!(θ, sp, x; NumDeriv=2);

# Print updated fields
sp.value
sp.gradient
sp.hessian

# Can also calculate individual term wihout full update 
sp.value[1] == getValue(sp)
sp.gradient[1] == get_dWdx(sp, 1) 
sp.gradient[2] == get_dWdx(sp, 2) 
sp.hessian[1,1] == get_d2Wdx2(sp, 1, 1) 
sp.hessian[1,2] == get_d2Wdx2(sp, 1, 2) 
sp.hessian[2,2] == get_d2Wdx2(sp, 2, 2) 
```


Further examples can be found in [Interpolation Example](./Examples/Interpolation_Example.jl).