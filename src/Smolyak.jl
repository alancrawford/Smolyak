VERSION < v"0.7.0-beta2.199" && __precompile__()

"""
## Description

Module contains code to define Smolyak Grids and construct corresponding Interpolating Smolyak 
Polynomials. Both Anisotrophic and Isotrophic Grids/Polynomials are supported and are constructed 
efficiently following the methodology outlined in Judd, Maliar, Maliar, Valero (2014). 

It supports construction of Smolyak Polynomials with Ordinary, Chebyshev and Spread Polynomial basis functions.

The code is designed for Julia version: 0.7.

"""
module Smolyak

using Base.Cartesian
using LinearAlgebra
using Combinatorics

# Type aliases 
VV{T} = Vector{Vector{T}}
VM{T} = Vector{Matrix{T}}

# Load Files
include("utils.jl")
include("smolkernel.jl")
include("smolyakHD.jl")
include("smolgrid.jl")
include("BasisFunctions.jl")
include("OrdinaryPolynomials.jl")
include("ChebyshevPolynomials.jl")
include("SpreadPolynomials.jl")
include("smolbasis.jl")
include("smolpoly.jl")

# Make functions defined in the module available =#
export 	SmolyakKernel, makeBasisIdx!,
		makeHDSmolIdx, makeNumGridPts,
		SmolyakGrid, makeGrid!, 
		lineartransform, x2z, z2x, dxdz, dzdx,
		OrdinaryPoly, ChebyShevPoly, SpreadPoly, makeBF!, 
		SmolyakBasis, new_x!, makeJacobian!, makeHessian!, makeSmolyakBasis!,
		SmolyakPoly, update_coef!, makeSmolyakPoly!, makeValue!, getValue, 
		get_dWdx, make_dWdx!, makeGradient!,
		get_d2Wdx2, make_d2Wdx2!, makeHessian!, makeSmolyakPoly!,
		VVtoMatrix, VV, VM

end
