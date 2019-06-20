#dev LInearVERSION < v"0.7.0-beta2.199" && __precompile__()

"""
## Description

Module contains code to define Smolyak Grids and construct corresponding Interpolating Smolyak 
Polynomials. Both Anisotrophic and Isotrophic Grids/Polynomials are supported and are constructed 
efficiently following the methodology outlined in Judd, Maliar, Maliar, Valero (2014). 

It supports construction of Smolyak Polynomials with Ordinary, Chebyshev and Spread Polynomial basis functions.

The code is designed for Julia version: 1.0.

"""
module Smolyak

using Base.Cartesian
using LinearAlgebra
using Combinatorics

using SparseArrays: findnz

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
export 	VVtoMatrix, VV, VM,
		SmolyakKernel, BasisIdx!,
		HDSmolIdx, NumGridPts,
		SmolyakKernel,
		SmolyakGrid, 
		grid!,
		grid, xgrid, 
		lineartransform, x2z, z2x, dxdz, dzdx,
		OrdinaryPoly, 
		ChebyshevPoly, 
		SpreadPoly,
		SmolyakBasis, 
		state!, BasisFunctions!, jacobian!, hessian!, SmolyakBasis!,
		state, BasisFunctions, jacobian, hessian, SBoutput,
		SmolyakPoly, 
		coef!, value!, dWdx!, gradient!, d2Wdx2!, SmolyakPoly!,		
		coef, value, dWdx, gradient, d2Wdx2, SPoutput

end
