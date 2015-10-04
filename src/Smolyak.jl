"""
## Description

Module contains code to define Smolyak Grids and construct corresponding Interpolating Smolyak 
Polynomials. Both Anisotrophic and Isotrophic Grids/Polynomials are supported and are constructed 
efficiently following the methodology outlined in Judd, Maliar, Maliar, Valero (2014). 

The code is designed for Julia version: 0.4.

#### Types 

3 main components with corresponding types:

- `SmolyakGrid` : Smolyak Grid
- `SmolyakBasis` : Smolyak Basis
- `SmolyakPoly` : Smolyak Polynomial

See help functions in REPL for details.

**Notes**: Utilises following typealiases:

```julia
typealias AA{T} Array{Array{T,1},1}
typealias AM{T} Array{Matrix{T},1}                                                                                 
typealias AAA{T} Array{Array{Array{T,1},1},1}
typealias AAAA{T} Array{Array{Array{Array{T,1},1},1},1}
typealias ScalarOrVec{T} Union{T,Vector{T}}
typealias VecOrAA{T} Union{Vector{T},AA{T}}
```
"""
module Smolyak

using Base.Cartesian, Base.LinAlg, Iterators
import Base.show 

typealias AA{T} Array{Array{T,1},1}
typealias AM{T} Array{Matrix{T},1}                                                                                 
typealias AAA{T} Array{Array{Array{T,1},1},1}
typealias AAAA{T} Array{Array{Array{Array{T,1},1},1},1}
typealias ScalarOrVec{T} Union{T,Vector{T}}
typealias VecOrAA{T} Union{Vector{T},AA{T}}

# Load Files
include("smolgrid.jl")
include("smolbasis.jl")
include("smolpoly.jl")

# Make functions defined in the module available =#
export 	SmolyakGrid, makeGrid!, x2z, z2x, x2z!, z2x!,
		SmolyakBasis, makeBasis!, new_x!, 
		SmolyakPoly, getValue!, getGrad!, getHess!, get_pinvBFt!, getCoef!

end
