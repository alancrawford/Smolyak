module Smolyak

using Base.Cartesian, Base.LinAlg, Iterators
import Base.show 

typealias ScalarOrVec{T} Union(T,Vector{T})
typealias VecOrArray{T} Union(Vector{T},Matrix{T})

# Load Files
include("smolgrid.jl")
include("smolbasis.jl")
include("smolpoly.jl")

# Make functions defined in the module available =#
export 	SmolyakGrid, xGrid!,x2z, z2x, x2z!,	z2x!,
		SmolyakBasis, makeBF!, new_x!, 
		SmolyakPoly, dfdx!, d2fdx2!, new_coef!, wgt_new_coef!

end
