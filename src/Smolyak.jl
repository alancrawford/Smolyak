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
export 	SmolyakGrid, x2z, z2x, x2z!, z2x!,
		SmolyakBasis, makeBF!, new_x!, 
		SmolyakPoly, getValue!, getGrad!, getHess!, get_pinvBFt!, getCoef!

end
