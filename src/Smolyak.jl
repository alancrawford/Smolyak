module Smolyak

using Base.Cartesian, Base.LinAlg, Iterators
import Base.show 

# Load Files
include("smolgrid.jl")
include("smolbasis.jl")
include("smolpoly.jl")

# Make functions defined in the module available =#
export 	SmolyakGrid,
		xGrid!,
		SmolyakBasis,
		x2z,
		z2x,
		x2z!,
		z2x!,		 
		SmolyakPoly,
		new_sp_coef!,
		wgt_sp_coef!,
		new_sp_f!
end
