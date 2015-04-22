module Smolyak

using Base.Cartesian, Base.LinAlg, Iterators
import Base.show 

# Load Files
include("smolgrid.jl")
include("smolbasis.jl")
include("smolpoly.jl")

# Make functions defined in the module available =#
export 	SmolyakGrid,
		SmolyakBasis,
		x2z,
		z2x, 
		SmolyakPoly,
		update_spoly_coef!,
		update_spoly_f!
end
