#= ***************** =#
#= Smolyak Grid Type =#
#= ***************** =#

struct SmolyakGrid{T<:Real} 
	sk 		:: SmolyakKernel	# SmolyakKernel
	N 		:: Int64 			# Number of Grid Points
	grid  	:: VV{T}			# Smolyak Grid on z <- x2z!() (i.e. if doing collocation, usine Chebyshev and [z_lb,z_ub]∈[-1,1])
end

# Call with Smolyak Kernal setup
function SmolyakGrid(sk::SmolyakKernel) where T<:SmolyakKernelType
	
	NGP = makeNumGridPts(sk.GridIdx,(sk.mu...,)) 
	grid = Vector{Float64}[Array{Float64}(undef,sk.D) for r in 1:NGP]
	makeGrid!(grid,sk.GridIdx,tuple(sk.mu...,))			# Make Grid on [-1,1]
	
	return SmolyakGrid(sk, NGP, grid)
end

# Set up Smolyak Kernel internally
function SmolyakGrid(D::Int64, mu_level::T; HD::Bool=false) where T<:Real
	
	sk = SmolyakKernel(D, mu_level; HD=HD)
	NGP = makeNumGridPts(sk.GridIdx,(sk.mu...,)) 
	grid = Vector{Float64}[Array{Float64}(undef,sk.D) for r in 1:NGP]
	makeGrid!(grid,sk.GridIdx,tuple(sk.mu...,))			# Make Grid on [-1,1]

	return SmolyakGrid(sk, NGP, grid)
end

function SmolyakGrid(mu::Vector{T}; HD::Bool=false) where T<:Real

	sk = SmolyakKernel(mu; HD=HD)
	NGP = makeNumGridPts(sk.GridIdx,(sk.mu...,)) 
	grid = Vector{Float64}[Array{Float64}(undef,sk.D) for r in 1:NGP]
	makeGrid!(grid,sk.GridIdx,tuple(sk.mu...,))			# Make Grid on [-1,1]
	
	return SmolyakGrid(sk, NGP, grid)
end

# Set up Smolyak Kernel internally
function SmolyakGrid(D::Int64, mu_level::T, xbnds::Matrix{S}; HD::Bool=false) where T<:Real where S<:Real
	
	sk = SmolyakKernel(D, mu_level, xbnds; HD=HD)
	NGP = makeNumGridPts(sk.GridIdx,(sk.mu...,)) 
	grid = Vector{Float64}[Array{Float64}(undef,sk.D) for r in 1:NGP]
	makeGrid!(grid,sk.GridIdx,tuple(sk.mu...,))			# Make Grid on [-1,1]

	return SmolyakGrid(sk, NGP, grid)
end

# Set up Smolyak Kernel internally
function SmolyakGrid(mu::Vector{T}, xbnds::Matrix{S}; HD::Bool=false) where T<:Real where S<:Real
	
	sk = SmolyakKernel(mu, xbnds; HD=HD)
	NGP = makeNumGridPts(sk.GridIdx,(sk.mu...,)) 
	grid = Vector{Float64}[Array{Float64}(undef,sk.D) for r in 1:NGP]
	makeGrid!(grid,sk.GridIdx,tuple(sk.mu...,))			# Make Grid on [-1,1]

	return SmolyakGrid(sk, NGP, grid)
end

#= ******************************************************** =#
#= 		Functions to make the Smolyak Grid 				    =#
#= ******************************************************** =#

@generated function makeGrid!(H::VV{Float64},inds::VV{Int64},mu::NTuple{N,T}) where N where T<:Real
	quote 
		zH = 0
		for i = 1:length(inds)
			for j in Iterators.product((@ntuple $N k->grid_A_i(inds[i][k]))...)
				zH += 1
				for w in 1:$N
					H[zH][w] = j[w]
				end
			end
		end
		return H
	end
end
