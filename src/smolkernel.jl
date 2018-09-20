

abstract type SmolyakKernelType end

struct SmolyakKernel{T} <: SmolyakKernelType
	D 			:: 	Int64 				# Dimensions
	mu 			::	Vector{Int64}		# Index of mu
	N 			:: 	Int64 				# Number of Grid Points
	xbnds		::	VV{T}			# Bounds of dimensions of x
	zbnds		::	VV{T}			# Bounds of dimensions of z
	GridIdx 	::  VV{Int64}             # Input for grid constriction
	BasisIdx	::  VV{Int64}				# Input to construct Basis Funs for set of grid points -> Will depend on mu.
end

# --------------------- Constructor functions -------------------------- #

function SmolyakKernel(D::Int64, mu_level::T; HD::Bool=false) where T<:Real
	
	muvec =  mu_level*ones(typeof(mu_level),D)
	if HD 
		GridIdx = HDSmolIdx(muvec)
		NGP = NumGridPts(GridIdx,(muvec...,))
	else
		NGP, GridIdx = SmolIdx(tuple(muvec...,))
	end

	xbnds = [[-1.0,1.0] for d in 1:D]
	zbnds = [[-1.0,1.0] for d in 1:D]
	BasisIdx = Vector{Int64}[Array{Int64}(undef,D) for r in 1:NGP]	
	BasisIdx!(BasisIdx,GridIdx,tuple(muvec...,)) # Basis Function Indices

	return SmolyakKernel(D, muvec, NGP, xbnds, zbnds, GridIdx , BasisIdx )
end

function SmolyakKernel(mu::Vector{T}; HD::Bool=false) where T<:Real
	
	# Setup
	D = length(mu)
	if HD 
		GridIdx = HDSmolIdx(mu)
		NGP = NumGridPts(GridIdx,(mu...,))
	else
		NGP, GridIdx = SmolIdx(tuple(mu...,))
	end
	
	# By default set both x and z bounds to [-1,1]
	xbnds = [[-1.0,1.0] for d in 1:D]
	zbnds = [[-1.0,1.0] for d in 1:D]
	BasisIdx = Vector{Int64}[Array{Int64}(undef,D) for r in 1:NGP]

	# Make Indices
	BasisIdx!(BasisIdx,GridIdx,tuple(mu...,)) # Basis Function Indices

	return SmolyakKernel(D, mu, NGP, xbnds, zbnds, GridIdx , BasisIdx )
end

function SmolyakKernel(D::Int64, mu_level::T, xbnds::VV{S}; HD::Bool=false) where T<:Real where S<:Real

	muvec = mu_level*ones(typeof(mu_level),D)
	
	if HD 
		GridIdx = HDSmolIdx(muvec)
		NGP = NumGridPts(GridIdx,(muvec...,))
	else
		NGP, GridIdx = SmolIdx(tuple(muvec...,))
	end
	
	# By default set both x and z bounds to [-1,1]
	zbnds = [[-one(S),one(S)] for d in 1:D]
	BasisIdx = Vector{Int64}[Array{Int64}(undef,D) for r in 1:NGP]

	# Make Indices
	BasisIdx!(BasisIdx,GridIdx,tuple(muvec...,)) # Basis Function Indices

	return SmolyakKernel(D, muvec, NGP, xbnds, zbnds, GridIdx , BasisIdx )
end

function SmolyakKernel(mu::Vector{T}, xbnds::VV{S}; HD::Bool=false) where T<:Real where S<:Real
	
	# Setup
	D = length(mu)
	if HD 
		GridIdx = HDSmolIdx(mu)
		NGP = NumGridPts(GridIdx,(mu...,))
	else
		NGP, GridIdx = SmolIdx(tuple(mu...,))
	end
	
	# By default set both x and z bounds to [-1,1]
	zbnds = [[-one(S),one(S)] for d in 1:D]
	BasisIdx = Vector{Int64}[Array{Int64}(undef,D) for r in 1:NGP]

	# Make Indices
	BasisIdx!(BasisIdx,GridIdx,tuple(mu...,)) # Basis Function Indices

	return SmolyakKernel(D, mu, NGP, xbnds, zbnds, GridIdx , BasisIdx )
end

function SmolyakKernel(D::Int64, mu_level::T,  xbnds::VV{S}, zbnds::VV{S}; HD::Bool=false) where T<:Real where S<:Real

	muvec = mu_level*ones(T,D)
	
	if HD 
		GridIdx = HDSmolIdx(muvec)
		NGP = NumGridPts(GridIdx,(muvec...,))
	else
		NGP, GridIdx = SmolIdx(tuple(muvec...,))
	end
	
	# By default set both x and z bounds to [-1,1]
	BasisIdx = Vector{Int64}[Array{Int64}(undef,D) for r in 1:NGP]

	# Make Indices
	BasisIdx!(BasisIdx,GridIdx,tuple(muvec...,)) # Basis Function Indices

	return SmolyakKernel(D, muvec, NGP, xbnds, zbnds, GridIdx , BasisIdx )
end

function SmolyakKernel(mu::Vector{T}, xbnds::VV{S}, zbnds::VV{S}; HD::Bool=false) where T<:Real where S<:Real
	
	# Setup
	D = length(mu)
	if HD 
		GridIdx = HDSmolIdx(mu)
		NGP = NumGridPts(GridIdx,(mu...,))
	else
		NGP, GridIdx = SmolIdx(tuple(mu...,))
	end
	
	# By default set both x and z bounds to [-1,1]
	BasisIdx = Vector{Int64}[Array{Int64}(undef,D) for r in 1:NGP]

	# Make Indices
	BasisIdx!(BasisIdx,GridIdx,tuple(mu...,)) # Basis Function Indices

	return SmolyakKernel(D, mu, NGP, xbnds, zbnds, GridIdx , BasisIdx )
end

#= ***************************************** =#
#= Sub Funs called when creating SmolyakGrid =#
#= ***************************************** =#

# This fn returns theta = cos.(x), where x are chebyshev nodes
function chebtheta(n::Int64)	
	===(n,1) ? 0.5pi : nothing 
	return [0:n-1;]pi./(n-1)
end	

function chebnodes(n::Int64)	
	return round.(cos.(chebtheta(n)); digits=14)
end

# This is the map from index to number of grid points in each dimension
function m_i(i::Int64)
	===(i,1) ? 1 : 2 .^(i-1)+1
end

# This is number of new points added by nested rule as function of index i ( = l+1 in Burkhardt)
function newpts(i::Int64)
	if <=(i,2)
		return i 
	else
		return 2 .^(i-2)
	end
end

# Disjoint Sets of Grid Points whose product are combined 
function grid_A_i(i::Int64)
	if ===(i,1)
		return 0.0
	elseif ===(i,2)
		return chebnodes(m_i(i))[1:2:m_i(i)]
	else
		return chebnodes(m_i(i))[2:2:m_i(i)] 
	end 	
end

#= ********************************* =#
#= Funs used to create SmolyakKernel =#
#= ********************************* =#

# Calculate Number of Grid Points & Construct Indices of Smolyak Grid
@generated function SmolIdx(mu::NTuple{N,T} ) where N where T<:Real
	quote
		max_mu=0
		for i in 1:$N
			max_mu = max(max_mu,mu[i])
		end
		iList= Array{Int64,1}[]
		ibar = Int64[]
		n = Int64[]
		NumGridPts = 0
		sum_i = 0 
		@nloops $N i j->1:mu[j]+1 begin
			@nexprs $N j -> push!(n,newpts(i_j))
			@nexprs $N j -> sum_i += i_j 
			if sum_i > $N + max_mu
				ibar = Int64[]
				n = Int64[]
				sum_i=0
				break
			end
			@nexprs $N j-> push!(ibar, i_j)
			push!(iList,ibar)
			NumGridPts += prod(n)
			ibar = Int64[]
			n = Int64[]
			sum_i=0
		end
		return NumGridPts, iList
	end
end

#= ******************************************* =#
#= Create Basis Funs Indices SmolyakGrid	 =#
#= ******************************************* =#

# Disjoint sets that define indexes for creation of Basis Indices
function A_pidx(ibar::Int64)
	A = []
	lb 	= [1;2]
	append!(lb, 2 .^(collect(3:ibar) .-1) .- 2 .^(collect(3:ibar) .-2) .+2)
	ub 	= [1]
	append!(ub, 2 .^(collect(2:ibar) .-1) .+1)
	for j = 1:ibar
		push!(A,collect(lb[j]:ub[j]))
	end
	return A
end

# --------------------- Updating functions -------------------------- #

# Make Basis Function Indexes
@generated function BasisIdx!(Binds::VV{Int64},GridIdx::VV{Int64},mu::NTuple{N,T}) where N where T<:Real
	quote
		max_mu=0
		for i in 1:$N
			max_mu = max(max_mu,mu[i])
		end
		A = A_pidx(max_mu+1)
		zT = 0
		for i = 1:length(GridIdx)
			for j in Iterators.product((@ntuple $N k->A[GridIdx[i][k]])...)
				zT += 1
				for w in 1:$N
					Binds[zT][w] = j[w]-1
				end
			end
		end
		return Binds
	end
end

