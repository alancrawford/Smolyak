#= 
	-------------------------------------
	Smolyak Grid for Julia version 0.4  
	-------------------------------------

This file contains code to define Smolyak Grid type. Both Anisotrophic and 
Isotrophic Grids are supported and they are constructed efficiently 
following the methodology outlined in JMMV (2014). The code is designed on 
latest Julia version: 0.4.

Key Refs: JMMV (2014), Burkhardt (2012)

=#

#= ***************** =#
#= Smolyak Grid Type =#
#= ***************** =#

type SmolyakGrid
	D 			::	Int64				# Dimensions
	mu 			::	ScalarOrVec{Int64}	# Index of mu
	NumGrdPts 	::	Int64				# Number of Grid Points
	lb 			::	Vector{Float64}		# Lower Bounds of dimensions
	ub 			::	Vector{Float64}		# Upper Bounds of dimensions
	zGrid 		::	AA{Float64}			# Smolyak Grid on z = [-1,1]
	xGrid 		::	AA{Float64}			# Smolyak Grid on original domain x in [lb,ub]
	Binds    	::  AA{Int64}			# Input to construct Basis Funs for set of grid points -> Will depend on mu.

	"Constructor function for Smolyak Grid"
	function SmolyakGrid(D::Int64, mu::ScalarOrVec{Int64},lb::Vector{Float64}=-1*ones(Float64,D), ub::Vector{Float64}=ones(Float64,D))
		
		# Setup
		NumGrdPts, Ginds = SmolIdx(tuple(mu...))
		z = Vector{Float64}[Array{Float64}(D) for r in 1:NumGrdPts]
		x = Vector{Float64}[Array{Float64}(D) for r in 1:NumGrdPts]
		Binds = Vector{Int64}[Array{Int64}(D) for r in 1:NumGrdPts]

		# Make Grid and Indices
		makeGrid!(z,Ginds,tuple(mu...))			# Make Grid on [-1,1]
		z2x!(z,x,lb,ub) 						# Grid on X
		makeBasisIdx!(Binds,Ginds,tuple(mu...)) # Basis Function Indices

		new(D, mu, NumGrdPts, lb, ub, z, x, Binds)
	end
end

function show(io::IO, sg::SmolyakGrid)
	if !=(minimum(sg.mu),maximum(sg.mu))
			mu_print = strip(string(sg.mu))
			msg = "Anisotrophic Smolyak Grid:\n"
			msg *= "\tD: $(sg.D)\n\tmu: $(mu_print)\n\tNum Grid Points: $(sg.NumGrdPts)"
	else
		msg = "Isotrophic Smolyak Grid:\n"
		msg *= "\tD: $(sg.D)\n\tmu: $(sg.mu[1])\n\tNum Grid Points: $(sg.NumGrdPts)"
	end
	print(io, msg)
end

#= ***************************************** =#
#= Sub Funs called when creating SmolyakGrid =#
#= ***************************************** =#

# This fn returns theta = cos(x), where x are chebyshev nodes
function chebtheta(n::Int64)	
	is(n,1) ? 0.5pi : nothing 
	return [0:n-1;]pi./(n-1)
end	

function chebnodes(n::Int64)	
	return round(cos(chebtheta(n)),14)
end

# This is the map from index to number of grid points in each dimension
function m_i(i::Int64)
	is(i,1) ? 1 : 2.^(i-1)+1
end

# This is number of new points added by nested rule as function of index i ( = l+1 in Burkhardt)
function newpts(i::Int64)
	if <=(i,2)
		return i 
	else
		return 2.^(i-2)
	end
end

# Disjoint Sets of Grid Points whose product are combined 
function grid_A_i(i::Int64)
	if is(i,1)
		return 0.0
	elseif is(i,2)
		return chebnodes(m_i(i))[1:2:m_i(i)]
	else
		return chebnodes(m_i(i))[2:2:m_i(i)] 
	end 	
end

#= ******************************* =#
#= Funs used to create SmolyakGrid =#
#= ******************************* =#

# Calculate Number of Grid Points & Construct Indices of Smolyak Grid
@generated function SmolIdx{N}(mu::NTuple{N,Int64} )
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
	
# Constructs the Smolyak Grid
@generated function makeGrid!{N}(H::AA{Float64},inds::AA{Int64},mu::NTuple{N,Int64})
	quote 
		zH = 0
		for i = 1:length(inds)
			for j in product((@ntuple $N k->grid_A_i(inds[i][k]))...)
				zH += 1
				for w in 1:$N
					H[zH][w] = j[w]
				end
			end
		end
		return H
	end
end

#= ******************************************* =#
#= Create Basis Funs Indices SmolyakGrid(d,mu) =#
#= ******************************************* =#

# Disjoint sets that define indexes for creation of Basis Indices
function A_pidx(ibar::Int64)
	A = []
	lb 	= [1;2;2.^(collect(3:ibar)-1) - 2.^(collect(3:ibar)-2)+2]
	ub 	= [1;2.^(collect(2:ibar)-1)+1]
	for j = 1:ibar
		push!(A,collect(lb[j]:ub[j]))
	end
	return A
end

# Make Basis Function Indexes
@generated function makeBasisIdx!{N}(Binds::AA{Int64},GridIdx::AA{Int64},mu::NTuple{N,Int64})
	quote
		max_mu=0
		for i in 1:$N
			max_mu = max(max_mu,mu[i])
		end
		A = A_pidx(max_mu+1)
		zT = 0
		for i = 1:length(GridIdx)
			for j in product((@ntuple $N k->A[GridIdx[i][k]])...)
				zT += 1
				for w in 1:$N
					Binds[zT][w] = j[w]-1
				end
			end
		end
		return Binds
	end
end

#= ******************************************************** =#
#= Functions switching between z in [-1,1] and x in [lb,ub] =#
#= ******************************************************** =#

# In place coordinate transform for vector: z→x
function z2x!(sg::SmolyakGrid)
	for n in eachindex(sg.zGrid), d in 1:sg.D
		sg.xGrid[n][d] = 0.5*( sg.ub[d] + sg.lb[d] + sg.zGrid[n][d]*(sg.ub[d] - sg.lb[d]) )
	end
end

# In place coordinate transform for vector: z→x
function z2x!(z::Vector{Float64},x::Vector{Float64},lb::Vector{Float64},ub::Vector{Float64})
	for d in 1:length(z)
		x[d] = 0.5*( ub[d] + lb[d] + z[d]*(ub[d] - lb[d]) )
	end
end

# In place coordinate transform for Array: z→x
function z2x!(z::AA{Float64},x::AA{Float64},lb::Vector{Float64},ub::Vector{Float64})
	for n in eachindex(z), d in eachindex(z[1])
		x[n][d] = 0.5*( ub[d] + lb[d] + z[n][d]*(ub[d] - lb[d]) )
	end
end

# In place coordinate transform for vector: x→z
function x2z!(sg::SmolyakGrid)
	for n in eachindex(sg.xGrid), d in 1:sg.D
		sg.zGrid[n][d] = (2sg.xGrid[n][d] - sg.ub[d] - sg.lb[d])/(sg.ub[d] - sg.lb[d])
	end
end

# In place coordinate transform for vector: x→z
function x2z!(x::Vector{Float64},z::Vector{Float64},lb::Vector{Float64},ub::Vector{Float64})
	for d in 1:length(x)
		z[d] = (2x[d] - ub[d] - lb[d])/(ub[d] - lb[d]) 
	end
end

# In place coordinate transform for vector: x→z
function x2z!(x::AA{Float64},z::AA{Float64},lb::Vector{Float64},ub::Vector{Float64})
	for n in eachindex(z), d in eachindex(x[1])
		z[n][d] = (2x[n][d] - ub[d] - lb[d])/(ub[d] - lb[d]) 
	end
end