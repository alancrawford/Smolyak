#= 
	-------------------------------------
	Smolyak Grid for Julia version 0.3.7  
	-------------------------------------

This file contains code to define Smolyak Grid type. Both Anisotrophic and 
Isotrophic Grids are supported and they are constructed efficiently 
following the methodology outlined in JMMV (2014). The code is designed on 
latest stable Julia version: 0.3.7.

Key Refs: JMMV (2014), Burkhardt (2012), github: ECONFORGE/Smolyak

=#

IntOrVec = Union(Int64,Vector{Int64},Float64,Vector{Float64})
VecOrArray = Union(Vector{Float64},Array{Float64,2}) 

#= ***************************************** =#
#= Sub Funs called when creating SmolyakGrid =#
#= ***************************************** =#

# This fn returns theta = cos(x), where x are chebyshev nodes
function chebtheta(n::Int64)	
	==(n,1) ? 0.5pi : nothing 
	return [0:n-1]pi./(n-1)
end	

function chebnodes(n::Int64)	
	return round(cos(chebtheta(n)),14)
end

# This is the map from index to number of grid points in each dimension
function m_i(i::Int64)
	==(i,1) ? 1 : 2.^(i-1)+1
end

# This is number of new points added by nested rule as function of index i ( = l+1 in Burkhardt)
function newpts(i::Int64)
	if i==1
		return 1 
	elseif i==2 
		return 2
	else
		return 2.^(i-2)
	end
end

# Disjoint Sets of Grid Points whose product are combined 
function grid_A_i(i::Int64)
	if ==(i,1)
		return 0.0
	elseif ==(i,2)
		return chebnodes(m_i(i))[1:2:m_i(i)]
	else
		return chebnodes(m_i(i))[2:2:m_i(i)] 
	end 	
end

#= ******************************* =#
#= Funs used to create SmolyakGrid =#
#= ******************************* =#

function µfun(d::Int64,µ::IntOrVec)
	if isa(µ,Int64)
		maxµ = µ
		µbar = tuple([µ*ones(Int64,d)]...)
	else
		maxµ = 0
		for m = 1:length(µ)
			maxµ = max(maxµ,µ[m])
		end
		µbar = tuple(µ...) 
	end
	@assert(d == length(µbar),"\n\tError: Dimension mismatch: \n
								\tNumber of dimensions!= length of precision index (µ)\n")
	return µbar , maxµ
end

# Calculate Number of Grid Points & Construct Indices of Smolyak Grid
@ngenerate N (Int64,Array{Array{Int64,1},1}) function SmolIdx(maxµ::Int64,µ::NTuple{N,Int64}...)
	µbar = Int64[]
	ibar = Int64[]
	n = Int64[]
	NumGridPts = 0
	sum_i = 0 
	iList=Array{Int64,1}[]	
	@nexprs N j-> push!(µbar,µ_j)
	@nloops N i j->1:µbar[j]+1 begin
		@nexprs N j -> push!(n,newpts(i_j))
		@nexprs N j -> sum_i += i_j 
		if sum_i > N + maxµ
			ibar = Int64[]
			n = Int64[]
			sum_i=0
			break
		end
		@nexprs N j-> push!(ibar, i_j)
		push!(iList,ibar)
		NumGridPts += prod(n)
		ibar = Int64[]
		n = Int64[]
		sum_i=0
	end
	return NumGridPts, iList
end

# Constructs the Smolyak Grid
@ngenerate N (Array{Float64,2}) function Make_Grid(NGP::Int64,inds::Array{Array{Int64,1},1},µ::NTuple{N,Int64}...)
	H = Array(Float64,N, NGP)
	zH = 0
	for i = 1:length(inds)
		for j in product((@ntuple N k->grid_A_i(inds[i][k]))...)
			zH += 1
			H[:,zH] = [j...]
		end
	end
	return H
end

#= ******************************************* =#
#= Create Basis Funs Indices SmolyakGrid(d,µ) =#
#= ******************************************* =#

# Disjoint sets that define indexes for creation of Basis Indices
function A_pidx(ibar::Int64)
	A = {}
	#= A = [] in version 0.4 =#
	lb 	= 	[	1,
				2,
				2.^([3:ibar]-1) - 2.^([3:ibar]-2)+2]
	ub 	= 	[	1,
				2.^([2:ibar]-1)+1]
	for j = 1:ibar
		push!(A,[lb[j]:ub[j]])
	end
	return A
end

@ngenerate N Array{Int64,2} function BasisIdx(NGP::Int64,GridIdx::Array{Array{Int64,1},1},maxµ::Int64,µ::NTuple{N,Int64}...)
	A = A_pidx(maxµ+1)
	T = Array(Int64,N,NGP)
	zT = 0
	for i = 1:length(GridIdx)
		for j in product((@ntuple N k->A[GridIdx[i][k]])...)
			zT += 1
			T[:,zT] = [j...]-1
		end
	end
	return T
end

#= ***************** =#
#= Smolyak Grid Type =#
#= ***************** =#

type SmolyakGrid
	D 			::	Int64				# Dimensions
	µ 			::	IntOrVec			# Index of µ
	NumGrdPts 	::	Int64				# Number of Grid Points
	lb 			::	Vector{Float64}		# Lower Bounds of dimensions
	ub 			::	Vector{Float64}		# Upper Bounds of dimensions
	zGrid 		::	Array{Float64,2}	# Smolyak Grid on z = [-1,1]
	xGrid 		::	Array{Float64,2}	# Smolyak Grid on original domain x in [lb,ub]
	Binds    	::  Array{Int64,2}		# Input to construct Basis Funs for set of grid points -> Will depend on µ.

	function SmolyakGrid(	D::Int64, µ::IntOrVec,
							lb::Vector{Float64}=-1*ones(Float64,D), ub::Vector{Float64}=ones(Float64,D))
		
		#= z in [-1,1]^D and x[d] in [lb[d],ub[d]] for d = 1,...,D =#
		
		@assert(is(size(lb,1),D) || is(size(ub,1),D),
				"\n\tError: [lb,ub] does not equal dimension of grid")
		#= @assert(.<=(lb,ub),
				"\n\tError: lb >= ub") =#

		µbar, maxµ = µfun(D,µ) 
		NumGrdPts, Ginds = SmolIdx(maxµ,µbar...)
		zGrid = Make_Grid(NumGrdPts,Ginds,µbar...)	
		xGrid  = z2x(zGrid,lb,ub)
		Binds = BasisIdx(NumGrdPts,Ginds,maxµ,µbar...)

		new(D, µ, NumGrdPts, lb, ub, zGrid, xGrid, Binds)
	end
end

function show(io::IO, sg::SmolyakGrid)
	if !=(minimum(sg.µ),maximum(sg.µ))
			µ_print = strip(string(sg.µ))
			msg = "Anisotrophic Smolyak Grid:\n"
			msg *= "\tD: $(sg.D)\n\tµ: $(µ_print)\n\tNum Grid Points: $(sg.NumGrdPts)"
	else
		msg = "Isotrophic Smolyak Grid:\n"
		msg *= "\tD: $(sg.D)\n\tµ: $(sg.µ[1])\n\tNum Grid Points: $(sg.NumGrdPts)"
	end
	print(io, msg)
end

#= Functions switching between z in [-1,1] and x in [lb,ub] 
	- See ECONFORGE/Smolyak for unadjusted code =#

function z2x(zpts::VecOrArray,lb::Vector{Float64},ub::Vector{Float64})
	centers = lb + (ub - lb)./2
    radii = (ub - lb)./2
    return centers .+ zpts.*radii
end

function x2z(xpts::VecOrArray,lb::Vector{Float64},ub::Vector{Float64})
 	centers = lb + (ub - lb)./2
    radii = (ub - lb)./2
    return (xpts .- centers)./radii
end

function z2x!(sg::SmolyakGrid)
    sg.xGrid = (sg.lb + (sg.ub - sg.lb)./2) .+ sg.zGrid.*((sg.ub - sg.lb)./2)
end

function x2z!(sg::SmolyakGrid)
    sg.zGrid = (sg.xGrid .- (sg.lb + (sg.ub - sg.lb)./2))./((sg.ub - sg.lb)./2)
end

function xGrid!(sg::SmolyakGrid, NewX::Array{Float64,2}, K::Int64=1, IsIbar::Bool=true)
	is(IsIbar,true) ?
		sg.xGrid[1:size(NewX,1),:] = NewX :
		sg.xGrid[K+1:K+size(NewX,1),:] = NewX
end


