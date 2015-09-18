#= 
	-----------------------------------------------------------------
	Smolyak Basis Functions on a Smolyak Grid for Julia version 0.3.7  
	-----------------------------------------------------------------

This file contains code to define Smolyak Constructtion of Chebyshev
basis functions. It is compatible with both Anisotrophic and Isotrophic 
Grids and are constructed efficiently following the methodology outlined
 in JMMV (2014). The code is designed on latest stable Julia version: v0.4.

Key Refs: JMMV (2014), Burkhardt (2012), github: ECONFORGE/Smolyak

=#


#= ******************* =#
#= Smolyak Basis type  =#
#= ******************* =#

type SmolyakBasis
	D 			:: Int64					# Dimensions
	mu 			:: ScalarOrVec{Int64}		# Index of mu
	lb 			:: Vector{Float64}			# Lower Bounds of dimensions
	ub 			:: Vector{Float64}			# Upper Bounds of dimensions
	Binds 		:: Array{Int64,2} 			# Basis Function Indices for Smolyak Interpolant, f(D,mu)
	NumPts  	:: Int64					# Number of points in = Num Rows BF
	NumBF		:: Int64					# Number of basis functions under D, mu = Num Cols BF
	NumDeriv	:: Int64					# Number of derivatives: {0,1,2}
	max_order	:: Int64 					# Maximum order of polynomial for T
	x 			:: VecOrArray{Float64} 		# Vector of coordinates at which SB is evaluated 
	z 			:: VecOrArray{Float64} 		# Transformed vector of coordinates into [-1,1]
	T 			:: Array{Float64,2} 		# 1-dim Chebyshev basis fn: level
	dT 			:: Array{Float64,2} 		# 1-dim Chebyshev basis fn: 1st derivative
	d2T 		:: Array{Float64,2} 		# 1-dim Chebyshev basis fn: 2nd derivative
	BF 			:: VecOrArray{Float64}		# Basis Funs
	pinvBF		:: Array{Float64,2} 		# Inverse Basis Funs (Moore-Penrose Inverse if BF not square)
	dBFdz 		:: Array{Matrix{Float64},1} # 1st derivative basis funs wrt z
	d2BFdz2 	:: Array{Matrix{Float64},2}	# 2nd derivative basis funs wrt z
	dzdx		:: Vector{Float64} 			# Gradient of transform z→x
	d2zdx2		:: Array{Float64,2}			# Hessian of transform z→x
	dBFdx 		:: Array{Matrix{Float64},1} # 1st derivative basis funs wrt x
	d2BFdx2 	:: Array{Matrix{Float64},2}	# 2nd derivative basis funs wrt x
	CalcInv 	:: Bool 					# Whether or not to calculate BF^-1 				
	is_sg 		:: Bool						# INdicator of what was used to generate type

	# Constructor function with conformable memory allocations. Need to makeBF!(sb) to fill it in.
	function SmolyakBasis(sg::SmolyakGrid,NumDeriv::Int64=2,CalcInv::Bool=false,is_sg::Bool=true)
		
		# Components for evaluation of Basis Functions
		NumBF = size(sg.Binds,2)
		max_order = maximum(sg.Binds)+1
		T = Array(Float64, sg.D, max_order)
		dT = similar(T)
		d2T = similar(T)
		
		# For Basis Functions and transformation back: Allocate memory for BF, pinvBF, and derivatives -> then makeBF!(sb)
		BF = [1.0 for i1 in 1:NumBF, i2 in 1:sg.NumGrdPts]
		pinvBF = similar(BF)
		dBFdz = [1.0 for i1 in 1:NumBF, i2 in 1:sg.NumGrdPts, i3 in 1:sg.D] 
		d2BFdz2 = [1.0 for i1 in 1:NumBF, i2 in 1:sg.NumGrdPts, i3 in 1:sg.D, i3 in 1:sg.D]
		dzdx = 2./(sg.ub - sg.lb)
		d2zdx2 = dzdx*dzdx'
		dBFdx = similar(dBFdz)
		d2BFdx2 = similar(d2BFdz2)
		
		new(sg.D, sg.mu,  sg.lb, sg.ub, sg.Binds, 
			sg.NumGrdPts, NumBF, NumDeriv, max_order,
			sg.xGrid, sg.zGrid, T, dT, d2T, 
			BF, pinvBF, dBFdz, d2BFdz2, dzdx, d2zdx2, dBFdx, d2BFdx2,
			CalcInv, is_sg)
	end

	function SmolyakBasis(x::VecOrArray{Float64},sg::SmolyakGrid,NumDeriv::Int64=2,CalcInv::Bool=false,is_sg::Bool=false)
		
		z = x2z(x,sg.lb,sg.ub) #= x should be D x NumPts =#
		NumPts = size(x,2) 

		# Components for evaluation of Basis Functions
		NumBF = size(sg.Binds,2)
		max_order = maximum(sg.Binds)+1
		T = Array(Float64, sg.D, max_order)
		dT = similar(T)
		d2T = similar(T)
		
		# For Basis Functions and transformation back: Allocate memory for BF, pinvBF, and derivatives -> then makeBF!(sb)
		BF = [1.0 for i1 in 1:NumBF, i2 in 1:sg.NumGrdPts]
		pinvBF = similar(BF)
		dBFdz = [1.0 for i1 in 1:NumBF, i2 in 1:sg.NumGrdPts, i3 in 1:sg.D] 
		d2BFdz2 = [1.0 for i1 in 1:NumBF, i2 in 1:sg.NumGrdPts, i3 in 1:sg.D, i3 in 1:sg.D]
		dzdx = 2./(sg.ub - sg.lb)
		d2zdx2 = dzdx*dzdx'
		dBFdx = similar(dBFdz)
		d2BFdx2 = similar(d2BFdz2)
		
		new(sg.D, sg.mu,  sg.lb, sg.ub, sg.Binds, 
			NumPts, NumBF, NumDeriv, max_order,
			x, z, T, dT, d2T, 
			BF, pinvBF, dBFdz, d2BFdz2, dzdx, d2zdx2, dBFdx, d2BFdx2,
			CalcInv, is_sg)
	end 
	
end

# -------- New state to evaluate Smolyak Interpolant f(D,mu) where sb.use ------

# This is only for sb.is_sg=false
function new_x!(sb::SmolyakBasis,x::VecOrArray{Float64})
	for d in sb.x = copy(x)
	sb.z = x2z!(sb.x,sb.lb,sb.ub)
end

# ----------- Chebyshev Polynomials & derivatives ----------- #

# Basis Function evaluated at a vector z
function Tn!(sb::SmolyakBasis)
	for n in 1:sb.max_order, d in 1:sb.D
		if ==(n,1)
			sb.T[d,n] = 1.0
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 0.0 : nothing
			is(sb.NumDeriv,2) ? sb.d2T[d,n] = 0.0 : nothing				
		elseif ==(n,2)
			sb.T[d,n] = sb.z[d]
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 1.0 : nothing
			is(sb.NumDeriv,2) ? sb.d2T[d,n] = 0.0 : nothing				
		else
			sb.T[d,n] = 2*sb.z[d]*sb.T[d,n-1] - sb.T[d,n-2]
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 2*sb.T[d,n-1] + 2*sb.z[d]*sb.dT[d,n-1] - sb.dT[d,n-2] : nothing
			is(sb.NumDeriv,2) ? sb.d2T[d,n] = 4*sb.dT[d,n-1] + 2*sb.z[d]*sb.d2T[d,n-1] - sb.d2T[d,n-2] : nothing
		end
	end
end

# Basis Function evaluated at a vector z[:,i]
function Tn!(sb::SmolyakBasis,i::Int64=1)
	for n in 1:sb.max_order, d in 1:sb.D
		if ==(n,1)
			sb.T[d,n] = 1.0
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 0.0 : nothing
			is(sb.NumDeriv,2) ? sb.d2T[d,n] = 0.0 : nothing				
		elseif ==(n,2)
			sb.T[d,n] = sb.z[d,i]
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 1.0 : nothing
			is(sb.NumDeriv,2) ? sb.d2T[d,n] = 0.0 : nothing				
		else
			sb.T[d,n] = 2*sb.z[d,i]*sb.T[d,n-1] - sb.T[d,n-2]
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 2*sb.T[d,n-1] + 2*sb.z[d,i]*sb.dT[d,n-1] - sb.dT[d,n-2] : nothing
			is(sb.NumDeriv,2) ? sb.d2T[d,n] = 4*sb.dT[d,n-1] + 2*sb.z[d,i]*sb.d2T[d,n-1] - sb.d2T[d,n-2] : nothing
		end
	end
end

# ----------- BF & derivaitves ----------- # 

# Construct Basis Function
function BF!(sb::SmolyakBasis, BFIdx::Int64, DimIdx::Int64, GridIdx::Int64=1)
	sb.BF[BFIdx,GridIdx] *= sb.T[DimIdx,sb.Binds[DimIdx,BFIdx]+1]   
end

# Construct dBFdz! constructs 1st derivative of BF wrt z ∈ [-1,1], the transformed domain of state vector.  
function dBFdz!(sb::SmolyakBasis, BFIdx::Int64, DimIdx::Int64, GridIdx::Int64=1)
	for d in 1:sb.D
		if ==(d,DimIdx)
			sb.dBFdz[BFIdx,GridIdx,d] *= sb.dT[DimIdx,sb.Binds[DimIdx,BFIdx]+1 ]
		else
			sb.dBFdz[BFIdx,GridIdx,d] *= sb.T[DimIdx,sb.Binds[DimIdx,BFIdx]+1 ]
		end
	end
end

# Construct d2BFdz2! constructs 1st derivative of BF wrt z ∈ [-1,1], the transformed domain of state vector.  
function d2BFdz2!(sb::SmolyakBasis, BFIdx::Int64, DimIdx::Int64, GridIdx::Int64=1)
	for j in 1:sb.D, i in 1:sb.D
		if DimIdx==i==j
			sb.d2BFdz2[BFIdx,GridIdx,i,j] *= sb.d2T[DimIdx,sb.Binds[DimIdx,BFIdx]+1]
		elseif DimIdx==j!=i  
			sb.d2BFdz2[BFIdx,GridIdx,i,j] *= sb.dT[DimIdx,sb.Binds[DimIdx,BFIdx]+1]
		elseif DimIdx==i!=j
			sb.d2BFdz2[BFIdx,GridIdx,i,j] *= sb.dT[DimIdx,sb.Binds[DimIdx,BFIdx]+1]
		else
			sb.d2BFdz2[BFIdx,GridIdx,i,j] *= sb.T[DimIdx,sb.Binds[DimIdx,BFIdx]+1]
		end
	end
end

# Derivative of Transformations:dzdx, d2zdx2

#= Derivative constant over grid points under linear transform =#
function dzdx!(sb::SmolyakBasis)
	for d in 1:sb.D
		sb.dzdx[d] = 2/(sb.ub[d] - sb.lb[d])
	end
end

function d2zdx2!(sb::SmolyakBasis)
	sb.d2zdx2 = sb.dzdx*sb.dzdx'
end

# Derivatives of Basis Function with respect to state vector, x
function dBFdx!(sb::SmolyakBasis)
	for d in 1:sb.D
		sb.dBFdx[:,:,d] = sb.dBFdz[:,:,d]'*sb.dzdx[d]
	end
end

function d2BFdx2!(sb::SmolyakBasis)
	for j in 1:sb.D, i in 1:sb.D 
		sb.d2BFdx2[:,:,i,j] = sb.d2BFdz2[:,:,i,j]'*sb.d2zdx2[i,j]
	end
end

#= --------------------------------------- =#
#= Construct Basis Functions & Derivatives =#
#= --------------------------------------- =#

function check_domain(sb::SmolyakBasis, X::VecOrArray{Float64})
	@assert(is(size(X,1),sb.D ),
		"\n\tError: Dimension mismatch between Smolyak Grid and Matrix of Grid Points")	# Check correct number of dims
	# Check x is in bounds
	minX = Inf*ones(Float64,sb.D)
	maxX = zeros(Float64,sb.D)
	for n in 1:size(X,2), d in 1:size(X,1) 
		minX[d] = min(minX[d],X[d,n])
		maxX[d] = max(maxX[d],X[d,n])
	end
	for d in 1:sb.D 										# Report if X not in Bounds
		>(maxX[d],sb.ub[d]) ? 
			println("\tWarning: x[$d] > Upper Bound") : nothing
		<(minX[d],sb.lb[d]) ?
			println("\tWarning: x[$d] < Lower Bound") : nothing
	end
end

# This uses the Tn!(sb)
function makeBF!(sb::SmolyakBasis)
	if is(sb.NumDeriv,2)
		sb.BF = ones(Float64, sb.NumBF, sb.NumPts)
		sb.dBFdz = ones(Float64, sb.NumBF, sb.NumPts, sb.D) 
		sb.d2BFdz2 = ones(Float64, sb.NumBF, sb.NumPts, sb.D, sb.D)		
		for i in 1:sb.NumPts
			Tn!(sb,i)
			for d in 1:sb.D, p in 1:sb.NumBF
				BF!(sb, p, d, i)
				dBFdz!(sb, p, d, i)
				d2BFdz2!(sb, p, d, i)
			end
		end
		dBFdx!(sb)
		d2BFdx2!(sb)
	elseif is(sb.NumDeriv,1)
		sb.BF = ones(Float64, sb.NumBF, sb.NumPts)
		sb.dBFdz = ones(Float64, sb.NumBF, sb.NumPts, sb.D) 
		for i in 1:sb.NumPts
			Tn!(sb,i)
			for d in 1:sb.D, p in 1:sb.NumBF
				BF!(sb, p, d, i)
				dBFdz!(sb, p, d, i)
			end
		end
		dBFdx!(sb)
	elseif is(sb.NumDeriv,0)
		sb.BF = ones(Float64, sb.NumBF, sb.NumPts)
		for i in 1:sb.NumPts
			Tn!(sb,i)
			for d in 1:sb.D, p in 1:sb.NumBF
				BF!(sb, p, d, i)
			end
		end
	else
		print("Warning: sb.NumDeriv∈{0,1,2}")
	end
	is(sb.CalcInv,true) ? sb.pinvBF = pinv(sb.BF) : nothing
end

function show(io::IO, sb::SmolyakBasis)
	msg = "\n\tCreated Smolyak Basis:\n"
	msg *= "\t- Dim: $(sb.D), mu: $(sb.mu)\n"
	msg *= "\t- NumPts: $(sb.NumPts)\n"
	msg *= "\t- Number of Basis Functions: $(sb.NumBF)\n"
	if ==(sb.NumDeriv,0) msg *= "\t- No Derivative supplied. Do not call dBFdx or d2BFdx2.\n" end
	if ==(sb.NumDeriv,1) msg *= "\t- with dBFdx. Do not call d2BFdx2.\n" end
	if ==(sb.NumDeriv,2) msg *= "\t- with dBFdx & d2BFdx2.\n" end
	print(io, msg)
end
