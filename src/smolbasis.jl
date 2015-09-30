#= 
	-----------------------------------------------------------------
	Smolyak Basis Functions on a Smolyak Grid for Julia Version 0.4  
	-----------------------------------------------------------------

This file contains code to define Smolyak Constructtion of Chebyshev
basis functions. It is compatible with both Anisotrophic and Isotrophic 
Grids and are constructed efficiently following the methodology outlined
 in JMMV (2014). 

Key Refs: JMMV (2014), Burkhardt (2012)

=#

#= ******************* =#
#= Smolyak Basis type  =#
#= ******************* =#

type SmolyakBasis
	D 			:: Int64					# Dimensions
	mu 			:: ScalarOrVec{Int64}		# Index of mu
	lb 			:: Vector{Float64}			# Lower Bounds of dimensions
	ub 			:: Vector{Float64}			# Upper Bounds of dimensions
	Binds 		:: AA{Int64} 				# Basis Function Indices for Smolyak Interpolant, f(D,mu)
	NumPts  	:: Int64					# Number of points in = Num Rows BF
	NumBF		:: Int64					# Number of basis functions under D, mu = Num Cols BF
	NumDeriv	:: Int64					# Number of derivatives: {0,1,2}
	NumDerivArgs:: Int64					# 1st NumDerivArgs used (i.e. 1:NumDerivArgs)
	max_order	:: Int64 					# Maximum order of polynomial for T
	x 			:: VecOrAA{Float64} 		# Vector of coordinates at which SB is evaluated 
	z 			:: VecOrAA{Float64} 		# Transformed vector of coordinates into [-1,1]
	T 			:: Matrix{Float64} 			# 1-dim Chebyshev basis fn: level
	dT 			:: Matrix{Float64} 			# 1-dim Chebyshev basis fn: 1st derivative
	d2T 		:: Matrix{Float64} 			# 1-dim Chebyshev basis fn: 2nd derivative
	BF 			:: AA{Float64}				# Basis Funs
	dBFdz 		:: AAA{Float64} 			# 1st derivative basis funs wrt z
	d2BFdz2 	:: AAAA{Float64}			# 2nd derivative basis funs wrt z
	dzdx		:: Vector{Float64} 			# Gradient of transform z→x
	d2zdx2		:: Vector{Float64}			# Diagonal of Hessian of transform z→x (0 in linear maps, so only useful in nonlinear mapping)
	dBFdx 		:: AAA{Float64} 			# 1st derivative basis funs wrt x
	d2BFdx2 	:: AAAA{Float64}			# 2nd derivative basis funs wrt x

	# Constructor function with conformable memory allocations. Need to makeBF!(sb) to fill it in.
	function SmolyakBasis(sg::SmolyakGrid,NumDeriv::Int64=2,NumDerivArgs::Int64=sg.D)
		
		# Components for evaluation of Basis Functions
		NumBF = length(sg.Binds)
		max_order = 0
		for i in eachindex(sg.Binds), j in 1:sg.D
			max_order = max(max_order,sg.Binds[i][j]+1)
		end

		# Preallocate memory for basis functions
		T = Array{Float64}(sg.D, max_order) 	
		dT = similar(T)
		d2T = similar(T)
		
		# For Basis Functions and transformation back: Allocate memory for BF, pinvBF, and derivatives -> then makeBF!(sb)
		BF = Vector{Float64}[ones(Float64,NumBF) 
					for n in 1:sg.NumGrdPts]			# BF[n][p] where n =1:NumGrdPts, p=1:NumBF
		
		dBFdz = AA{Float64}[[ones(Float64,NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:sg.NumGrdPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		d2BFdz2 = AAA{Float64}[[[ones(Float64,NumBF) 
					for i in k:NumDerivArgs]
					for k in 1:NumDerivArgs]
					for n in 1:sg.NumGrdPts] 			#= d2BFdz2[n][i][k][p] where n =1:NumGrdPts,
														 i is position of 1st derivative, 
														 k = j - i + 1 where j in position of 2nd derivative, and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:NumDerivArgs
			push!(dzdx,2/(sg.ub[i] - sg.lb[i]))
		end
		d2zdx2 = zeros(Float64,NumDerivArgs) 
		
		dBFdx = AA{Float64}[[ones(Float64,NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:sg.NumGrdPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		d2BFdx2 = AAA{Float64}[[[ones(Float64,NumBF) 
					for i in k:NumDerivArgs]
					for k in 1:NumDerivArgs]
					for n in 1:sg.NumGrdPts] 			#= d2BFdz2[n][i][k][p] where n =1:NumGrdPts,
														i is position of 1st derivative, 
														k = j - i + 1 where j in position of 2nd derivative, and p=1:NumBF =#	
		new(sg.D, sg.mu,  sg.lb, sg.ub, sg.Binds, 
			sg.NumGrdPts, NumBF, NumDeriv, NumDerivArgs, max_order,
			sg.xGrid, sg.zGrid, T, dT, d2T, 
			BF, dBFdz, d2BFdz2, dzdx, d2zdx2, dBFdx, d2BFdx2)
	end

	function SmolyakBasis(x::VecOrAA{Float64},sg::SmolyakGrid,NumDeriv::Int64=2,NumDerivArgs::Int64=sg.D)
		
		z = x2z(x,sg.lb,sg.ub) #= x should be D x NumPts =#
		NumPts = size(x,2) 

		# Components for evaluation of Basis Functions
		NumBF = length(sg.Binds)
		max_order = 0
		for i in eachindex(sg.Binds), j in 1:sg.D
			max_order = max(max_order,sg.Binds[i][j]+1)
		end

		# Preallocate memory for basis functions
		T = Array{Float64}(sg.D, max_order) 	
		dT = similar(T)
		d2T = similar(T)

		# For Basis Functions and transformation back: Allocate memory for BF, pinvBF, and derivatives -> then makeBF!(sb)
		BF = Vector{Float64}[ones(Float64,NumBF) 
					for n in 1:NumPts]			# BF[n][p] where n =1:NumGrdPts, p=1:NumBF
		
		dBFdz = AA{Float64}[[ones(Float64,NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		d2BFdz2 = AAA{Float64}[[[ones(Float64,NumBF) 
					for i in k:NumDerivArgs]
					for k in 1:NumDerivArgs]
					for n in 1:NumPts] 			#= d2BFdz2[n][i][k][p] where n =1:NumGrdPts,
														 i is position of 1st derivative, 
														 k = j - i + 1 where j in position of 2nd derivative, and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:NumDerivArgs
			push!(dzdx,2/(sg.ub[i] - sg.lb[i]))
		end
		d2zdx2 = zeros(Float64,NumDerivArgs) 
		
		dBFdx = AA{Float64}[[ones(Float64,NumBF) 
					for i in 1:sg.D] 
					for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		d2BFdx2 = AAA{Float64}[[[ones(Float64,NumBF) 
					for i in k:NumDerivArgs]
					for k in 1:NumDerivArgs]
					for n in 1:NumPts] 			#= d2BFdz2[n][i][k][p] where n =1:NumGrdPts,
														i is position of 1st derivative, 
														k = j - i + 1 where j in position of 2nd derivative, and p=1:NumBF =#
		new(sg.D, sg.mu,  sg.lb, sg.ub, sg.Binds, 
			NumPts, NumBF, NumDeriv, NumDerivArgs, max_order,
			x, z, T, dT, d2T, 
			BF, dBFdz, d2BFdz2, dzdx, d2zdx2, dBFdx, d2BFdx2)
	end 
	
end

# -------- New state to evaluate Smolyak Interpolant f(D,mu) where sb.use ------

# New state vector -> applicable when sb.x is a vector
function new_x!(sb::SmolyakBasis,x::Vector{Float64})
	for d in 1:sb.D
		sb.x[d] = x[d]
	end
	x2z!(sb.x,sb.z,sb.lb,sb.ub)
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

# Basis Function evaluated at a vector z at a grid point i = 1:NumPts
function Tn!(sb::SmolyakBasis,i::Int64=1)
	for n in 1:sb.max_order, d in 1:sb.D
		if ==(n,1)
			sb.T[d,n] = 1.0
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 0.0 : nothing
			is(sb.NumDeriv,2) ? sb.d2T[d,n] = 0.0 : nothing				
		elseif ==(n,2)
			sb.T[d,n] = sb.z[i][d]
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 1.0 : nothing
			is(sb.NumDeriv,2) ? sb.d2T[d,n] = 0.0 : nothing				
		else
			sb.T[d,n] = 2*sb.z[i][d]*sb.T[d,n-1] - sb.T[d,n-2]
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 2*sb.T[d,n-1] + 2*sb.z[i][d]*sb.dT[d,n-1] - sb.dT[d,n-2] : nothing
			is(sb.NumDeriv,2) ? sb.d2T[d,n] = 4*sb.dT[d,n-1] + 2*sb.z[i][d]*sb.d2T[d,n-1] - sb.d2T[d,n-2] : nothing
		end
	end
end


# ----------- BF & derivaitves ----------- # 

# Construct Basis Function
function BF!(sb::SmolyakBasis, BFIdx::Int64, DimIdx::Int64, GridIdx::Int64=1)
	sb.BF[GridIdx][BFIdx] *= sb.T[DimIdx,sb.Binds[BFIdx][DimIdx]+1]   
end

# Construct dBFdz! constructs 1st derivative of BF wrt z ∈ [-1,1], the transformed domain of state vector.   for first n arguments of state vector 
function dBFdz!(sb::SmolyakBasis, BFIdx::Int64, DimIdx::Int64, GridIdx::Int64=1, N::Int64=sb.NumDerivArgs)
	for d in 1:N 		# n specifies 1:n derivatives to avoid unnecessary computations is derivatives of first n arguments required
		if ==(d,DimIdx)
			sb.dBFdz[GridIdx][d][BFIdx] *= sb.dT[DimIdx,sb.Binds[BFIdx][DimIdx]+1 ]
		else
			sb.dBFdz[GridIdx][d][BFIdx] *= sb.T[DimIdx,sb.Binds[BFIdx][DimIdx]+1 ]
		end
	end
end

# Construct d2BFdz2! constructs 1st derivative of BF wrt z ∈ [-1,1], the transformed domain of state vector for first n arguments of state vector 
function d2BFdz2!(sb::SmolyakBasis, BFIdx::Int64, DimIdx::Int64, GridIdx::Int64=1, N::Int64=sb.NumDerivArgs)
	for i in 1:N, j in i:N
		k = j-i+1 				# Translate position of argument to index used to access the jth argument's derivative.
		if DimIdx==i==j
			sb.d2BFdz2[GridIdx][i][k][BFIdx] *= sb.d2T[DimIdx,sb.Binds[BFIdx][DimIdx]+1]
		elseif DimIdx==j!=i  
			sb.d2BFdz2[GridIdx][i][k][BFIdx] *= sb.dT[DimIdx,sb.Binds[BFIdx][DimIdx]+1]
		elseif DimIdx==i!=j
			sb.d2BFdz2[GridIdx][i][k][BFIdx] *= sb.dT[DimIdx,sb.Binds[BFIdx][DimIdx]+1]
		else
			sb.d2BFdz2[GridIdx][i][k][BFIdx] *= sb.T[DimIdx,sb.Binds[BFIdx][DimIdx]+1]
		end
	end
end

# Derivative of Transformations:dzdx, d2zdx2

#= Derivative constant over grid points under linear transform =#
function dzdx!(sb::SmolyakBasis,N::Int64=sb.NumDerivArgs)
	for d in 1:N
		sb.dzdx[d] = 2/(sb.ub[d] - sb.lb[d])  
	end
end

# 1st derivatives for first n arguments of state vector
function dBFdx!(sb::SmolyakBasis,N::Int64=sb.NumDerivArgs)
	for n in 1:sb.NumPts, d in 1:N, p in 1:sb.NumBF
		sb.dBFdx[n][d][p] = sb.dBFdz[n][d][p]*sb.dzdx[d] # matrix x scalar
	end
end

# Second derivatives for first n arguments of state vector
function d2BFdx2!(sb::SmolyakBasis,N::Int64=sb.NumDerivArgs)
	for n in 1:sb.NumPts, i in 1:N, j in i:N, p in 1:sb.NumBF
		k = j-i+1 
		sb.d2BFdx2[n][i][k][p] = sb.d2BFdz2[n][i][k][p]*sb.dzdx[i]*sb.dzdx[j] # matrix x scalar x scalar
	end
end

#= --------------------------------------- =#
#= Construct Basis Functions & Derivatives =#
#= --------------------------------------- =#

# Initialise Basis Functions with 1's
function initBF!(sb::SmolyakBasis,N::Int64=sb.NumDerivArgs)
	if is(sb.NumDeriv,2)
		sb.BF = Vector{Float64}[ones(Float64,sb.NumBF) 
				for n in 1:sb.NumPts]			# BF[n][p] where n =1:NumGrdPts, p=1:NumBF		
		sb.dBFdz = AA{Float64}[[ones(Float64,sb.NumBF) 
				for i in 1:N] 
				for n in 1:sb.NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, i is position of 1st derivative,and p=1:NumBF =#
		sb.d2BFdz2 = AAA{Float64}[[[ones(Float64,sb.NumBF) 
				for i in k:N]
				for k in 1:N]
				for n in 1:sb.NumPts] 			#= d2BFdz2[n][i][k][p] where n =1:NumGrdPts, i is position of 1st derivative, k = j - i + 1 where j in position of 2nd derivative, and p=1:NumBF =#
	elseif is(sb.NumDeriv,1)
				sb.BF = Vector{Float64}[ones(Float64,sb.NumBF) 
				for n in 1:sb.NumPts]			# BF[n][p] where n =1:NumGrdPts, p=1:NumBF		
		sb.dBFdz = AA{Float64}[[ones(Float64,sb.NumBF) 
				for i in 1:N] 
				for n in 1:sb.NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, i is position of 1st derivative,and p=1:NumBF =#
	else
		sb.BF = Vector{Float64}[ones(Float64,sb.NumBF) 
				for n in 1:sb.NumPts]			# BF[n][p] where n =1:NumGrdPts, p=1:NumBF
	end
end


# Makes Basis Functions with sb.NumDeriv derivatives of the first n arguments of state vector
function makeBF!(sb::SmolyakBasis,N::Int64=sb.NumDerivArgs)
	isa(sb.x,Vector{Float64}) ? sb.NumPts = length(sb.x) : nothing
	initBF!(sb,N) # Need to start with BF and derivatives as 1 because will take product over loops
	if is(sb.NumDeriv,2)
		for i in 1:sb.NumPts
			Tn!(sb,i)
			for d in 1:N, p in 1:sb.NumBF
				BF!(sb, p, d, i)
				dBFdz!(sb, p, d, i, N)
				d2BFdz2!(sb, p, d, i, N) 			# Hess
			end
		end
		dBFdx!(sb, N)
		d2BFdx2!(sb, N)
	elseif is(sb.NumDeriv,1)
		for i in 1:sb.NumPts
			Tn!(sb,i)
			for d in 1:N, p in 1:sb.NumBF
				BF!(sb, p, d, i)
				dBFdz!(sb, p, d, i, N) 				# Jac
			end
		end
		dBFdx!(sb, N)
	elseif is(sb.NumDeriv,0)
		for i in 1:sb.NumPts 						# Loop over each grid point.
			Tn!(sb,i) 								# Evaluate Basis Function Value in each Dimension at each grid point, i: This is D x NumBF 
			for d in 1:n, p in 1:sb.NumBF 			# Make Basis Function levels & derivatives for grid points i:
				BF!(sb, p, d, i) 					# Now need to multiply over dimension to create a NBF-vector of Basis Function at grid point i. END OF LOOP → NumBF x NumPts for all grid points. 
			end 									
		end
	else
		print("Warning: sb.NumDeriv∈{0,1,2}")
	end
end

function show(io::IO, sb::SmolyakBasis)
	msg = "\n\tCreated Smolyak Basis:\n"
	msg *= "\t- Dim: $(sb.D), mu: $(sb.mu)\n"
	msg *= "\t- NumPts: $(sb.NumPts)\n"
	msg *= "\t- Number of Basis Functions: $(sb.NumBF)\n"
	if ==(sb.NumDeriv,0) msg *= "\t- No Derivative supplied. Do not call dBFdx or d2BFdx2.\n" end
	if ==(sb.NumDeriv,1) msg *= "\t- with dBFdx up to first $(sb.NumDerivArgs) arguments. Do not call d2BFdx2.\n" end
	if ==(sb.NumDeriv,2) msg *= "\t- with dBFdx & d2BFdx2 up to first $(sb.NumDerivArgs) arguments.\n" end
	print(io, msg)
end
