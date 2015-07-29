#= 
	-----------------------------------------------------------------
	Smolyak Basis Functions on a Smolyak Grid for Julia version 0.3.7  
	-----------------------------------------------------------------

This file contains code to define Smolyak Constructtion of Chebyshev
basis functions. It is compatible with both AnisotroΨc and IsotroΨc 
Grids and are constructed efficiently following the methodology outlined
 in JMMV (2014). The code is designed on latest stable Julia version: 0.3.7.

Key Refs: JMMV (2014), Burkhardt (2012), github: ECONFORGE/Smolyak

=#

#= ---------------- =#
#= Within loop funs =#
#= ---------------- =#

VecOrArray = Union(Vector{Float64},Array{Float64,2}) 

function Tn!(T::Array{Float64,2},zpts::VecOrArray, MaxChebIdx::Int64, Dims::Int64)
	for n in 1:MaxChebIdx, d in 1:Dims
		if ==(n,1)
			T[d,n] = 1.0
		elseif ==(n,2)
			T[d,n] = zpts[d]
		else
			T[d,n] = 2*zpts[d]*T[d,n-1] - T[d,n-2]
		end
	end
end


function Tn!(T::Array{Float64,2}, ∂T::Array{Float64,2}, zpts::VecOrArray, MaxChebIdx::Int64, Dims::Int64)
	for n in 1:MaxChebIdx, d in 1:Dims
		if ==(n,1)
			T[d,n] = 1.0
			∂T[d,n] = 0.0
		elseif ==(n,2)
			T[d,n] = zpts[d]
			∂T[d,n] = 1.0
		else
			T[d,n] = 2*zpts[d]*T[d,n-1] - T[d,n-2]
			∂T[d,n] = 2*T[d,n-1] + 2*zpts[d]*∂T[d,n-1] - ∂T[d,n-2]
		end
	end
end

function Tn!(T::Array{Float64,2}, ∂T::Array{Float64,2}, ∂2T::Array{Float64,2}, zpts::VecOrArray, MaxChebIdx::Int64, Dims::Int64)
	for n in 1:MaxChebIdx, d in 1:Dims
		if ==(n,1)
			T[d,n] = 1.0
			∂T[d,n] = 0.0
			∂2T[d,n] = 0.0				
		elseif ==(n,2)
			T[d,n] = zpts[d]
			∂T[d,n] = 1.0
			∂2T[d,n] = 0.0				
		else
			T[d,n] = 2*zpts[d]*T[d,n-1] - T[d,n-2]
			∂T[d,n] = 2*T[d,n-1] + 2*zpts[d]*∂T[d,n-1] - ∂T[d,n-2]
			∂2T[d,n] = 4*∂T[d,n-1] + 2*zpts[d]*∂2T[d,n-1] - ∂2T[d,n-2]
		end
	end
end

function Ψ!(Ψ::Vector{Float64}, T::Array{Float64,2}, BFIdx::Int64, DimIdx::Int64, sg::SmolyakGrid)
	ChebIdx = sg.Binds[DimIdx,BFIdx]+1 
	Ψ[BFIdx] *= T[DimIdx,ChebIdx] 
end

function Ψ!(Ψ::Array{Float64,2}, T::Array{Float64,2}, 
				GridIdx::Int64, BFIdx::Int64, DimIdx::Int64, sg::SmolyakGrid)
	ChebIdx = sg.Binds[DimIdx,BFIdx]+1 
	Ψ[BFIdx,GridIdx] *= T[DimIdx,ChebIdx] 
end


function ∂Ψ∂z!(∂Ψ∂z::Array{Float64,3}, T::Array{Float64,2}, ∂T::Array{Float64,2}, 
					GridIdx::Int64, BFIdx::Int64, DimIdx::Int64, sg::SmolyakGrid)
	ChebIdx = sg.Binds[DimIdx,BFIdx]+1 
	for i in 1:sg.D
		if ==(i,DimIdx)
			∂Ψ∂z[BFIdx,GridIdx,i] *= ∂T[DimIdx,ChebIdx]
		else
			∂Ψ∂z[BFIdx,GridIdx,i] *= T[DimIdx,ChebIdx]
		end
	end
end

function ∂2Ψ∂z2!(∂2Ψ∂z2::Array{Float64,4}, T::Array{Float64,2}, ∂T::Array{Float64,2}, ∂2T::Array{Float64,2},
					 GridIdx::Int64, BFIdx::Int64, DimIdx::Int64, sg::SmolyakGrid)
	ChebIdx = sg.Binds[DimIdx,BFIdx]+1
	for j in 1:sg.D, i in 1:sg.D
		if DimIdx==i==j
			∂2Ψ∂z2[BFIdx,GridIdx,i,j] *= ∂2T[DimIdx,ChebIdx]
		elseif DimIdx==j!=i  
			∂2Ψ∂z2[BFIdx,GridIdx,i,j] *= ∂T[DimIdx,ChebIdx]
		elseif DimIdx==i!=j
			∂2Ψ∂z2[BFIdx,GridIdx,i,j] *= ∂T[DimIdx,ChebIdx]
		else
			∂2Ψ∂z2[BFIdx,GridIdx,i,j] *= T[DimIdx,ChebIdx]
		end
	end
end

#= --------------------------------------- =#
#= Construct Basis Functions & Derivatives =#
#= --------------------------------------- =#

function check_domain(X::Array{Float64,2}, sg::SmolyakGrid)
	D, N = size(X) 	# Size of Array of Grid points
	@assert(is(D,sg.D ),
		"\n\tError: Dimension mismatch between 
			Smolyak Grid and Matrix of Grid Points")	# Check correct number of dims
	# Check x is in bounds
	minX = Inf*ones(Float64,D)
	maxX = zeros(Float64,D)
	for i in 1:N, d in 1:D 
		minX[d] = min(minX[d],X[d,i])
		maxX[d] = max(maxX[d],X[d,i])
	end
	for d in 1:D 										# Report if X not in Bounds
		>(maxX[d],sg.ub[d]) ? 
			println("\tWarning: X[$d] > Upper Bound") : nothing
		<(minX[d],sg.lb[d]) ?
			println("\tWarning: X[$d] < Lower Bound") : nothing
	end
end

function Ψ_fun(sg::SmolyakGrid, k::Int64, SpOut::Bool)
	NBF = size(sg.Binds,2)
	M = maximum(sg.Binds)+1
	if ==(k,2)
		Ψ = ones(Float64, NBF, sg.NumGrdPts)
		∂Ψ∂z = ones(Float64, NBF, sg.NumGrdPts, sg.D)
		∂2Ψ∂z2 = ones(Float64, NBF, sg.NumGrdPts, sg.D, sg.D)
		T = Array(Float64, sg.D, M)
		∂T = similar(T)
		∂2T = similar(T)
		for i in 1:sg.NumGrdPts
			Tn!(T, ∂T, ∂2T, sg.zGrid[:,i], M, sg.D)
			for d in 1:sg.D, p in 1:NBF
				Ψ!(Ψ, T, i, p, d, sg)
				∂Ψ∂z!(∂Ψ∂z, T, ∂T, i, p, d, sg)
				∂2Ψ∂z2!(∂2Ψ∂z2, T, ∂T, ∂2T, i, p, d, sg)
			end
		end
		if is(SpOut,true) 								# Return Sparse Basis Matrix
			Z1idx = sub2ind(size(∂Ψ∂z),findn(∂Ψ∂z)...)
			Z2idx = sub2ind(size(∂2Ψ∂z2),findn(∂2Ψ∂z2)...) 
			return Ψ, sparsevec(Z1idx,∂Ψ∂z[Z1idx]), sparsevec(Z2idx,∂2Ψ∂z2[Z2idx])
		else										# Return Full Basis Matrix
			return Ψ, ∂Ψ∂z, ∂2Ψ∂z2
		end
	elseif ==(k,1)
		Ψ = ones(Float64, NBF, sg.NumGrdPts)
		∂Ψ∂z = ones(Float64, NBF, sg.NumGrdPts, sg.D)
		T = Array(Float64, sg.D, M)
		∂T = similar(T)
		for i in 1:sg.NumGrdPts
			Tn!(T, ∂T, sg.zGrid[:,i], M, sg.D)
			for d in 1:sg.D, p in 1:NBF
				Ψ!(Ψ, T, i, p, d, sg)
				∂Ψ∂z!(∂Ψ∂z, T, ∂T, i, p, d, sg)
			end
		end
		if is(SpOut,true) 								# Return Sparse Basis Matrix
			Z1idx = sub2ind(size(∂Ψ∂z),findn(∂Ψ∂z)...)
			return Ψ, sparsevec(Z1idx,∂Ψ∂z[Z1idx]), Array(Float64,1,1,1,1)
		else										# Return Full Basis Matrix
			return Ψ, ∂Ψ∂z, Array(Float64,1,1,1,1)
		end
	elseif ==(k,0) 
		Ψ = ones(Float64, NBF, sg.NumGrdPts)
		T = Array(Float64, sg.D, M)
		for i in 1:sg.NumGrdPts
			Tn!(T, sg.zGrid[:,i], M, sg.D)
			for d in 1:sg.D, p in 1:NBF
				Ψ!(Ψ, T, i, p, d, sg)
			end
		end
		return Ψ, Array(Float64,1,1,1), Array(Float64,1,1,1,1)
	else
		print("Warning: You µst specify number of derivatives to be 0, 1, or 2")
	end
end


function Ψ_fun(X::Vector{Float64}, sg::SmolyakGrid, k::Int64, SpOut::Bool)

	# Convert to z in [-1,1]
	#check_domain(X,sg)
	zX = x2z(X,sg.lb,sg.ub)					# Convert to z in [-1,1]

	D = length(X) 							# Size of Array of Grid points
	NBF = size(sg.Binds,2)
	M = maximum(sg.Binds)+1
	if ==(k,2)
		Ψ = ones(Float64,NBF)
		∂Ψ∂z = ones(Float64, NBF, 1, sg.D)
		∂2Ψ∂z2 = ones(Float64, NBF, 1, sg.D, sg.D)
		T = Array(Float64, sg.D, M)
		∂T = similar(T)
		∂2T = similar(T)	
		Tn!(T, ∂T, ∂2T, zX, M, sg.D)
		for d in 1:sg.D, p in 1:NBF
			Ψ!(Ψ, T, p, d, sg)
			∂Ψ∂z!(∂Ψ∂z, T, ∂T, 1, p, d, sg)
			∂2Ψ∂z2!(∂2Ψ∂z2, T, ∂T, ∂2T, 1, p, d, sg)
		end
		if is(SpOut,true) 								# Return Sparse Basis Matrix
			Z1idx = sub2ind(size(∂Ψ∂z),findn(∂Ψ∂z)...)
			Z2idx = sub2ind(size(∂2Ψ∂z2),findn(∂2Ψ∂z2)...) 
			return Ψ, sparsevec(Z1idx,∂Ψ∂z[Z1idx]), sparsevec(Z2idx,∂2Ψ∂z2[Z2idx])
		else										# Return Full Basis Matrix
			return Ψ, ∂Ψ∂z, ∂2Ψ∂z2
		end
	elseif ==(k,1)
		Ψ = ones(Float64,  NBF)
		∂Ψ∂z = ones(Float64, NBF, 1, sg.D)
		T = Array(Float64, sg.D, M)
		∂T = similar(T)
		Tn!(T, ∂T, zX, M, sg.D)
		for d in 1:sg.D, p in 1:NBF
			Ψ!(Ψ, T, p, d, sg)
			∂Ψ∂z!(∂Ψ∂z, T, ∂T, 1, p, d, sg)
		end
		if is(SpOut,true) 								# Return Sparse Basis Matrix
			Z1idx = sub2ind(size(∂Ψ∂z),findn(∂Ψ∂z)...)
			return Ψ, sparsevec(Z1idx,∂Ψ∂z[Z1idx]), Array(Float64,1,1,1,1)
		else										# Return Full Basis Matrix
			return Ψ, ∂Ψ∂z, Array(Float64,1,1,1,1)
		end
	elseif ==(k,0) 
		Ψ = ones(Float64,NBF)
		T = Array(Float64,sg.D,M)
		Tn!(T, zX, M, sg.D)
		for d in 1:sg.D, p in 1:NBF
			Ψ!(Ψ, T, p, d, sg)
		end
		return Ψ, Array(Float64,1,1,1), Array(Float64,1,1,1,1)
	else
		print("Warning: You µst specify number of derivatives to be 0, 1, or 2")
	end
end

function Ψ_fun(X::Array{Float64,2}, sg::SmolyakGrid, k::Int64, SpOut::Bool)

	# Convert to z in [-1,1]
	#check_domain(X,sg)
	zX = x2z(X,sg.lb,sg.ub)					# Convert to z in [-1,1]

	D, N = size(X) 							# Size of Array of Grid points
	NBF = size(sg.Binds,2)
	M = maximum(sg.Binds)+1
	if ==(k,2)
		Ψ = ones(Float64,NBF, N)
		∂Ψ∂z = ones(Float64, NBF, N,sg.D)
		∂2Ψ∂z2 = ones(Float64, NBF, N, sg.D, sg.D)
		T = Array(Float64, sg.D, M)
		∂T = similar(T)
		∂2T = similar(T)
		for i in 1:N
			Tn!(T, ∂T, ∂2T, zX[:,i], M, sg.D)
			for d in 1:sg.D, p in 1:NBF
				Ψ!(Ψ, T, i, p, d, sg)
				∂Ψ∂z!(∂Ψ∂z, T, ∂T, i, p, d, sg)
				∂2Ψ∂z2!(∂2Ψ∂z2, T, ∂T, ∂2T, i, p, d, sg)
			end
		end
		if is(SpOut,true) 								# Return Sparse Basis Matrix
			Z1idx = sub2ind(size(∂Ψ∂z),findn(∂Ψ∂z)...)
			Z2idx = sub2ind(size(∂2Ψ∂z2),findn(∂2Ψ∂z2)...) 
			return Ψ, sparsevec(Z1idx,∂Ψ∂z[Z1idx]), sparsevec(Z2idx,∂2Ψ∂z2[Z2idx])
		else										# Return Full Basis Matrix
			return Ψ, ∂Ψ∂z, ∂2Ψ∂z2
		end
	elseif ==(k,1)
		Ψ = ones(Float64,  NBF, N)
		∂Ψ∂z = ones(Float64, NBF, N, sg.D)
		T = Array(Float64, sg.D, M)
		∂T = similar(T)
		for i in 1:N
			Tn!(T, ∂T, zX[:,i], M, sg.D)
			for d in 1:sg.D, p in 1:NBF
				Ψ!(Ψ, T, i, p, d, sg)
				∂Ψ∂z!(∂Ψ∂z, T, ∂T, i, p, d, sg)
			end
		end
		if is(SpOut,true) 								# Return Sparse Basis Matrix
			Z1idx = sub2ind(size(∂Ψ∂z),findn(∂Ψ∂z)...)
			return Ψ, sparsevec(Z1idx,∂Ψ∂z[Z1idx]), Array(Float64,1,1,1,1)
		else										# Return Full Basis Matrix
			return Ψ, ∂Ψ∂z, Array(Float64,1,1,1,1)
		end
	elseif ==(k,0) 
		Ψ = ones(Float64,NBF, N)
		T = Array(Float64,sg.D,M)
		for i in 1:N
			Tn!(T, zX[:,i], M, sg.D)
			for d in 1:sg.D, p in 1:NBF
				Ψ!(Ψ, T, i, p, d, sg)
			end
		end
		return Ψ, Array(Float64,1,1,1), Array(Float64,1,1,1,1)
	else
		print("Warning: You µst specify number of derivatives to be 0, 1, or 2")
	end
end



#= Derivative constant over grid points under linear transform =#
function ∂z∂x_fun(sg::SmolyakGrid)
	∂z∂x = 2./(sg.ub - sg.lb)
	∂2z∂x2 = ∂z∂x*∂z∂x'
	return ∂z∂x, ∂2z∂x2 
end



#= ************** =#
#= PolyBasis type =#
#= ************** =#

IntOrVec = Union(Int64,Vector{Int64},Float64,Vector{Float64})
SpArray = Union(SparseMatrixCSC,Array)

type SmolyakBasis
	D 			:: Int64				# Dimensions
	µ 			:: IntOrVec				# Index of µ
	NumPts  	:: Int64				# Number of points in = Num Rows Ψ
	NumBasisFun	:: Int64				# Number of basis functions under D, µ = Num Cols Ψ
	Ψ 			:: VecOrArray 			# Basis Funs
	pinvΨ		:: Array{Float64,2} 	# Inverse Basis Funs
	∂Ψ∂z 		:: SpArray 				# 1st derivative basis funs
	∂2Ψ∂z2 		:: SpArray			 	# 2nd derivative basis funs
	∂z∂x		:: Vector{Float64} 		# Gradient of transform z2x()
	∂2z∂x2		:: Array{Float64,2}		# Hessian of transform z2x()
	CalcInv 	:: Bool 				# Whether or not to calculate Ψ^-1 				
	NumDeriv	:: Int64				# Number of derivatives: {0,1,2}
	SpOut 		:: Bool					# Sparse Output indicator == 1 if Sparse, 0 Otherwise

	function SmolyakBasis(sg::SmolyakGrid, CalcInv::Bool=false, NumDeriv::Int64=1, SpOut::Bool=false)
		Ψ, ∂Ψ∂z, ∂2Ψ∂z2 = Ψ_fun(sg, NumDeriv, SpOut)
		NBF, NGP = size(Ψ)
		is(CalcInv,true) ? invΨ = inv(Ψ) : invΨ = Array(Float64,1,1)
		if >(NumDeriv,0)
			∂z∂x, ∂2z∂x2 = ∂z∂x_fun(sg) 
		else 
			∂z∂x = Array(Float64,1)
			∂2z∂x2 = Array(Float64,1,1) 
		end
		new(sg.D, sg.µ, NGP, NBF, Ψ, invΨ, ∂Ψ∂z, ∂2Ψ∂z2, ∂z∂x, ∂2z∂x2, CalcInv, NumDeriv, SpOut)
	end
	
	function SmolyakBasis(X::Vector{Float64}, sg::SmolyakGrid, CalcInv::Bool=false, NumDeriv::Int64=1, SpOut::Bool=false)
		Ψ, ∂Ψ∂z, ∂2Ψ∂z2 = Ψ_fun(X, sg, NumDeriv, SpOut)
		NBF = length(Ψ)
		NGP = 1
		is(CalcInv,true) ? invΨ = pinv(Ψ) : invΨ = Array(Float64,1,1)	
		if >(NumDeriv,0)
			∂z∂x, ∂2z∂x2 = ∂z∂x_fun(sg) 
		else 
			∂z∂x = Array(Float64,1)
			∂2z∂x2 = Array(Float64,1,1) 
		end		
		new(sg.D, sg.µ, NGP, NBF, Ψ, invΨ, ∂Ψ∂z, ∂2Ψ∂z2, ∂z∂x, ∂2z∂x2, CalcInv, NumDeriv, SpOut )	
	end

	function SmolyakBasis(X::Array{Float64,2}, sg::SmolyakGrid, CalcInv::Bool=false, NumDeriv::Int64=1, SpOut::Bool=false)
		Ψ, ∂Ψ∂z, ∂2Ψ∂z2 = Ψ_fun(X, sg, NumDeriv, SpOut)
		NBF, NGP = size(Ψ)
		is(CalcInv,true) ? invΨ = pinv(Ψ) : invΨ = Array(Float64,1,1)	
		if >(NumDeriv,0)
			∂z∂x, ∂2z∂x2 = ∂z∂x_fun(sg) 
		else 
			∂z∂x = Array(Float64,1)
			∂2z∂x2 = Array(Float64,1,1) 
		end		
		new(sg.D, sg.µ, NGP, NBF, Ψ, invΨ, ∂Ψ∂z, ∂2Ψ∂z2, ∂z∂x, ∂2z∂x2, CalcInv, NumDeriv, SpOut )	
	end
	
end

function show(io::IO, sb::SmolyakBasis)
	msg = "\n\tCreated Smolyak Basis:\n"
	msg *= "\t- Dim: $(sb.D), µ: $(sb.µ)\n"
	msg *= "\t- NumPts: $(sb.NumPts)\n"
	msg *= "\t- Number of Basis Functions: $(sb.NumBasisFun)\n"
	if ==(sb.NumDeriv,0) msg *= "\t- No Derivative supplied. Do not call ∂Ψ∂z or ∂2Ψ∂z2!\n" end
	if ==(sb.NumDeriv,1) msg *= "\t- with ∂Ψ∂z. Do not call ∂2Ψ∂z2.\n" end
	if ==(sb.NumDeriv,2) msg *= "\t- with ∂Ψ∂z & ∂2Ψ∂z2\n" end
	print(io, msg)
end

function sparse2full(sb::SmolyakBasis)
	# Automatically Return NumDeriv in sb::SmolyakBasis 
	@assert(sb.SpOut==false,"Derivatives Not Sparse")	 
	if ==(sb.NumDeriv,1) 
		inds,vals = findnz(sb.∂Ψ∂z) 
		dims = [sb.NumPts,sb.NumBasisFun,sb.D]
		maxind = inds[end]
		numel = prod(dims)
		return reshape([full(sb.∂Ψ∂z),zeros(numel - maxind)],dims...)
	end
	if ==(sb.NumDeriv,2) 
		inds,vals = findnz(sb.∂2Ψ∂z2)
		dims = [sb.NumPts,sb.NumBasisFun,sb.D,sb.D] 
		maxind = inds[end]
		numel = prod(dims)
		return reshape([full(sb.∂2Ψ∂z2),zeros(numel - maxind)],dims...)
	end
end	

function sparse2full(sb::SmolyakBasis,k::Int64)
	# k specifies which sparse derivative to return to full 
	@assert(sb.SpOut==false,"Derivatives Not Sparse")	 
	if ==(k,1) 
		inds,vals = findnz(sb.∂Ψ∂z) 
		dims = [sb.NumPts,sb.NumBasisFun,sb.D]
		maxind = inds[end]
		numel = prod(dims)
		return reshape([full(sb.∂Ψ∂z),zeros(numel - maxind)],dims...)
	end
	if ==(k,2) 
		inds,vals = findnz(sb.∂2Ψ∂z2)
		dims = [sb.NumPts,sb.NumBasisFun,sb.D,sb.D] 
		maxind = inds[end]
		numel = prod(dims)
		return reshape([full(sb.∂2Ψ∂z2),zeros(numel - maxind)],dims...)
	end
end	
