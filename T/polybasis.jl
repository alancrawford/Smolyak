#= 
	-----------------------------------------------------------------
	Smolyak Basis Functions on a Smolyak Grid for Julia version 0.3.7  
	-----------------------------------------------------------------

This file contains code to define Smolyak Constructtion of Chebyshev
basis functions. It is compatible with both Anisotrophic and Isotrophic 
Grids and are constructed efficiently following the methodology outlined
 in JMMV (2014). The code is designed on latest stable Julia version: 0.3.7.

Key Refs: JMMV (2014), Burkhardt (2012), github: ECONFORGE/Smolyak

=#

#= ---------------- =#
#= Within loop funs =#
#= ---------------- =#

function Binc!(B::Array{Float64,2},Grid::Array{Float64,2},binds::Array{Int64,2})
	NBF,D = size(binds)
	for m in 1:NBF
		for n in 1:D
			@inbounds B[m,n] = abs(cos(binds[m,n]*Grid[1,n])) < 1e-14 ?
								0.0 : cos(binds[m,n]*Grid[1,n])
		end
	end
end

function dBinc!(dB::Array{Float64,2},Grid::Array{Float64,2},binds::Array{Int64,2})
	NBF,D = size(binds)
	for m in 1:NBF
		for n in 1:D
			@inbounds dB[m,n] = abs(-binds[m,n]*sin(binds[m,n]*Grid[1,n])) < 1e-14 ?
								0.0 : -binds[m,n]*sin(binds[m,n]*Grid[1,n])
		end
	end
end

function d2Binc!(d2B::Array{Float64,2},Grid::Array{Float64,2},binds::Array{Int64,2})
	NBF,D = size(binds)
	for m in 1:NBF
		for n in 1:D
			@inbounds d2B[m,n] = abs(-binds[m,n]*binds[m,n]*cos(binds[m,n]*Grid[1,n])) < 1e-14 ?
								0.0 : -binds[m,n]*binds[m,n]*cos(binds[m,n]*Grid[1,n])
		end
	end
end
 
function Phiinc!(Phi::Array{Float64,2},B::Array{Float64,2},binds::Array{Int64,2},i::Int64)
	NBF,D = size(binds)
	for n = 1:NBF
		for d in 1:D
			@inbounds Phi[i,n] *= B[n,d]
			if d==D && abs(Phi[i,n]) < 1e-14 
				Phi[i,n] = 0.0
			end
		end
	end
end 

function Phiinc!(Phi::Array{Float64,3},B::Array{Float64,2},binds::Array{Int64,2},i::Int64,j::Int64)
	NBF,D = size(binds)
	for n = 1:NBF
		for d in 1:D
			@inbounds Phi[i,n,j] *= B[n,d]
			if d==D && abs(Phi[i,n,j]) < 1e-14 
				Phi[i,n,j] = 0.0
			end
		end
	end
end 

function Phiinc!(Phi::Array{Float64,4},B::Array{Float64,2},binds::Array{Int64,2},i::Int64,j::Int64,k::Int64)
	NBF,D = size(binds)
	for n = 1:NBF
		for d in 1:D
			@inbounds Phi[i,n,j,k] *= B[n,d]
			if d==D && abs(Phi[i,n,j,k]) < 1e-14 
				Phi[i,n,j,k] = 0.0
			end
		end
	end
end 

function check_domain(X::Array{Float64,2}, sg::SmolyakGrid)
	N,D = size(X) 	# Size of Array of Grid points
	@assert(is(D,sg.D ),
		"\n\tError: Dimension mismatch between 
			Smolyak Grid and Matrix of Grid Points")	# Check correct number of dims
	# Check x is in bounds
	minX = Inf*ones(Float64,D)
	maxX = zeros(Float64,D)
	for i in 1:N
		for d in 1:D
			minX[d] = min(minX[d],X[i,d])
			maxX[d] = max(maxX[d],X[i,d])
		end
	end
	for d in 1:D 										# Report if X not in Bounds
		>(maxX[d],sg.ub[d]) ? 
			println("\tWarning: X[$d] > Upper Bound") : nothing
		<(minX[d],sg.lb[d]) ?
			println("\tWarning: X[$d] < Lower Bound") : nothing
	end
end

#= --------------------------------------- =#
#= Construct Basis Functions & Derivatives =#
#= --------------------------------------- =#

# Basis Functions for original points on Smolyak Grid
function Phifun(sg::SmolyakGrid, SpOut::Int64)
	Phi = ones(Float64,sg.NumGrdPts,sg.NumGrdPts)	# Allocate matrix of Basis functions of Polynomial
	B = Array(Float64,sg.NumGrdPts,sg.D )			# Allocate matrix of Chebyshev Basis Functions
	for i in 1:sg.NumGrdPts
		Binc!(B,sg.thetaGrid[i,:],sg.Binds)
		Phiinc!(Phi,B,sg.Binds,i)
	end
	return Phi
end	

# Basis Functions for D-dimensional matrix of N grid points in original domain x
function Phifun(X::Array{Float64,2}, sg::SmolyakGrid, DomCheck::Int64, SpOut::Int64)	
		
	# Convert to z in [-1,1]
	is(DomCheck,1) ? check_domain(X,sg) : nothing	
	thetaX = x2theta(X,sg.lb,sg.ub)					# Convert to z in [-1,1]

	N,D = size(X) 							# Size of Array of Grid points
	Phi = ones(Float64,N,sg.NumGrdPts)		# Allocate matrix of Basis functions of Polynomial
	B = Array(Float64,sg.NumGrdPts,sg.D )	# Allocate matrix of Chebyshev Basis Functions 
	for i in 1:N
		Binc!(B,thetaX[i,:],sg.Binds)
		Phiinc!(Phi,B,sg.Binds,i)
	end
	return Phi
end	

# 1st Derivatives for Collocationon z in [-1, 1]
function ∂Phifun(sg::SmolyakGrid, SpOut::Int64)
	∂Phi =  ones(Float64,sg.NumGrdPts,sg.NumGrdPts,sg.D )	# Allocate matrix of Basis functions of Polynomial
	B = Array(Float64,sg.NumGrdPts,sg.D );					# Allocate matrix of Chebyshev Basis Functions 
	dB = similar(B);										# Allocate matrix of 1st deerivatives Chebyshev Basis Functions 
	for i in 1:sg.NumGrdPts
		dBinc!(dB,sg.thetaGrid[i,:],sg.Binds)
		for j in 1:sg.D 
			Binc!(B,sg.thetaGrid[j,:],sg.Binds) 
			B[:,j] = dB[:,j]
			Phiinc!(∂Phi,B,sg.Binds,i,j)
		end
	end
	if ==(SpOut,1) 						# Return Sparse Matrix
		Zidx = sub2ind(size(∂Phi),findn(∂Phi)...) 
		return sparsevec(Zidx,∂Phi[Zidx])
	else											# Return Full Basis Matrix
		return ∂Phi
	end
end

# 1st Derivatives for D-dimensional matrix of N grid points in original domain x
function ∂Phifun(X::Array{Float64,2}, sg::SmolyakGrid, DomCheck::Int64, SpOut::Int64)

	# Convert to z in [-1,1]
	is(DomCheck,1) ? check_domain(X,sg) : nothing
	thetaX = x2theta(X,sg.lb,sg.ub)						# Convert to z in [-1,1]

	N,D = size(X) 								# Size of Array of Grid points
	∂Phi =  ones(Float64,N,sg.NumGrdPts,sg.D )	# Allocate Matrix of derivatives of Basis Funs
	B = Array(Float64,sg.NumGrdPts,sg.D )		# Allocate matrix of Chebyshev Basis Functions
	dB = similar(B)								# Allocate matrix of 1st derivatives of Chebyshev Basis Functions
	for i in 1:N
		dBinc!(dB,thetaX[i,:],sg.Binds)
		for j in 1:D
			Binc!(B,thetaX[i,:],sg.Binds) 
			B[:,j] = dB[:,j]
			Phiinc!(∂Phi,B,sg.Binds,i,j)
		end
	end
	if ==(SpOut,1) 					# Return Sparse Basis Matrix
		Zidx = sub2ind(size(∂Phi),findn(∂Phi)...) 
		return sparsevec(Zidx,∂Phi[Zidx])
	else										# Return Full Basis Matrix
		return ∂Phi
	end
end

# Second Derivatives for collocation 
function ∂2Phifun(sg::SmolyakGrid,SpOut::Int64)
	∂2Phi = ones(Float64,sg.NumGrdPts,sg.NumGrdPts,sg.D ,sg.D )	# Allocate Memory for 2nd Derivatives 
	B = Array(Float64,sg.NumGrdPts,sg.D )						# Allocate matrix of Chebyshev Basis Functions
	dB = similar(B)												# Allocate matrix of 1st deriv. Chebyshev Basis Functions
	d2B = similar(B)											# Allocate matrix of 2nd deriv. Chebyshev Basis Functions
	for i in 1:sg.NumGrdPts
		dBinc!(dB,sg.thetaGrid[i,:],sg.Binds)
		d2Binc!(dB,sg.thetaGrid[i,:],sg.Binds)
		for j in 1:sg.D 
			for k = 1:sg.D 
				Binc!(B,sg.thetaGrid[i,:],sg.Binds)
				B[:,j] = dB[:,j]
				j==k ? B[:,k] = d2B[:,k] : B[:,k] = dB[:,k] 				
				Phiinc!(∂2Phi,B,sg.Binds,i,j,k)
			end
		end
	end
	if ==(SpOut,1) 									# Return Sparse Basis Matrix
		Zidx = sub2ind(size(∂2Phi),findn(∂2Phi)...) 
		return sparsevec(Zidx,∂2Phi[Zidx])
	else														# Return Full Basis Matrix
		return ∂2Phi
	end
end

# Second derivatives for D-dimensional matrix of N grid points in original domain x
function ∂2Phifun(X::Array{Float64,2}, sg::SmolyakGrid, DomCheck::Int64, SpOut::Int64)

	is(DomCheck,1) ? check_domain(X,sg) : nothing		# Check the bounds of Grid Points in original domain
	thetaX = x2theta(X,sg.lb,sg.ub)							# Convert to theta in [0,π]

	N,D = size(X) 									# Size of Array of Grid points
	∂2Phi = ones(Float64,N,sg.NumGrdPts,sg.D ,sg.D )	# Allocate Memory for 2nd Derivatives 
	B = Array(Float64,sg.NumGrdPts,sg.D )			# Allocate matrix of Chebyshev Basis Functions
	dB = similar(B)									# Allocate matrix of 1st deriv. Chebyshev Basis Functions
	d2B = similar(B)								# Allocate matrix of 2nd deriv. Chebyshev Basis Functions
	for i in 1:N
		dBinc!(dB,thetaX[i,:],sg.Binds)
		d2Binc!(dB,thetaX[i,:],sg.Binds)
		for j in 1:D
			for k = 1:D
				Binc!(B,thetaX[i,:],sg.Binds)
				B[:,j] = dB[:,j]
				j==k ? B[:,k] = d2B[:,k] : B[:,k] = dB[:,k] 				
				Phiinc!(∂2Phi,B,sg.Binds,i,j,k)
			end
		end
	end
	if ==(SpOut,1) 								# Return Sparse Basis Matrix
		Zidx = sub2ind(size(∂2Phi),findn(∂2Phi)...) 
		return sparsevec(Zidx,∂2Phi[Zidx])
	else										# Return Full Basis Matrix
		return ∂2Phi
	end
end


#= ************** =#
#= PolyBasis type =#
#= ************** =#

IntOrVec = Union(Int64,Vector{Int64},Float64,Vector{Float64})
SpArray = Union(SparseMatrixCSC,Array)

type SmolyakBasis
	D 			:: Int64				# Dimensions
	mu 			:: IntOrVec				# Index of mu
	NumPts  	:: Int64				# Number of points in = Num Rows Phi
	NumBasisFun	:: Int64				# Number of basis functions under D, mu = Num Cols Phi
	Phi 		:: Array{Float64,2} 	# Basis Funs
	pinvPhi		:: Array{Float64,2} 	# Inverse Basis Funs
	∂Phi 		:: SpArray 				# 1st derivative basis funs
	∂2Phi 		:: SpArray 				# 2nd derivative basis funs
	SpOut		:: Int64				# Sparse Output indicator == 1 if Sparse, 0 Otherwise
	NumDeriv	:: Int64				# Number of derivatives: {0,1,2}

	function SmolyakBasis(sg::SmolyakGrid, NumDeriv::Int64, SpOut::Int64)
		if ==(NumDeriv,0)
			Phi = Phifun(sg, SpOut)
			pinvPhi = inv(Phi)
			∂Phi = Array(Float64,1,1,1)
			∂2Phi = Array(Float64,1,1,1,1)	
		end
		if ==(NumDeriv,1)
			Phi = Phifun(sg, SpOut )
			pinvPhi = inv(Phi)
			∂Phi = ∂Phifun(sg, SpOut )
			∂2Phi = Array(Float64,1,1,1,1)
		end
		if ==(NumDeriv,2)
			Phi = Phifun(sg, SpOut )
			pinvPhi = inv(Phi)
			∂Phi = ∂Phifun(sg, SpOut )
			∂2Phi = ∂2Phifun(sg, SpOut )
		end
		NGP, NBF = size(Phi)
		new(sg.D, sg.mu, NGP, NBF, Phi, pinvPhi, ∂Phi, ∂2Phi, SpOut, NumDeriv)
	end
	
	function SmolyakBasis(X::Array{Float64,2}, sg::SmolyakGrid, NumDeriv::Int64, DomCheck::Int64, SpOut::Int64 )
		if ==(NumDeriv,0)
			Phi = Phifun(X, sg, DomCheck, SpOut )
			pinvPhi = pinv(Phi)
			∂Phi = Array(Float64,1,1,1)
			∂2Phi = Array(Float64,1,1,1,1)	
		end
		if ==(NumDeriv,1)
			Phi = Phifun(X, sg, DomCheck, SpOut )
			pinvPhi = pinv(Phi)			
			∂Phi = ∂Phifun(X, sg, DomCheck, SpOut )
			∂2Phi = Array(Float64,1,1,1,1)	
		end
		if ==(NumDeriv,2)
			Phi = Phifun(X, sg, DomCheck, SpOut )
			pinvPhi = pinv(Phi)			
			∂Phi = ∂Phifun(X, sg, DomCheck, SpOut )
			∂2Phi = ∂2Phifun(X, sg, DomCheck, SpOut )
		end
		NGP, NBF = size(Phi)
		new(sg.D, sg.mu, NGP, NBF, Phi, pinvPhi, ∂Phi, ∂2Phi, SpOut, NumDeriv)	
	end
	
end

function show(io::IO, sb::SmolyakBasis)
	msg = "\n\tCreated Smolyak Basis:\n"
	msg *= "\t- Dim: $(sb.D), mu: $(sb.mu)\n"
	msg *= "\t- NumPts: $(sb.NumPts)\n"
	msg *= "\t- Number of Basis Functions: $(sb.NumBasisFun)\n"
	if ==(sb.NumDeriv,0) msg *= "\t- No Derivative supplied. Do not call ∂Phi or ∂2Phi!\n" end
	if ==(sb.NumDeriv,1) msg *= "\t- Sparse ∂Phi = $((sb.SpOut)==1)\n" end
	if ==(sb.NumDeriv,2) msg *= "\t- Sparse ∂Phi & ∂2Phi = $((sb.SpOut)==1)\n" end
	print(io, msg)
end

function sparse2full(sb::SmolyakBasis,k::Int64)
	#= k specifies which sparse derivative to return to full =#
	@assert(sb.SpOut>0,"Derivatives Not Sparse")	 
	if ==(k,1) 
		inds,vals = findnz(sb.∂Phi) 
		dims = [sb.NumPts,sb.NumBasisFun,sb.D]
		maxind = inds[end]
		numel = prod(dims)
		return reshape([full(sb.∂Phi),zeros(numel - maxind)],dims...)
	end
	if ==(k,2) 
		inds,vals = findnz(sb.∂2Phi)
		dims = [sb.NumPts,sb.NumBasisFun,sb.D,sb.D] 
		maxind = inds[end]
		numel = prod(dims)
		return reshape([full(sb.∂2Phi),zeros(numel - maxind)],dims...)
	end
end	
