# *************** ORDINARY POLYNOMIALS *************** #

struct OrdinaryPoly{T<:Real} <: UnivariateBasisFunction{T}
	x :: Vector{T} 	# State Vector
	z :: Vector{T} 	# State vector transformed into domain on basis functions
	D :: Int
	N :: Vector{Int} # Accuracy/Order of Basis Functions 
	BF :: VV{T} 	# Basis functions  
	dBF :: VV{T} 	# Matrix of basis functions diff once : D x N+1
	d2BF :: VV{T} 	# Matrix of basis functions diff twice : D x N+1
	xbnds :: VV{T}	# Bounds on state variable
	zbnds :: VV{T}	# Bounds on domain of basis function
end

# ************* #
# Constructors  #
# ************* #

# Isotrophic order, default bounds, original state space
OrdinaryPoly(x::Vector{T}, N::Int64, D::Int64=length(x)) where T<:Real = 
		OrdinaryPoly(x, 
					similar(x), 
					D,
					N*ones(Int,D), 
					[ones(eltype(x), N+1) for d in 1:D],
					[ones(eltype(x), N+1) for d in 1:D], 
					[ones(eltype(x), N+1) for d in 1:D],
					[[-Inf,Inf] for d in 1:D],
					[[-Inf,Inf] for d in 1:D])

# Anisotrophic order, default bounds, original state space
OrdinaryPoly(x::Vector{T}, N::Vector{Int64}, D::Int64=length(x)) where T<:Real = 
		OrdinaryPoly(x, 
					similar(x), 
					D,
					N, 
					[ones(eltype(x), N[d]+1) for d in 1:D],
					[ones(eltype(x), N[d]+1) for d in 1:D], 
					[ones(eltype(x), N[d]+1) for d in 1:D],
					[[-Inf,Inf] for d in 1:D],
					[[-Inf,Inf] for d in 1:D])

# Same Order for all dims, default bounds, original state space
OrdinaryPoly(x::Vector{T}, N::Int64, xbnds::VV{T}, zbnds::VV{T}, D::Int64=length(x)) where T<:Real = 
		OrdinaryPoly(x, 
					similar(x), 
					D,
					N*ones(Int,D), 
					[ones(eltype(x), N+1) for d in 1:D],
					[ones(eltype(x), N+1) for d in 1:D], 
					[ones(eltype(x), N+1) for d in 1:D],
					xbnds,
					zbnds)

# Anisotrophic default bounds, original state space
OrdinaryPoly(x::Vector{T}, N::Vector{Int64}, xbnds::VV{T}, zbnds::VV{T}, D::Int64=length(x)) where T<:Real = 
		OrdinaryPoly(x, 
					similar(x), 
					D,
					N, 
					[ones(eltype(x), N[d]+1) for d in 1:D],
					[ones(eltype(x), N[d]+1) for d in 1:D], 
					[ones(eltype(x), N[d]+1) for d in 1:D],
					xbnds,
					zbnds)

# ******************** #
# Make basis functions #
# ******************** #

# Use Recursive Definitiions - consistent method with other polynomial types i am defining
function makeBF!(op::OrdinaryPoly; NumDeriv::Int=0)
	if op.xbnds==[[-Inf,Inf] for d in 1:op.D]
		copyto!(op.z, op.x)
	else
		copyto!(op.z, x2z(op))
	end
	if NumDeriv==0
		for (d,z) in enumerate(op.z)
			op.BF[d][1] = one(eltype(op.z))
			for n in 2:op.N[d]
				op.BF[d][n] = z*op.BF[d][n-1]
			end
		end
	elseif NumDeriv==1
		for (d,z) in enumerate(op.z)
			op.BF[d][1] = one(eltype(op.z))
			op.dBF[d][1] = zero(eltype(op.z))
			for n in 2:op.N[d]
				op.BF[d][n] = z*op.BF[d][n-1]
				op.dBF[d][n] = op.BF[d][n-1] + z*op.dBF[d][n-1]
			end
		end
	elseif NumDeriv==2
		for (d,z) in enumerate(op.z)
			op.BF[d][1] = one(eltype(op.z))
			op.dBF[d][1] = zero(eltype(op.z))
			op.d2BF[d][1] = zero(eltype(op.z))
			for n in 2:op.N[d]
				op.BF[d][n] = z*op.BF[d][n-1]
				op.dBF[d][n] = op.BF[d][n-1] + z*op.dBF[d][n-1]
				op.d2BF[d][n] = (1+z)*op.dBF[d][n] + z*op.d2BF[d][n-1]
			end
		end
	else
		throw("Max Number of Derivatives is 2. Consider ForwardDiff if want more.")	
	end
end

