# *************** ChEBYSHEV POLYNOMIALS *************** #

struct ChebyshevPoly{T<:Real} <: UnivariateBasisFunction{T}
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
ChebyshevPoly(x::Vector{T}, N::Int64, xbnds::VV{T}, D::Int64=length(x)) where T<:Real = 
		ChebyshevPoly(x, 
					similar(x), 
					D,
					N*ones(Int,D), 
					[ones(eltype(x), N+1) for d in 1:D],
					[ones(eltype(x), N+1) for d in 1:D], 
					[ones(eltype(x), N+1) for d in 1:D],
					xbnds,
					[[-one(eltype(x)),one(eltype(x))] for d in 1:D])

# Anisotrophic order, default bounds, original state space
ChebyshevPoly(x::Vector{T}, N::Vector{Int64}, xbnds::VV{T}, D::Int64=length(x)) where T<:Real = 
		ChebyshevPoly(x, 
					similar(x), 
					D,
					N, 
					[ones(eltype(x), N[d]+1) for d in 1:D],
					[ones(eltype(x), N[d]+1) for d in 1:D], 
					[ones(eltype(x), N[d]+1) for d in 1:D],
					xbnds,
					[[-one(eltype(x)),one(eltype(x))] for d in 1:D])

# Same Order for all dims, default bounds, original state space
ChebyshevPoly(x::Vector{T}, N::Int64, xbnds::VV{T}, zbnds::VV{T}, D::Int64=length(x)) where T<:Real = 
		ChebyshevPoly(x, 
					similar(x), 
					D,
					N*ones(Int,D), 
					[ones(eltype(x), N+1) for d in 1:D],
					[ones(eltype(x), N+1) for d in 1:D], 
					[ones(eltype(x), N+1) for d in 1:D],
					xbnds,
					zbnds)

# Anisotrophic default bounds, original state space
ChebyshevPoly(x::Vector{T}, N::Vector{Int64}, xbnds::VV{T}, zbnds::VV{T}, D::Int64=length(x)) where T<:Real = 
		ChebyshevPoly(x, 
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

# Recursive Definitions for cumulative definitions for n>1
Tn(z::T,Tnm1::T,Tnm2::T) where {T} = 2*z*Tnm1 - Tnm2
dTn(z::T,Tnm1::T,dTnm1::T,dTnm2::T) where {T} = 2*Tnm1 + 2*z*dTnm1 - dTnm2
d2Tn(z::T,dTnm1::T,d2Tnm1::T,d2Tnm2::T) where {T} = 4*dTnm1 + 2*z*d2Tnm1 - d2Tnm2

# Use Recursive Definitiions - consistent method with other polynomial types i am defining
function BasisFunctions!(cheb::ChebyshevPoly; NumDeriv::Int=0)
	copyto!(cheb.z, x2z(cheb))
	if NumDeriv==0
		for (d,z) in enumerate(cheb.z)
			cheb.BF[d][1] = one(eltype(cheb.z))
			cheb.BF[d][2] = z
			@inbounds for n in 3:cheb.N[d]+1
				cheb.BF[d][n] = Tn(z, cheb.BF[d][n-1], cheb.BF[d][n-2])
			end
		end
	elseif NumDeriv==1
		for (d,z) in enumerate(cheb.z)
			cheb.BF[d][1] = one(eltype(cheb.z))
			cheb.BF[d][2] = z
			cheb.dBF[d][1] = zero(eltype(cheb.z))
			cheb.dBF[d][2] = one(eltype(cheb.z))
			@inbounds for n in 3:cheb.N[d]+1
				cheb.BF[d][n] = Tn(z, cheb.BF[d][n-1], cheb.BF[d][n-2])
				cheb.dBF[d][n] = dTn(z, cheb.BF[d][n-1], cheb.dBF[d][n-1], cheb.dBF[d][n-2])
			end
		end
	elseif NumDeriv==2
		for (d,z) in enumerate(cheb.z)
			cheb.BF[d][1] = one(eltype(cheb.z))
			cheb.BF[d][2] = z
			cheb.dBF[d][1] = zero(eltype(cheb.z))
			cheb.dBF[d][2] = one(eltype(cheb.z))
			cheb.d2BF[d][1] = zero(eltype(cheb.z))
			cheb.d2BF[d][2] = zero(eltype(cheb.z))
			@inbounds for n in 3:cheb.N[d]+1
				cheb.BF[d][n] = Tn(z, cheb.BF[d][n-1], cheb.BF[d][n-2])
				cheb.dBF[d][n] = dTn(z, cheb.BF[d][n-1], cheb.dBF[d][n-1], cheb.dBF[d][n-2])
				cheb.d2BF[d][n] = d2Tn(z, cheb.dBF[d][n-1], cheb.d2BF[d][n-1], cheb.d2BF[d][n-2])
			end
		end
	else
		throw("Max Number of Derivatives is 2. Consider ForwardDiff if want more.")	
	end
end

