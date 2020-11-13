# *************** SPREAD POLYNOMIALS *************** #

struct SpreadPoly{T<:Real} <: UnivariateBasisFunction{T}
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
SpreadPoly(x::Vector{T}, N::Int64, xbnds::VV{T}, D::Int64=length(x)) where T<:Real = 
		SpreadPoly(x, 
					similar(x), 
					D,
					N*ones(Int,D), 
					[ones(eltype(x), N+1) for d in 1:D],
					[ones(eltype(x), N+1) for d in 1:D], 
					[ones(eltype(x), N+1) for d in 1:D],
					xbnds,
					[[zero(eltype(x)),one(eltype(x))] for d in 1:D])

# Anisotrophic order, default bounds, original state space
SpreadPoly(x::Vector{T}, N::Vector{Int64}, xbnds::VV{T}, D::Int64=length(x)) where T<:Real = 
		SpreadPoly(x, 
					similar(x), 
					D,
					N, 
					[ones(eltype(x), N[d]+1) for d in 1:D],
					[ones(eltype(x), N[d]+1) for d in 1:D], 
					[ones(eltype(x), N[d]+1) for d in 1:D],
					xbnds,
					[[zero(eltype(x)),one(eltype(x))] for d in 1:D])

# Same Order for all dims, default bounds, original state space
SpreadPoly(x::Vector{T}, N::Int64, xbnds::VV{T}, zbnds::VV{T}, D::Int64=length(x)) where T<:Real = 
		SpreadPoly(x, 
					similar(x), 
					D,
					N*ones(Int,D), 
					[ones(eltype(x), N+1) for d in 1:D],
					[ones(eltype(x), N+1) for d in 1:D], 
					[ones(eltype(x), N+1) for d in 1:D],
					xbnds,
					zbnds)

# Anisotrophic default bounds, original state space
SpreadPoly(x::Vector{T}, N::Vector{Int64}, xbnds::VV{T}, zbnds::VV{T}, D::Int64=length(x)) where T<:Real = 
		SpreadPoly(x, 
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
Sn(s::T,Snm1::T,Snm2::T) where {T} = 2(1-2s)*Snm1 - Snm2 + 2s
dSn(s::T,Snm1::T,dSnm1::T,dSnm2::T) where {T} = 2(1-2s)*dSnm1 - 4Snm1 - dSnm2 + 2
d2Sn(s::T,dSnm1::T,d2Snm1::T,d2Snm2::T) where {T} = 2(1-2s)*d2Snm1 - 8dSnm1 - d2Snm2

# Use Recursive Definitiions
#= Note: I add a constant to S0 after other basis function formed using  
the spread polynomial definiton if `constant==true`. Benefit of this is 
tensor products in Smolyak polynomial do not all go zero as z→0, but those
that depend on z do go to zero. =#
function BasisFunctions!(spd::SpreadPoly; NumDeriv::Int=0, constant::Bool=true)
	copyto!(spd.z, x2z(spd))
	if NumDeriv==0
		for (d,z) in enumerate(spd.z)
			spd.BF[d][1] = zero(eltype(spd.z))
			spd.BF[d][2] = z
			@inbounds for n in 3:spd.N[d]+1
				spd.BF[d][n] = Sn(z, spd.BF[d][n-1], spd.BF[d][n-2])
			end
			constant ? spd.BF[d][1] = one(eltype(spd.z)) : nothing
		end
	elseif NumDeriv==1
		for (d,z) in enumerate(spd.z)
			spd.BF[d][1] = zero(eltype(spd.z))
			spd.BF[d][2] = z
			spd.dBF[d][1] = zero(eltype(spd.z))
			spd.dBF[d][2] = one(eltype(spd.z))
			@inbounds for n in 3:spd.N[d]+1
				spd.BF[d][n] = Sn(z, spd.BF[d][n-1], spd.BF[d][n-2])
				spd.dBF[d][n] = dSn(z, spd.BF[d][n-1],  spd.dBF[d][n-1], spd.dBF[d][n-2])
			end
			constant ? spd.BF[d][1] = one(eltype(spd.z)) : nothing
		end
	elseif NumDeriv==2
		for (d,z) in enumerate(spd.z)
			spd.BF[d][1] = zero(eltype(spd.z))
			spd.BF[d][2] = z
			spd.dBF[d][1] = zero(eltype(spd.z))
			spd.dBF[d][2] = one(eltype(spd.z))
			spd.d2BF[d][1] = zero(eltype(spd.z))
			spd.d2BF[d][2] = zero(eltype(spd.z))
			@inbounds for n in 3:spd.N[d]+1
				spd.BF[d][n] = Sn(z, spd.BF[d][n-1], spd.BF[d][n-2])
				spd.dBF[d][n] = dSn(z, spd.BF[d][n-1],  spd.dBF[d][n-1], spd.dBF[d][n-2])
				spd.d2BF[d][n] = d2Sn(z, spd.dBF[d][n-1],  spd.d2BF[d][n-1], spd.d2BF[d][n-2])
			end
			constant ? spd.BF[d][1] = one(eltype(spd.z)) : nothing
		end
	else
		throw("Max Number of Derivatives is 2. Consider ForwardDiff if want more.")	
	end
end
#=
# Use Recursive Definitiions - shifting order back so S[0] = S[-1].
function ShiftedBasisFunctions!(spd::SpreadPoly; NumDeriv::Int=0)
	copyto!(spd.z, x2z(spd))
	if NumDeriv==0
		for (d,z) in enumerate(spd.z)
			spd.BF[d][1] = z
			spd.BF[d][2] = Sn(z, spd.BF[d][1], zero(eltype(spd.z)))
			@inbounds for n in 3:spd.N[d]+1
				spd.BF[d][n] = Sn(z, spd.BF[d][n-1], spd.BF[d][n-2])
			end
		end
	elseif NumDeriv==1
		for (d,z) in enumerate(spd.z)
			spd.BF[d][1] = z
			spd.BF[d][2] = Sn(z, spd.BF[d][1], zero(eltype(spd.z)) )
			spd.dBF[d][1] = zero(eltype(spd.z))
			spd.dBF[d][2] = dSn(z, spd.BF[d][1],  spd.dBF[d][1], zero(eltype(spd.z)) )
			@inbounds for n in 3:spd.N[d]+1
				spd.BF[d][n] = Sn(z, spd.BF[d][n-1], spd.BF[d][n-2])
				spd.dBF[d][n] = dSn(z, spd.BF[d][n-1],  spd.dBF[d][n-1], spd.dBF[d][n-2])
			end
		end
	elseif NumDeriv==2
		for (d,z) in enumerate(spd.z)
			spd.BF[d][1] = z
			spd.BF[d][2] = Sn(z, spd.BF[d][1], zero(eltype(spd.z)) )
			spd.dBF[d][1] = zero(eltype(spd.z))
			spd.dBF[d][2] = dSn(z, spd.BF[d][1],  spd.dBF[d][1], zero(eltype(spd.z)) )
			spd.d2BF[d][1] = zero(eltype(spd.z))
			spd.d2BF[d][2] = d2Sn(z, spd.dBF[d][1],  spd.d2BF[d][1], zero(eltype(spd.z)) )
			@inbounds for n in 3:spd.N[d]+1
				spd.BF[d][n] = Sn(z, spd.BF[d][n-1], spd.BF[d][n-2])
				spd.dBF[d][n] = dSn(z, spd.BF[d][n-1],  spd.dBF[d][n-1], spd.dBF[d][n-2])
				spd.d2BF[d][n] = d2Sn(z, spd.dBF[d][n-1],  spd.d2BF[d][n-1], spd.d2BF[d][n-2])
			end
		end
	else
		throw("Max Number of Derivatives is 2. Consider ForwardDiff if want more.")	
	end
end
=#
