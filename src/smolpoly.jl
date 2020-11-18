#= ************************ =#
#= Smolyak Polynomial type  =#
#= ************************ =#


struct SmolyakPoly{T<:Real}
	sb 			:: SmolyakBasis{T}
	coef 	 	:: Vector{T}
	value 		:: Vector{T}
	gradient	:: Vector{T}
	hessian	 	:: Matrix{T}
	gradidx 	:: Dict{Int64, Vector{Int64}}
	hessidx 	:: Dict{Vector{Int64}, Vector{Int64}}
end

# --------------------- Constructor functions -------------------------- #

function SmolyakPoly(x::Vector{T},  ubf_type::Symbol, sk::SmolyakKernel{T}; 
							NumDeriv::Int64=0, NumDerivArgs::Int64=sk.D) where {T<:Real}
	sb = SmolyakBasis(x, ubf_type, sk; NumDeriv=NumDeriv, NumDerivArgs=NumDerivArgs)
	NumBF = length(sb.BF)
	if NumDeriv==0
		coef = zeros(T, NumBF)
		gradient = zeros(T, 0)
		hessian =  zeros(T, 0,0)
	elseif NumDeriv==1
		coef = zeros(T, NumBF)
		gradient = zeros(T, NumDerivArgs)
		hessian =  zeros(T, 0,0)
	elseif NumDeriv==2
		coef = zeros(T, NumBF)
		gradient = zeros(T, NumDerivArgs)
		hessian =  zeros(T, NumDerivArgs, NumDerivArgs)
	else
		throw("Max deriv 2")
	end

	gradidx, hessidx = getSparseIdx(sk; NumDeriv=NumDeriv)

	return SmolyakPoly(sb, coef, [zero(T)], gradient, hessian, gradidx, hessidx )
end

function SmolyakPoly( ubf_type::Symbol, sk::SmolyakKernel{T}; 
							NumDeriv::Int64=0, NumDerivArgs::Int64=sk.D) where {T<:Real}
	sb = SmolyakBasis(ubf_type, sk; NumDeriv=NumDeriv, NumDerivArgs=NumDerivArgs)
	NumBF = length(sb.BF)
	if NumDeriv==0
		coef = zeros(T, NumBF)
		gradient = zeros(T, 0)
		hessian =  zeros(T, 0,0)
	elseif NumDeriv==1
		coef = zeros(T, NumBF)
		gradient = zeros(T, NumDerivArgs)
		hessian =  zeros(T, 0,0)
	elseif NumDeriv==2
		coef = zeros(T, NumBF)
		gradient = zeros(T, NumDerivArgs)
		hessian =  zeros(T, NumDerivArgs, NumDerivArgs)
	else
		throw("Max deriv 2")
	end

	gradidx, hessidx = getSparseIdx(sk; NumDeriv=NumDeriv)

	return SmolyakPoly(sb, coef, [zero(T)], gradient, hessian, gradidx, hessidx )
end

function SmolyakPoly(sb::SmolyakBasis{T}; NumDeriv::Int64=0, NumDerivArgs::Int64=sb.sk.D) where {T<:Real}

	NumBF = length(sb.BF)
	if NumDeriv==0
		coef = zeros(T, NumBF)
		gradient = zeros(T, 0)
		hessian =  zeros(T, 0,0)
	elseif NumDeriv==1
		coef = zeros(T, NumBF)
		gradient = zeros(T, NumDerivArgs)
		hessian =  zeros(T, 0,0)
	elseif NumDeriv==2
		coef = zeros(T, NumBF)
		gradient = zeros(T, NumDerivArgs)
		hessian =  zeros(T, NumDerivArgs, NumDerivArgs)
	else
		throw("Max deriv 2")
	end

	gradidx, hessidx = getSparseIdx(sb.sk; NumDeriv=NumDeriv)

	return SmolyakPoly(sb, coef, [zero(T)], gradient, hessian, gradidx, hessidx )
end

function getSparseIdx(sk::SmolyakKernel{T}; NumDeriv::Int64=0) where {T<:Real}
	if NumDeriv==0
		gradidx = Dict{Int64,Vector{Int64}}()
		hessidx = Dict{Vector{Int64},Vector{Int64}}()
	elseif NumDeriv==1
		gradidx = Dict{Int64,Vector{Int64}}()
		for d in 1:sk.D
		    gradidx[d] = findall((VVtoMatrix(sk.BasisIdx)[:,d].>0))
		end
		hessidx = Dict{Vector{Int64},Vector{Int64}}()
	elseif NumDeriv==2
		gradidx =  Dict{Int64,Vector{Int64}}()
		for d in 1:sk.D
		    gradidx[d] = findall((VVtoMatrix(sk.BasisIdx)[:,d].>0))
		end
		hessidx = Dict{Vector{Int64},Vector{Int64}}()
		for d_i in 1:sk.D, d_j in d_i:sk.D
		    if d_i==d_j
		       hessidx[[d_i, d_i]] = findall(VVtoMatrix(sk.BasisIdx)[:,d_i].>1)
		    else 
		       hessidx[[d_i, d_j]] = intersect(gradidx[d_i],gradidx[d_j])
		    end
		end
	else 
		throw("Max number of derivatives is 2")
	end
	gradidx, hessidx
end

# --------------------- Updating functions -------------------------- #

# Add methods to manipulate basis functions with SmolyakPolynomial Types
state!(x::Vector{T}, sp::SmolyakPoly{T}) where T<:Real  = state!(x, sp.sb)

# Update Coefficient Vector
coef!(coef::Vector{T}, sp::SmolyakPoly{T}) where T<:Real = copyto!(sp.coef, coef)

# In-place update of value
function value!(sp::SmolyakPoly{T}) where T<:Real
	sp.value[1] = dot(sp.sb.BF, sp.coef)
end

# In-place update of derivative wrt d_i
function dWdx!(sp::SmolyakPoly{T}, d_i::Int64) where T<:Real
	sp.gradient[d_i] = zero(T)
	for n in sp.gradidx[d_i]
		sp.gradient[d_i]  += sp.sb.jacobian[n][d_i]*sp.coef[n]
	end
end

# In-place update of gradient
function gradient!(sp::SmolyakPoly{T}) where T<:Real
	fill!(sp.gradient, zero(T))
	@inbounds for d in 1:sp.sb.sk.D, n in sp.gradidx[d]
		sp.gradient[d] += sp.sb.jacobian[n][d]*sp.coef[n]
	end
end

# Get cross-derivatives wrt state variable d_i, d_j
function d2Wdx2!(sp::SmolyakPoly{T}, d_i::Int64, d_j::Int64) where T<:Real
	sp.hessian[d_i, d_j] = zero(T)
	for n in sp.hessidx[[d_i,d_j]]
		sp.hessian[d_i, d_j]  += sp.sb.hessian[n][d_i,d_j]*sp.coef[n]
	end
	sp.hessian[d_j, d_i] = sp.hessian[d_i, d_j]
end

# In place update of the hessian
function hessian!(sp::SmolyakPoly{T}) where T<:Real
	fill!(sp.hessian, zero(T))
	@inbounds for d_i in 1:sp.sb.sk.D, d_j in d_i:sp.sb.sk.D
		for n in sp.hessidx[[d_i,d_j]]
			sp.hessian[d_i,d_j] += sp.sb.hessian[n][d_i,d_j]*sp.coef[n]
		end 
		sp.hessian[d_j,d_i] = sp.hessian[d_i,d_j]
	end
end

# Update Smolyak Polynomial at a exising state vector
function SmolyakPoly!( sp::SmolyakPoly{T}; NumDeriv::Int64=0) where T<:Real
	SmolyakBasis!(sp.sb; NumDeriv=NumDeriv)
	if NumDeriv==0
		value!(sp)
	elseif NumDeriv==1
		value!(sp)
		gradient!(sp)
	elseif NumDeriv==2
		value!(sp)
		gradient!(sp)
		hessian!(sp)
	else 
		throw("Maximum number of hand coded derivatives is 2. Consider ForwardDiff for more")
	end
end

# Update Smolyak Polynomial at a new state vector
function SmolyakPoly!( x::Vector{T}, sp::SmolyakPoly{T}; NumDeriv::Int64=0) where T<:Real
	SmolyakBasis!(x, sp.sb; NumDeriv=NumDeriv)
	if NumDeriv==0
		value!(sp)
	elseif NumDeriv==1
		value!(sp)
		gradient!(sp)
	elseif NumDeriv==2
		value!(sp)
		gradient!(sp)
		hessian!(sp)
	else 
		throw("Maximum number of hand coded derivatives is 2. Consider ForwardDiff for more")
	end
end

# Some other methods
SmolyakBasis!(x::Vector{T}, sp::SmolyakPoly{T}; NumDeriv::Int64=0) where T<:Real = SmolyakBasis!(x, sp.sb; NumDeriv=NumDeriv)
SmolyakBasis!(sp::SmolyakPoly{T}; NumDeriv::Int64=0) where T<:Real = SmolyakBasis!(sp.sb; NumDeriv=NumDeriv)

# -------------------------- Extraction functions ---------------------------- #

# State Vector
state( sp::SmolyakPoly{T}) where T<:Real  = deepcopy(sp.sb.state)

# State Vector
coef( x::Vector{T}, sp::SmolyakPoly{T}) where T<:Real  = deepcopy(sp.coef)

# Value of Smolyak Polynomial at state vector x
value(sp::SmolyakPoly{T}) where T<:Real = dot(sp.sb.BF, sp.coef)

# Value of Smolyak Polynomial at new state vector x
function value(x::Vector{T}, sp::SmolyakPoly{T}) where T<:Real
	SmolyakPoly!(x, sp; NumDeriv=0)
	return dot(sp.sb.BF, sp.coef)
end

# Value of Smolyak Polynomial at vector of new state vectors xx = [x_1; x_2, ... x_N]
value(xx::VV{T}, sp::SmolyakPoly{T}) where T<:Real = [value(x_n, sp) for x_n in xx]

# Derivative wrt state variable d_i
function dWdx(sp::SmolyakPoly{T}, d_i::Int64) where T<:Real
	dWdx = zero(T)
	for n in sp.gradidx[d_i]
		dWdx += sp.sb.jacobian[n][d_i]*sp.coef[n]
	end
	dWdx
end

# Derivative wrt state variable d_i new x
function dWdx(x::Vector{T}, sp::SmolyakPoly{T}, d_i::Int64) where T<:Real
	SmolyakBasis!(x, sp.sb; NumDeriv=1)
	return dWdx(sp, d_i)
end
dWdx(xx::VV{T}, sp::SmolyakPoly{T}, d_i::Int64) where T<:Real = [dWdx(x_n, sp, d_i) for x_n in xx]

# Get gradient
gradient(sp::SmolyakPoly{T}) where T<:Real = deepcopy(sp.gradient)

# Gradient at new x
function gradient(x::Vector{T}, sp::SmolyakPoly{T}) where T<:Real
	SmolyakPoly!(x, sp; NumDeriv=1)
	return gradient(sp)
end
gradient(xx::VV{T}, sp::SmolyakPoly{T}) where T<:Real = [gradient(x_n, sp) for x_n in xx]
	
# Cross-derivatives wrt state variable d_i, d_j
function d2Wdx2(sp::SmolyakPoly{T}, d_i::Int64, d_j::Int64) where T<:Real
	D2W = zero(T)
	for n in sp.hessidx[[d_i,d_j]]
		D2W += sp.sb.hessian[n][d_i,d_j]*sp.coef[n]
	end
	D2W
end

# Cross-derivatives wrt state variable d_i, d_j at a new x
function d2Wdx2(x::Vector{T}, sp::SmolyakPoly{T}, d_i::Int64, d_j::Int64) where T<:Real
	SmolyakBasis!(x, sp.sb; NumDeriv=2)
	return  d2Wdx2(sp, d_i, d_j)

end
d2Wdx2(xx::VV{T}, sp::SmolyakPoly{T}, d_i::Int64, d_j::Int64) where T<:Real = [d2Wdx2(x_n, sp, d_i, d_j) for x_n in xx]

# Hessian at new x
function hessian(x::Vector{T}, sp::SmolyakPoly{T}) where T<:Real
	SmolyakPoly!(x, sp; NumDeriv=2)
	return  deepcopy(sp.hessian)
end
hessian(xx::VV{T}, sp::SmolyakPoly{T}) where T<:Real = [hessian( x_n, sp) for x_n in xx]

# Extract information on basis function at x
function SPoutput(x::Vector{T}, sp::SmolyakPoly{T}; NumDeriv::Int64=0) where T<:Real
	SmolyakPoly!(x, sp; NumDeriv=NumDeriv)
	if NumDeriv==0
		return deepcopy(sp.value)
	elseif NumDeriv==1 
		return deepcopy(tuple(sp.value, sp.gradient))
	elseif NumDeriv==2 
		return deepcopy(tuple(sp.value, sp.gradient, sp.hessian))
	else 
		throw("Number of Derivative <= 2")
	end
end
SPoutput(xx::VV{T}, sp::SmolyakPoly{T}; NumDeriv::Int64=0) where T<:Real = [SPoutput(x_n, sp; NumDeriv=NumDeriv) for x_n in xx]
