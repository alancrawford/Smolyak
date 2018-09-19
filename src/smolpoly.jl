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

# Add methods to manipulate basis functions with SmolyakPolynomial Types
makeSmolyakBasis!(sp::SmolyakPoly{T}; NumDeriv::Int64=0) where T<:Real = makeSmolyakBasis!(sp.sb; NumDeriv=NumDeriv)
makeSmolyakBasis!(sp::SmolyakPoly{T}, x::Vector{T}; NumDeriv::Int64=0) where T<:Real = makeSmolyakBasis!(sp.sb, x; NumDeriv=NumDeriv)

# Update Coefficient Vector
update_coef!(sp::SmolyakPoly{T}, coef::Vector{T}) where T<:Real = copyto!(sp.coef, coef)

# Value of Smolyak Polynomial
getValue(sp::SmolyakPoly{T}) where T<:Real = dot(sp.sb.BF, sp.coef)

function makeValue!(sp::SmolyakPoly{T}) where T<:Real
	sp.value[1] = dot(sp.sb.BF, sp.coef)
end

# Make derivative wrt state variable d_i
function get_dWdx(sp::SmolyakPoly{T}, d_i::Int64) where T<:Real
	dWdx = zero(T)
	for n in sp.gradidx[d_i]
		dWdx += sp.sb.jacobian[n][d_i]*sp.coef[n]
	end
	dWdx
end

function make_dWdx!(sp::SmolyakPoly{T}, d_i::Int64) where T<:Real
	sp.gradient[d_i] = zero(T)
	for n in sp.gradidx[d_i]
		sp.gradient[d_i]  += sp.sb.jacobian[n][d_i]*sp.coef[n]
	end
end

function makeGradient!(sp::SmolyakPoly{T}) where T<:Real
	fill!(sp.gradient, zero(T))
	@inbounds for d in 1:sp.sb.sk.D, n in sp.gradidx[d]
		sp.gradient[d] += sp.sb.jacobian[n][d]*sp.coef[n]
	end
end

# make cross-derivatives wrt state variable d_i, d_j
function get_d2Wdx2(sp::SmolyakPoly{T}, d_i::Int64, d_j::Int64) where T<:Real
	D2W = zero(T)
	for n in sp.hessidx[[d_i,d_j]]
		D2W += sp.sb.hessian[n][d_i,d_j]*sp.coef[n]
	end
	D2W
end

# Get cross-derivatives wrt state variable d_i, d_j
function make_d2Wdx2!(sp::SmolyakPoly{T}, d_i::Int64, d_j::Int64) where T<:Real
	sp.hessian[d_i, d_j] = zero(T)
	for n in sp.hessidx[[d_i,d_j]]
		sp.hessian[d_i, d_j]  += sp.sb.hessian[n][d_i,d_j]*sp.coef[n]
	end
	sp.hessian[d_j, d_i] = sp.hessian[d_i, d_j]
end

# make cross-derivatives wrt state variable d_i, d_j
function makeHessian!(sp::SmolyakPoly{T}) where T<:Real
	fill!(sp.hessian, zero(T))
	@inbounds for d_i in 1:sp.sb.sk.D, d_j in d_i:sp.sb.sk.D
		for n in sp.hessidx[[d_i,d_j]]
			sp.hessian[d_i,d_j] += sp.sb.hessian[n][d_i,d_j]*sp.coef[n]
		end 
		sp.hessian[d_j,d_i] = sp.hessian[d_i,d_j]
	end
end

# Update Smolyak Polynomial : new coef new x
function makeSmolyakPoly!(θ::Vector{T}, sp::SmolyakPoly{T}, x::Vector{T}; NumDeriv::Int64=0) where T<:Real
	update_coef!(sp, θ)
	makeSmolyakBasis!(sp.sb, x; NumDeriv=NumDeriv)
	if NumDeriv==0
		makeValue!(sp)
	elseif NumDeriv==1
		makeValue!(sp)
		makeGradient!(sp)
	elseif NumDeriv==2
		makeValue!(sp)
		makeGradient!(sp)
		makeHessian!(sp)
	else 
		throw("Maximum number of hand coded derivatives is 2. Consider ForwardDiff for more")
	end
end

# Update Smolyak Polynomial : new coef same x
function makeSmolyakPoly!(θ::Vector{T}, sp::SmolyakPoly{T}; NumDeriv::Int64=0) where T<:Real
	update_coef!(sp, θ)
	makeSmolyakBasis!(sp.sb; NumDeriv=NumDeriv)
	if NumDeriv==0
		makeValue!(sp)
	elseif NumDeriv==1
		makeValue!(sp)
		makeGradient!(sp)
	elseif NumDeriv==2
		makeValue!(sp)
		makeGradient!(sp)
		makeHessian!(sp)
	else 
		throw("Maximum number of hand coded derivatives is 2. Consider ForwardDiff for more")
	end
end


# Update Smolyak Polynomial : same coef, new x
function makeSmolyakPoly!(sp::SmolyakPoly{T}, x::Vector{T}; NumDeriv::Int64=0) where T<:Real
	makeSmolyakBasis!(sp.sb, x; NumDeriv=NumDeriv)
	if NumDeriv==0
		makeValue!(sp)
	elseif NumDeriv==1
		makeValue!(sp)
		makeGradient!(sp)
	elseif NumDeriv==2
		makeValue!(sp)
		makeGradient!(sp)
		makeHessian!(sp)
	else 
		throw("Maximum number of hand coded derivatives is 2. Consider ForwardDiff for more")
	end
end
