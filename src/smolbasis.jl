#= ******************* =#
#= Smolyak Basis type  =#
#= ******************* =#

"""
## Description

Smolyak Basis type. Both Anisotrophic and Isotrophic Grids are supported 
and they are constructed efficiently following the methodology outlined in
Judd, Maliar, Maliar, Valero (2014). The code is designed for Julia v0.7.

"""
struct SmolyakBasis{T<:Real}
	sk 			:: SmolyakKernel 			# Smolyak kernel
	ubf 		:: UnivariateBasisFunction  # Basis in each dimension for tensor prods for Smol Basis 			
	state 		:: Vector{T} 				# State Vector
	BF 			:: Vector{T}				# Basis Functions 
	jacobian 	:: VV{T} 					# Jacobian of basis funs wrt x
	hessian 	:: VM{T}					# Hessian of basis funs wrt x
end

# Create SmolyakBasis having defined bounds of x & z in Smolyak Kernel
function SmolyakBasis(x::Vector{T},  ubf_type::Symbol, sk::SmolyakKernel; 
							NumDeriv::Int64=2, NumDerivArgs::Int64=sk.D, DefaultDomain::Bool=true) where T<:Real

	@assert length(x) == sk.D

	# Note that defined in this way changes in sk.(x/z)bnds are linked directly to poly bounds fields 
	# If user doesn't want this then define new methods with deepcopy(sk.zbnds) etc. as arguments.
	if ubf_type == :ordinary
		DefaultDomain ? copyto!(sk.zbnds, sk.xbnds) : nothing  
		ubf = OrdinaryPoly(x, getOrder(sk), sk.xbnds, sk.zbnds) 
	elseif ubf_type == :chebyshev
		DefaultDomain ? copyto!(sk.zbnds, [[-one(T),one(T)] for d in 1:sk.D]) : nothing  
		ubf = ChebyshevPoly(x, getOrder(sk), sk.xbnds, sk.zbnds)
	elseif ubf_type == :spread
		DefaultDomain ? copyto!(sk.zbnds , [[zero(T),one(T)] for d in 1:sk.D]) : nothing  
		ubf = SpreadPoly(x, getOrder(sk), sk.xbnds, sk.zbnds)
	else 
		throw("Unsupported basis function, call ubf_type ∈ {:ordinary, :chebyshev, :spread}")
	end

	# Allocate mem for Basis Functionsm Derivatives
	NumBF = length(sk.BasisIdx)
	BF = ones(T, NumBF)

	if NumDeriv==0
		jacobian =  [Vector{T}(undef, 0) ]
		hessian = [Matrix{T}(undef, 0, 0)]
	elseif NumDeriv==1
		jacobian =  [zeros(T, NumDerivArgs) for n in 1:NumBF]
		hessian = [Matrix{T}(undef, 0, 0)]
	elseif NumDeriv==2
		jacobian =  [zeros(T, NumDerivArgs) for n in 1:NumBF]
		hessian =  [zeros(T, NumDerivArgs, NumDerivArgs) for n in 1:NumBF]
	else 
		throw("Up to 2 derivatives hand coded. If want more set NumDeriv=0 and use ForwardDiff.jl functionality")
	end

	SmolyakBasis(sk, ubf, x, BF, jacobian, hessian)
end

# Create SmolyakBasis having defined bounds of x & z in Smolyak Kernel
function SmolyakBasis( ubf_type::Symbol, sk::SmolyakKernel; 
					NumDeriv::Int64=2, NumDerivArgs::Int64=sk.D,DefaultDomain::Bool=true)

	x = rand(sk.D)
	T = eltype(x)

	@assert length(x) == sk.D

	if ubf_type == :ordinary
		DefaultDomain ? copyto!(sk.zbnds, sk.xbnds) : nothing  
		ubf = OrdinaryPoly(x, getOrder(sk), sk.xbnds, sk.zbnds)
	elseif ubf_type == :chebyshev
		DefaultDomain ? copyto!(sk.zbnds, [[-one(T),one(T)] for d in 1:sk.D]) : nothing  
		ubf = ChebyshevPoly(x, getOrder(sk), sk.xbnds, sk.zbnds)
	elseif ubf_type == :spread
		DefaultDomain ? copyto!(sk.zbnds , [[zero(T),one(T)] for d in 1:sk.D]) : nothing  
		ubf = SpreadPoly(x, getOrder(sk), sk.xbnds, sk.zbnds)
	else 
		throw("Unsupported basis function, call ubf_type ∈ {:ordinary, :chebyshev, :spread}")
	end

	# Allocate mem for Basis Functionsm Derivatives
	NumBF = length(sk.BasisIdx)
	BF = ones(T, NumBF)

	if NumDeriv==0
		jacobian =  [Vector{T}(undef, 0) ]
		hessian = [Matrix{T}(undef, 0, 0)]
	elseif NumDeriv==1
		jacobian =  [zeros(T, NumDerivArgs) for n in 1:NumBF]
		hessian = [Matrix{T}(undef, 0, 0)]
	elseif NumDeriv==2
		jacobian =  [zeros(T, NumDerivArgs) for n in 1:NumBF]
		hessian =  [zeros(T, NumDerivArgs, NumDerivArgs) for n in 1:NumBF]
	else 
		throw("Up to 2 derivatives hand coded. If want more set NumDeriv=0 and use ForwardDiff.jl functionality")
	end

	SmolyakBasis(sk, ubf, x, BF, jacobian, hessian)
end

# Get the order of the basis functions
function getOrder(BasisIdx::VV{T}, D::Int64) where T<:Real
	maxorder = zeros(Int64, D)
	for bfi in BasisIdx, d in 1:D
		maxorder[d] = max(maxorder[d] , bfi[d])
	end
	return maxorder
end

getOrder(sk::SmolyakKernel) = getOrder(sk.BasisIdx, sk.D)


# -------- New state to evaluate Smolyak Interpolant f(D,mu) where sb.use ------

# New state vector -> applicable when sb.x is a vector

new_state!(sb::SmolyakBasis,x::Vector{Float64}) = copyto!(sb.state, x)

# ----------- BF & derivaitves ----------- # 

# Construct Smolyak Basis Functions
function makeBF!(sb::SmolyakBasis)
	fill!(sb.BF, one(eltype(sb.ubf.x)))
	for (bf_num, bfi) in enumerate(sb.sk.BasisIdx), (d, nm1) in enumerate(bfi)
		sb.BF[bf_num] *= sb.ubf.BF[d][nm1+1]
	end
end

# Construct jacobian of Smolyak Basis Functions
function makeJacobian!(sb::SmolyakBasis)
	[fill!(sb.jacobian[n], one(eltype(sb.ubf.x))) for n in eachindex(sb.jacobian)]
	# di indexed variable that derivate is wrt , dk are the other variables
	for (bf_num, bfi) in enumerate(sb.sk.BasisIdx), (di, ni_m1) in enumerate(bfi), (dk, nk_m1) in enumerate(bfi)
		if di==dk
			sb.jacobian[bf_num][di] *= sb.ubf.dBF[di][ni_m1+1]*dzdx(sb.ubf, di)
		else
			sb.jacobian[bf_num][di] *= sb.ubf.BF[dk][nk_m1+1]
		end	
	end
end

# Make the Hessian
function makeHessian!(sb::SmolyakBasis)
	[fill!(sb.hessian[n], one(eltype(sb.ubf.x))) for n in eachindex(sb.BF)]
	# di and dj index variable that 2nd derivate is wrt , dk are the other variables
	for (bf_num, bfi) in enumerate(sb.sk.BasisIdx), (di, ni_m1) in enumerate(bfi), (dj, nj_m1) in enumerate(bfi), (dk, nk_m1) in enumerate(bfi)
		if dj>=di
			if dk==di==dj
				sb.hessian[bf_num][di, dj] *= sb.ubf.d2BF[di][ni_m1+1]*dzdx(sb.ubf, di)
			elseif dk==di!=dj
				sb.hessian[bf_num][di, dj] *= sb.ubf.dBF[di][ni_m1+1]*dzdx(sb.ubf, di)
			elseif di!=dj==dk
				sb.hessian[bf_num][di, dj] *= sb.ubf.dBF[dj][nj_m1+1]*dzdx(sb.ubf, dj)
			else
				sb.hessian[bf_num][di, dj] *= sb.ubf.BF[dk][nk_m1+1]
			end 
		else 
			sb.hessian[bf_num][di, dj] = sb.hessian[bf_num][dj, di]
		end
	end
end

# Smolyak Basis at existing x
function makeSmolyakBasis!(sb::SmolyakBasis; NumDeriv::Int64=0) where T<:Real
	if NumDeriv==0
		makeBF!(sb.ubf; NumDeriv=NumDeriv)
		makeBF!(sb)
	elseif NumDeriv==1
		makeBF!(sb.ubf; NumDeriv=NumDeriv)
		makeBF!(sb)
		makeJacobian!(sb)
	elseif NumDeriv==2
		makeBF!(sb.ubf; NumDeriv=NumDeriv)
		makeBF!(sb)
		makeJacobian!(sb)
		makeHessian!(sb)
	else
		throw("Max number hand coded deriv = 2. For more consider ForwardDiff")
	end
end

# Smolyak Basis at new x
function makeSmolyakBasis!(sb::SmolyakBasis, x::Vector{T}; NumDeriv::Int64=0) where T<:Real
	if NumDeriv==0
		new_state!(sb, x)
		makeBF!(sb.ubf; NumDeriv=NumDeriv)
		makeBF!(sb)
	elseif NumDeriv==1
		new_state!(sb, x)
		makeBF!(sb.ubf; NumDeriv=NumDeriv)
		makeBF!(sb)
		makeJacobian!(sb)
	elseif NumDeriv==2
		new_state!(sb, x)
		makeBF!(sb.ubf; NumDeriv=NumDeriv)
		makeBF!(sb)
		makeJacobian!(sb)
		makeHessian!(sb)
	else
		throw("Max number hand coded deriv = 2. For more consider ForwardDiff")
	end
end

#=
function show(io::IO, sb::SmolyakBasis)
	msg = "\n\tCreated Smolyak Basis:\n"
	msg *= "\t- Dim: $(sb.sk.D), mu: $(sb.sk.mu)\n"
	msg *= "\t- NumPts: $(sb.sk.N)\n"
	msg *= "\t- Number of Basis Functions: $(length(sb.BF))\n"
	if ==(sb.NumDeriv,0) msg *= "\t- No Derivative supplied. Do not call dBFdx or d2BFdx2.\n" end
	if ==(sb.NumDeriv,1) msg *= "\t- with dBFdx up to first $(sb.NumDerivArgs) arguments. Do not call d2BFdx2.\n" end
	if ==(sb.NumDeriv,2) msg *= "\t- with dBFdx & d2BFdx2 up to first $(sb.NumDerivArgs) arguments.\n" end
	print(io, msg)
end
=#
