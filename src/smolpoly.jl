
#= Smolyak Interpolant =#

function make_dfdx(sb::SmolyakBasis, coef::Vector{Float64},n::Int64=sb.D)
	Out = [Array(Float64,sb.D) for i in 1:n]
	for d in 1:n
		Out[d] = At_mul_B(sb.dBFdx[d],coef)
	end
	return Out
end

function make_d2fdx2(sb::SmolyakBasis, coef::Vector{Float64},n::Int64=sb.D)
	Out = [Array(Float64,sb.D) for j in 1:n, i in 1:n]
	for j in 1:n, i in 1:n
		Out[i,j] = At_mul_B(sb.d2BFdx2[i,j],coef)
	end
	return Out
end

type SmolyakPoly
	coef 	 :: Vector{Float64}
	f  		 :: ScalarOrVec{Float64} 		# Value of Interpolant at each of the sb.NumPts
	dfdx	 :: Array{Vector{Float64},1} 	# Jacobian  where each element a vector of sb.NumPts x 1 wrt [i]
	d2fdx2	 :: Array{Vector{Float64},2} 	# Hessian where each element a vector of sb.NumPts x 1 for cross derivative [i,j] where i,j ∈{1,...,sb.D} 
	NumDeriv :: Int64 						# Number of Derivative
	function SmolyakPoly(sb::SmolyakBasis, coef::Vector, NumDeriv::Int64=sb.NumDeriv)		
		if is(NumDeriv,2)
			f = At_mul_B(sb.BF,coef) 
			dfdx = make_dfdx(sb,coef)
			d2fdx2 = make_d2fdx2(sb,coef)
		elseif is(NumDeriv,1)
			f = At_mul_B(sb.BF,coef) 
			dfdx = make_dfdx(sb,coef)
			d2fdx2 = [Array(Float64,sb.NumPts) for i in 1:sb.D, j in sb.D]
		elseif is(NumDeriv,0)
			f = At_mul_B(sb.BF,coef) 
			dfdx = [Array(Float64,sb.NumPts) for i in 1:sb.D]
			d2fdx2 = [Array(Float64,sb.NumPts) for i in 1:sb.D, j in sb.D]
		else
			println("Require: 0,1, or, 2 derivatives")
		end
		new(coef, f, dfdx, d2fdx2, NumDeriv)
	end

end

# In place fn Update
function f!(sp::SmolyakPoly,sb::SmolyakBasis,coef::Vector{Float64}=sp.coef)
	At_mul_B!(sp.f,sb.BF,coef) 
end

# In place 1st Derivative Update
function dfdx!(sp::SmolyakPoly,sb::SmolyakBasis,n::Int64=sb.D,coef::Vector{Float64}=sp.coef)
	for d in 1:n
		At_mul_B!(sp.dfdx[d],sb.dBFdx[d],coef)
	end
end

# In place 2nd Derivative Update
function d2fdx2!(sp::SmolyakPoly,sb::SmolyakBasis,n::Int64=sb.D,coef::Vector{Float64}=sp.coef)
	for j in 1:n, i in 1:n
		At_mul_B!(sp.d2fdx2[i,j],sb.d2BFdx2[i,j],coef)
	end
end

#= New Coefficient =#
function coef!(sp::SmolyakPoly, sb::SmolyakBasis, f::Vector{Float64}=sp.f)
	At_mul_B!(sp.coef,sb.pinvBF,f) 
end

# Update Smolyak Poly f, dfdx, d2fdx2
function SmolyakPoly!(sp::SmolyakPoly,sb::SmolyakBasis,n::Int64=sb.D,coef::Vector{Float64}=sp.coef,NumDeriv::Int64=0)
	if is(NumDeriv,2)
		f!(sp,sb,coef) 
		dfdx!(sp,sb,n,coef)
		d2fdx2!(sp,sb,n,coef)
	elseif is(NumDeriv,1)
		f!(sp,sb,coef) 
		dfdx!(sp,sb,n,coef)
	else
		f!(sp,sb,coef) 
	end
end

function show(io::IO, sp::SmolyakPoly)
	msg = "\n\tCreated Smolyak Interpolant:\n"
	if ==(sp.NumDeriv,0) msg *= "\t- No Derivative supplied. Do not call dfdx or d2fdx2.\n" end
	if ==(sp.NumDeriv,1) msg *= "\t- with dfdx. Do not call d2fdx2.\n" end
	if ==(sp.NumDeriv,2) msg *= "\t- with dfdx & d2dx2.\n" end
	print(io, msg)
end





