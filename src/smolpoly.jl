
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
	coef 	:: Vector{Float64}
	f  		:: ScalarOrVec{Float64}
	dfdx	:: Array{Vector{Float64},1}
	d2fdx2	:: Array{Vector{Float64},2}
	d2fdx2	:: Array{Vector{Float64},2}

	function SmolyakPoly(sb::SmolyakBasis, coef::Vector)		
		f = At_mul_B(sb.BF,coef) 
		dfdx = make_dfdx(sb,coef)
		d2fdx2 = make_d2fdx2(sb,coef)
		new(coef, f, dfdx, d2fdx2)
	end

end

# In place fn Update
function f!(sp::SmolyakPoly,sb::SmolyakBasis,coef::Vector{Float64}=sp.coef)
	At_mul_B!(sp.f,sb.BF,coef) 
end

# In place 1st Derivative Update
function dfdx!(sp::SmolyakPoly,sb::SmolyakBasis,coef::Vector{Float64}=sp.coef,n::Int64=sb.D)
	for d in 1:n
		At_mul_B!(sb.dfdx[d],sb.dBFdx[d],coef)
	end
end

# In place 2nd Derivative Update
function d2fdx2!(sp::SmolyakPoly,sb::SmolyakBasis,coef::Vector{Float64}=sp.coef,n::Int64=sb.D)
	for j in 1:n, i in 1:n
		At_mul_B!(sb.dfdx[i,j],sb.d2BFdx2[i,j],coef)
	end
end

#= New Coefficient =#
function coef!(sp::SmolyakPoly, sb::SmolyakBasis, f::Vector{Float64}=sp.f)
	At_mul_B!(sp.coef,sb.pinvBF,f) 
end

# Update Smolyak Poly f, dfdx, d2fdx2
function SmolyakPoly!(sp::SmolyakPoly,sb::SmolyakBasis,coef::Vector{Float64}=sp.coef,NumDeriv::Int64=0,n::Int64=sb.D)
	if is(NumDeriv,2)
		f!(sp,sb,coef) 
		dfdx!(sp,sb,coef,n)
		d2fdx2!(sp,sb,coef,n)
	elseif is(NumDeriv,1)
		f!(sp,sb,coef) 
		dfdx!(sp,sb,coef,n)
	else
		f!(sp,sb,coef) 
	end
end







