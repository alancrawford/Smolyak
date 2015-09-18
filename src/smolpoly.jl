
#= Smolyak Interpolant =#

function dfdx!(sb::SmolyakBasis, coef::Vector{Float64},Out::Array{Float64,2}=Array(Float64,sb.NumPts,sb.D))
	for d in 1:sb.D
		Out[:,d] = sb.dBFdx[:,:,d]'*coef
	end
	return Out
end

function d2fdx2!(sb::SmolyakBasis, coef::Vector{Float64},Out::Array{Float64,3}=Array(Float64,sb.NumPts,sb.D,sb.D))
	for j in 1:sb.D, i in 1:sb.D 
		Out[:,i,j] = sb.d2BFdx2[:,:,i,j]'*coef
	end
	return Out
end

type SmolyakPoly
	coef 	:: Vector{Float64}
	f  		:: ScalarOrVec{Float64}
	dfdx	:: Array{Float64,2}
	d2fdx2	:: Array{Float64,3}

	#= coefmap :: Vector # This records to link to larger coefficient vector =#

	function SmolyakPoly(sb::SmolyakBasis, coef::Vector)
		
		f = sb.BF'*coef 
		dfdx = dfdx!(sb,coef)
		d2fdx2 = d2fdx2!(sb,coef)
		
		new(coef, f, dfdx, d2fdx2)
	end

end

function dfdx!(sb::SmolyakBasis,sp::SmolyakPoly)
	for d in 1:sb.D
		sp.dfdx[:,d] = sb.dBFdx[:,:,d]'*coef
	end
end


function d2fdx2!(sb::SmolyakBasis,sp::SmolyakPoly)
	for j in 1:sb.D, i in 1:sb.D 
		sp.d2fdx2[:,i,j] = sb.d2BFdx2[:,:,i,j]'*coef
	end
end

#= Full Update =#
function new_coef!(sp::SmolyakPoly, sb::SmolyakBasis, f::Vector{Float64})
	sp.coef = sb.pinvBF'*f 
end

#= Partial Updating : wgt is the amount of weight to put on the new estimate =#
function wgt_new_coef!(sp::SmolyakPoly, sb::SmolyakBasis, f::Vector{Float64}, wgt::Float64)
	sp.coef = wgt*sb.pinvBF'*f  + (1.0-wgt)*sp.coef
end




