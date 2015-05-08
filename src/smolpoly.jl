
#= Smolyak Interpolation =#

#= To do: Write a µltiple dispatch to pick a subset of indices for derivatives? =#

#= write for full, them look at sparse implementation =#

function ∂f∂x_fun(sb::SmolyakBasis, coef::Vector{Float64})
	Out = Array(Float64,sb.NumPts,sb.D)
	if ==(sb.SpOut,true)
		∂Ψ∂z = sparse2full(sb,1)
		for d in 1:sb.D
			@inbounds Out[:,d] = ∂Ψ∂z[:,:,d]'*sb.∂z∂x[d]*coef
		end
	else
		for d in 1:sb.D
			@inbounds Out[:,d] = sb.∂Ψ∂z[:,:,d]'*sb.∂z∂x[d]*coef
		end
	end
	return Out
end

function ∂2f∂x2_fun(sb::SmolyakBasis, coef::Vector{Float64})
	Out = Array(Float64,sb.NumPts,sb.D,sb.D)
	if ==(sb.SpOut,true)
		∂2Ψ∂z2 = sparse2full(sb,2)
		for j in 1:sb.D, i in 1:sb.D
			@inbounds Out[:,i,j] = ∂2Ψ∂z2[:,:,i,j]'*sb.∂2z∂x2[i,j]*coef
		end		
	else
		for j in 1:sb.D, i in 1:sb.D 
			@inbounds Out[:,i,j] = sb.∂2Ψ∂z2[:,:,i,j]'*sb.∂2z∂x2[i,j]*coef
		end
	end
	return Out
end

type SmolyakPoly
	coef 	:: Vector{Float64}
	f  		:: Vector{Float64}
	∂f∂x	:: Array{Float64,2}
	∂2f∂x2	:: Array{Float64,3}

	#= coefmap :: Vector # This records to link to larger coefficient vector =#

	function SmolyakPoly(sb::SmolyakBasis, coef::Vector)
		if ==(sb.NumDeriv, 0)
			f = sb.Ψ'*coef 
			∂f∂x = Array(Float64,1,1)
			∂2f∂x2 = Array(Float64,1,1,1)
		elseif ==(sb.NumDeriv,1)
			f = sb.Ψ'*coef 
			∂f∂x = ∂f∂x_fun(sb, coef)
			∂2f∂x2 = Array(Float64,1,1,1)
		else
			f = sb.Ψ'*coef 
			∂f∂x = ∂f∂x_fun(sb, coef)
			∂2f∂x2 = ∂2f∂x2_fun(sb, coef)
		end
		new(coef, f, ∂f∂x, ∂2f∂x2)
	end

	#= Constructor allowing user to specify how many derivatives in spoly k = {0,1,2} =#
	function SmolyakPoly(sb::SmolyakBasis, coef::Vector, k::Int64)
		@assert(0<=k<=2,"Number of derivatives of f µst be 0, 1, or 2")
		if ==(k, 0)
			f = sb.Ψ'*coef 
			∂f∂x = Array(Float64,1,1)
			∂2f∂x2 = Array(Float64,1,1,1)
		elseif ==(k,1)
			f = sb.Ψ'*coef 
			∂f∂x = ∂f∂x_fun(sb, coef)
			∂2f∂x2 = Array(Float64,1,1,1)
		else
			f = sb.Ψ'*coef 
			∂f∂x = ∂f∂x_fun(sb, coef)
			∂2f∂x2 = ∂2f∂x2_fun(sb, coef)
		end
		new(coef, f, ∂f∂x, ∂2f∂x2)
	end

end

#= Full Update =#
function new_sp_coef!(sp::SmolyakPoly, sb::SmolyakBasis, f::Vector{Float64})
	sp.coef = sb.pinvΨ'*f 
end

#= Partial Updating : wgt is the amount of weight to put on the new estimate =#
function wgt_new_sp_coef!(sp::SmolyakPoly, sb::SmolyakBasis, f::Vector{Float64}, wgt::Float64)
	sp.coef = wgt*sb.pinvΨ'*f  + (1.0-wgt)*sp.coef
end

# Calculate 1st derivative in place
function ∂f∂x_fun!(sp::SmolyakPoly, sb::SmolyakBasis)
	if ==(sb.SpOut,true)
		∂Ψ∂z = sparse2full(sb,1)
		for d in 1:sb.D
			@inbounds sp.∂f∂x[:,d] = sb.∂Ψ∂z[:,:,d]'*sb.∂z∂x[d]*sp.coef
		end
	else
		for d in 1:sb.D
			@inbounds sp.∂f∂x[:,d] = sb.∂Ψ∂z[:,:,d]'*sb.∂z∂x[d]*sp.coef
		end
	end
end
# Calculate 2nd derivative in place
function ∂2f∂x2_fun!(sp::SmolyakPoly, sb::SmolyakBasis)
	if ==(sb.SpOut,true)
		∂2Ψ∂z2 = sparse2full(sb,2)
		for j in 1:sb.D, i in 1:sb.D
			@inbounds sp.∂2f∂x2[:,i,j] = sb.∂2Ψ∂z2[:,:,i,j]'*sb.∂2z∂x2[i,j]*sp.coef
		end		
	else
		for j in 1:sb.D, i in 1:sb.D 
			@inbounds sp.∂2f∂x2[:,i,j] = sb.∂2Ψ∂z2[:,:,i,j]'*sb.∂2z∂x2[i,j]*sp.coef
		end
	end
end

#= Update f & derivatives as specified in sb  =#
function new_sp_f!(sp::SmolyakPoly, sb::SmolyakBasis)
	if ==(sb.NumDeriv, 0)
		sp.f = sb.Ψ'*sp.coef 
	elseif ==(sb.NumDeriv, 1)
		sp.f = sb.Ψ'*sp.coef 
		∂f∂x_fun!(sp, sb)
	else
		sp.f = sb.Ψ'*sp.coef 
		∂f∂x_fun!(sp, sb)
		∂2f∂x2_fun!(sp, sb)
	end
end

#= Can also specify the number of derivatives to update =#
function new_sp_f!(sp::SmolyakPoly, sb::SmolyakBasis, k::Int64)
	@assert(0<=k<=2,"Number of derivatives of f µst be 0, 1, or 2")
	if ==(k,2) & >=(sb.NumDeriv, k)
		sp.f = sb.Ψ'*sp.coef 
		∂f∂x_fun!(sp, sb)
		∂2f∂x2_fun!(sp, sb)
	elseif ==(k,1) & >=(sb.NumDeriv, k)
		sp.f = sb.Ψ'*sp.coef 
		∂f∂x_fun!(sp, sb)
		>(sb.NumDeriv, k) ? print("Warning: You have chosen not to update $(sp).∂2f∂x2") : nothing
	else
		sp.f = sb.Ψ'*sp.coef
		==(sb.NumDeriv-k,2) ? print("Warning: You have chosen not to update $(sp).∂f∂x or $(sp).∂2f∂x2") : nothing
		==(sp.NumDeriv-k,1) ? print("Warning: You have chosen not to update $(sp).∂2f∂x2") : nothing
	end
end



