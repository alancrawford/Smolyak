
#= Smolyak Interpolation =#

#= To do: Write a multiple dispatch to pick a subset of indices for derivatives? =#

#= write for full, them look at sparse implementation =#

function ∂f∂x_fun(sb::SmolyakBasis, coef::Vector{Float64})
	Out = Array(Float64,sb.NumPts,sb.D)
	if ==(sb.SpOut,1)
		∂Psi∂z = sparse2full(sb,1)
		for d in 1:sb.D
			@inbounds Out[:,d] = ∂Psi∂z[:,:,d]*sb.∂z∂x[d]*coef
		end
	else
		for d in 1:sb.D
			@inbounds Out[:,d] = sb.∂Psi∂z[:,:,d]*sb.∂z∂x[d]*coef
		end
	end
	return Out
end

function ∂2f∂x2_fun(sb::SmolyakBasis, coef::Vector{Float64})
	Out = Array(Float64,sb.NumPts,sb.D,sb.D)
	if ==(sb.SpOut,1)
		∂2Psi∂z2 = sparse2full(sb,2)
		for i in 1:sb.D
			for j in 1:sb.D
				@inbounds Out[:,i,j] = ∂2Psi∂z2[:,:,i,j]*sb.∂2z∂x2[i,j]*coef
			end
		end		
	else
		for i in 1:sb.D
			for j in 1:sb.D
				@inbounds Out[:,i,j] = sb.∂2Psi∂z2[:,:,i,j]*sb.∂2z∂x2[i,j]*coef
			end
		end
	end
	return Out
end

type SmolyakPoly
	sb 		:: SmolyakBasis
	coef 	:: Vector{Float64}
	∂kf		:: Int64
	f  		:: Vector{Float64}
	∂f∂x	:: Array{Float64,2}
	∂2f∂x2	:: Array{Float64,3}

	#= coefmap :: Vector # This records to link to larger coefficient vector =#

	function SmolyakPoly(sb::SmolyakBasis, coef::Vector)
		if ==(sb.NumDeriv, 0)
			f = sb.Psi*coef 
			∂f∂x = Array(Float64,1,1)
			∂2f∂x2 = Array(Float64,1,1,1)
		elseif ==(sb.NumDeriv,1)
			f = sb.Psi*coef 
			∂f∂x = ∂f∂x_fun(sb, coef)
			∂2f∂x2 = Array(Float64,1,1,1)
		else
			f = sb.Psi*coef 
			∂f∂x = ∂f∂x_fun(sb, coef)
			∂2f∂x2 = ∂2f∂x2_fun(sb, coef)
		end
		new(sb, coef, sb.NumDeriv, f, ∂f∂x, ∂2f∂x2)
	end

	#= Constructor allowing user to specify how many derivatives in spoly k = {0,1,2} =#
	function SmolyakPoly(sb::SmolyakBasis, coef::Vector, k::Int64)
		@assert(0<=k<=2,"Number of derivatives of f must be 0, 1, or 2")
		if ==(k, 0)
			f = sb.Psi*coef 
			∂f∂x = Array(Float64,1,1)
			∂2f∂x2 = Array(Float64,1,1,1)
		elseif ==(k,1)
			f = sb.Psi*coef 
			∂f∂x = ∂f∂x_fun(sb, coef)
			∂2f∂x2 = Array(Float64,1,1,1)
		else
			f = sb.Psi*coef 
			∂f∂x = ∂f∂x_fun(sb, coef)
			∂2f∂x2 = ∂2f∂x2_fun(sb, coef)
		end
		new(sb, coef, k, f, ∂f∂x, ∂2f∂x2)
	end

end

#= Full Update =#
function update_spoly_coef!(spoly::SmolyakPoly, f::Vector{Float64})
	spoly.coef = spoly.sb.pinvPsi*f 
end

#= Partial Updating : wgt is the amount of weight to put on the new estimate =#
function update_spoly_coef!(spoly::SmolyakPoly, f::Vector{Float64}, wgt::Float64)
	spoly.coef = wgt*spoly.sb.pinvPsi*f  + (1.0-wgt)*spoly.coef
end

#= Update f & derivatives as specified in sp.sb (where sb::SmolyakBasis)
	 when coefficients have changed =#
function update_spoly_f!(sp::SmolyakPoly)
	if ==(sp.∂kf, 0)
		sp.f = sp.sb.Psi*sp.coef 
	elseif ==(sp.∂kf, 1)
		sp.f = sp.sb.Psi*sp.coef 
		sp.∂f∂x = sp.∂f∂x_fun(sp.sb, sp.coef)
	else
		sp.f = sp.sb.Psi*sp.coef 
		sp.∂f∂x = ∂f∂x_fun(sp.sb, sp.coef)
		sp.∂2f∂x2 = ∂2f∂x2_fun(sp.sb, sp.coef)
	end
end

#= Can also specify the number of derivatives to update =#
function update_spoly_f!(sp::SmolyakPoly, k::Int64)
	@assert(0<=k<=2,"Number of derivatives of f must be 0, 1, or 2")
	if ==(k,2) & >=(sp.∂kf, k)
		sp.f = sp.sb.Psi*sp.coef 
		sp.∂f∂x = sp.∂f∂x_fun(sp.sb, sp.coef)
		sp.∂2f∂x2 = ∂2f∂x2_fun(sp.sb, sp.coef)
	elseif ==(k,1) & >=(sp.∂fk, k)
		sp.f = sp.sb.Psi*sp.coef 
		sp.∂f∂x = sp.∂f∂x_fun(sp.sb, sp.coef)
		>(sp.∂kf, k) ? print("Warning: You have chosen not to update $(sp).∂2f∂x2") : nothing
	else
		sp.f = sp.sb.Psi*sp.coef
		==(sp.∂kf-k,2) ? print("Warning: You have chosen not to update $(sp).∂f∂x or $(sp).∂2f∂x2") : nothing
		==(sp.∂kf-k,1) ? print("Warning: You have chosen not to update $(sp).∂2f∂x2") : nothing
	end
end



