# Basis Function Types & aliases 

abstract type UnivariateBasisFunction{T<:Real} end

# ********* GENERIC FUNS ********** #

# Linear coordinate transforms
lineartransform(x::T,old_lb::T,old_ub::T,new_lb::T,new_ub::T) where {T<:Real} = 
	((x .- old_lb)./ (old_ub .- old_lb) ) .* (new_ub .- new_lb ) .+ new_lb 

lineartransform(x::T,old_bnds::Vector{T},new_bnds::Vector{T}) where {T<:Real} = 
	lineartransform(x,old_bnds[1],old_bnds[2],new_bnds[1],new_bnds[2]) 

dlt(old_lb::T,old_ub::T,new_lb::T,new_ub::T) where {T<:Real} = (new_ub .- new_lb )./ (old_ub .- old_lb)

dlt(old_bnds::Vector{T},new_bnds::Vector{T}) where {T<:Real} = dlt(old_bnds[1],old_bnds[2],new_bnds[1],new_bnds[2])

# Basis Functions

# linear transforms x to z (If did nonlinear transform would need d2xdz2 and so on.)
x2z(bf::T) where {T<:UnivariateBasisFunction} = lineartransform.(bf.x, bf.xbnds, bf.zbnds)
dxdz(bf::T) where {T<:UnivariateBasisFunction} = dlt.(bf.xbnds, bf.zbnds)

x2z(bf::T, d::Int64) where {T<:UnivariateBasisFunction} = lineartransform(bf.x[d], bf.xbnds[d], bf.zbnds[d])
dxdz(bf::T, d::Int64) where {T<:UnivariateBasisFunction} = dlt(bf.xbnds[d], bf.zbnds[d])

# z to x 
z2x(bf::T) where {T<:UnivariateBasisFunction} = lineartransform.(bf.z, bf.zbnds, bf.xbnds)
dzdx(bf::T) where {T<:UnivariateBasisFunction} = dlt.(bf.zbnds, bf.xbnds)

z2x(bf::T,  d::Int64) where {T<:UnivariateBasisFunction} = lineartransform(bf.z[d], bf.zbnds[d], bf.xbnds[d])
dzdx(bf::T, d::Int64) where {T<:UnivariateBasisFunction} = dlt(bf.zbnds[d], bf.xbnds[d])

#= ******************************************************** =#
#= Functions switching between z in [-1,1] and x in [lb,ub] =#
#= ******************************************************** =#

# Smolyak Grids
function x2z!(sg::SmolyakGrid) 
	for (n,x) in enumerate(sg.x)
		copyto!( 
			sg.z[n], lineartransform.(x, view(sg.xbnds,:,1), view(sg.xbnds,:,2), 
				view(sg.zbnds,:,1), view(sg.zbnds,:,2))
		)
	end
end

function z2x!(sg::SmolyakGrid) 
	for (n,z) in enumerate(sg.z)
		copyto!( 
			sg.x[n],lineartransform.(z, view(sg.zbnds,:,1), view(sg.zbnds,:,2), 
					view(sg.xbnds,:,1), view(sg.xbnds,:,2))
		)
	end
end

