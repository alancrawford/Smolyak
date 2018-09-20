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

# Vectors of points from arguments of polynomials / grid points to x-coodinates
z2x(z::Vector{T}, sk::SmolyakKernel{T}) where T<:Real = map(z_n->lineartransform.(z_n, sk.zbnds, sk.xbnds), z)
z2x(zz::VV{T}, sk::SmolyakKernel{T}) where T<:Real = map(z_n->lineartransform.(z_n, sk.zbnds, sk.xbnds), zz)
z2x(Z::Matrix{T}, sk::SmolyakKernel{T}) where T<:Real = z2x([z[n] for n in 1:size(Z,1)], sk)
z2x(sg::SmolyakGrid{T}) where T<:Real = z2x(sg.grid, sg.sk)


# x-coordinates to VV of z-coordinates
x2z(x::Vector{T}, sk::SmolyakKernel{T}) where T<:Real = map(x_n->lineartransform.(x_n, sk.xbnds, sk.zbnds), x)
x2z(xx::VV{T}, sk::SmolyakKernel{T}) where T<:Real = map(x_n->lineartransform.(x_n, sk.xbnds, sk.zbnds), xx)
x2z(X::Matrix{T}, sk::SmolyakKernel{T}) where T<:Real = x2z([x[n] for n in 1:size(X,1)], sk)
