#= ******************* =#
#= Smolyak Basis type  =#
#= ******************* =#

"""
## Description

Smolyak Basis type. Both Anisotrophic and Isotrophic Grids are supported 
and they are constructed efficiently following the methodology outlined in
Judd, Maliar, Maliar, Valero (2014). The code is designed for Julia v0.5.

#### Fields

- `D :: Int64` : Dimensions
- `mu :: ScalarOrVec{Int64}` : Index of mu
- `lb :: Vector{Float64}`	: Lower Bounds of dimensions
- `ub :: Vector{Float64}`	: Upper Bounds of dimensions
- `Binds :: AA{Int64}` : Basis Function Indices for Smolyak Interpolant, f(D,mu)
- `NumPts :: Int64` : Number of points in = Num Rows BF
- `NumBF :: Int64` : Number of basis functions under D, mu = Num Cols BF
- `NumDeriv :: Int64`	: Number of derivatives ∈ {0,1,2}
- `NumDerivArgs :: Int64`	: 1st NumDerivArgs used (i.e. 1:NumDerivArgs)
- `max_order :: Int64` : Maximum order of polynomial for T
- `x :: Array` : Vector of coordinates at which SB is evaluated 
- `z :: Array`	: Transformed vector of coordinates into [-1,1]
- `T :: Matrix{Float64}` : 1-dim Chebyshev basis fns: level
- `dT :: Matrix{Float64}` : 1-dim Chebyshev basis fns: 1st derivative
- `d2T :: Matrix{Float64}` : 1-dim Chebyshev basis fns: 2nd derivative
- `BF :: AA{Float64}`	: Basis Funs
- `dBFdz :: AAA{Float64}` : 1st derivative basis funs wrt z
- `d2BFdz2 :: AAAA{Float64}` : 2nd derivative basis funs wrt z
- `dzdx :: Vector{Float64}` : Gradient of transform z→x
- `d2zdx2 :: Vector{Float64}`	: Diagonal of Hessian of transform z→x (0 in linear maps, so only useful in nonlinear mapping)
- `dBFdx :: AAA{Float64}` : 1st derivative basis funs wrt x
- `d2BFdx2 :: AAAA{Float64}` : 2nd derivative basis funs wrt x

**Notes**: 'AA{T}' is type alias for 'Array{Array{T,1},1}' and 'AAA{T}' is defined analogously, etc.
 See ?Smolyak for list of typealias used in the package.

## Constructor functions

The constructor function creates the fields to contain the Smolyak Basis. 
To construct the Basis Functions as described by `mu`, `lb`, `ub`, `NumDeriv`, 
`NumDerivArgs` then call makeBasis!(sb). 

#### Using Smolyak Grid

To construct the container for SmolyakBasis call:

`sb = SmolyakBasis(sg,NumDeriv,NumDerivArgs)`

where:

- `sg :: SmolyakGrid`
- `NumDeriv :: Int64=2`
- `NumDerivArgs :: Int64=sg.D`

Then need to make basis functions call: 

`makeBasis!(sb)`

#### Using Smolyak Grid at new point, x

To construct the container for SmolyakBasis on a new (set of) grid point(s), x, 
and a Smolyak Grid is defined call:

`sb = SmolyakBasis(x,sg,NumDeriv,NumDerivArgs)`

where:

- `x :: Array`
- `sg :: SmolyakGrid`
- `NumDeriv :: Int64=2`
- `NumDerivArgs :: Int64=sg.D`

#### Construct without Smolyak grid

To construct the container for SmolyakBasis on a new (set of) grid point(s), x, 
without specifying a Smolyak Grid call:

`sb = SmolyakBasis(x,mu,lb,ub,NumDeriv,NumDerivArgs,D)`
		
where:

- `x :: Array` 
- `mu :: ScalarOrVec{Int64}`
- `lb :: Vector{Float64} = -1*ones(Float64,D)`
- `ub :: Vector{Float64} = ones(Float64,D)`
- `NumDeriv :: Int64=2`
- `NumDerivArgs :: Int64=D`
- `D :: Int64 = length(mu)`

## Examples

If Smolyak Grid already defined and want to evaluate Smolyak Basis functions on the grid:
```julia
using Smolyak
mu = [2,2,2]
lb = -2.*ones(length(mu))
ub = 3.*ones(length(mu))
sg = SmolyakGrid(mu,lb,ub)
sb = SmolyakBasis(sg)
makeBasis!(sb)
```

If Smolyak Grid already defined and want to evaluate Smolyak Basis functions on new grid point(s), x:
```julia
using Smolyak
mu = [2,2,2]
lb = -2.*ones(length(mu))
ub = 3.*ones(length(mu))
sg = SmolyakGrid(mu,lb,ub)
NumObs = 1000
x = Vector{Float64}[Array{Float64}(length(mu)) for i = 1:NumObs]
for i in 1:NumObs
	x[i] = [ sg.lb[d]+( sg.ub[d]- sg.lb[d])*rand() for d in 1:sg.D]
end
sb = SmolyakBasis(x,sg)
makeBasis!(sb)
```

Without a Smolyak Grid:
```julia
using Smolyak
mu = [2,2,2]
lb = -2.*ones(length(mu))
ub = 3.*ones(length(mu))
NumObs = 1000
x = Vector{Float64}[Array{Float64}(length(mu)) for i = 1:NumObs]
for i in 1:NumObs
	x[i] = [ lb[d]+( ub[d]- lb[d])*rand() for d in 1:length(mu)]
end
sb = SmolyakBasis(x,mu,lb,ub)
makeBasis!(sb)
```
"""
mutable struct SmolyakBasis{T}
	D 			:: Int64					# Dimensions
	mu 			:: ScalarOrVec{Int64}		# Index of mu
	lb 			:: Vector{Float64}			# Lower Bounds of dimensions
	ub 			:: Vector{Float64}			# Upper Bounds of dimensions
	Binds 		:: AA{Int64} 				# Basis Function Indices for Smolyak Interpolant, f(D,mu)
	NumPts  	:: Int64					# Number of points in = Num Rows BF
	NumBF		:: Int64					# Number of basis functions under D, mu = Num Cols BF
	NumDeriv	:: Int64					# Number of derivatives: {0,1,2}
	NumDerivArgs:: Int64					# 1st NumDerivArgs used (i.e. 1:NumDerivArgs)
	max_order	:: Int64 					# Maximum order of polynomial for T
	x 			:: Vector{T}				# Vector of coordinates at which SB is evaluated 
	z 			:: Vector{T} 				# Transformed vector of coordinates into [-1,1]
	T 			:: Matrix{Float64} 			# 1-dim Chebyshev basis fn: level
	dT 			:: Matrix{Float64} 			# 1-dim Chebyshev basis fn: 1st derivative
	d2T 		:: Matrix{Float64} 			# 1-dim Chebyshev basis fn: 2nd derivative
	BF 			:: AA{Float64}				# Basis Funs
	dBFdz 		:: AAA{Float64} 			# 1st derivative basis funs wrt z
	d2BFdz2 	:: AAAA{Float64}			# 2nd derivative basis funs wrt z
	dzdx		:: Vector{Float64} 			# Gradient of transform z→x
	d2zdx2		:: Vector{Float64}			# Diagonal of Hessian of transform z→x (0 in linear maps, so only useful in nonlinear mapping)
	dBFdx 		:: AAA{Float64} 			# 1st derivative basis funs wrt x
	d2BFdx2 	:: AAAA{Float64}			# 2nd derivative basis funs wrt x

end

# Constructor function with conformable memory allocations. Need to makeBasis!(sb) to fill it in.
function SmolyakBasis(sg::SmolyakGrid; NumDeriv::Int64=2,NumDerivArgs::Int64=sg.D)
	
	# Components for evaluation of Basis Functions
	NumBF = length(sg.Binds)
	max_order = 0
	for i in eachindex(sg.Binds), j in 1:sg.D
		max_order = max(max_order,sg.Binds[i][j]+1)
	end

	# Preallocate memory for basis functions
	T = zeros(sg.D, max_order) 	
	dT = zeros(sg.D, max_order)
	d2T = zeros(sg.D, max_order)
	
	# For Basis Functions and transformation back: Allocate memory for BF, pinvBF, and derivatives -> then makeBasis!(sb)
	NumDerivArgs = max(1,NumDerivArgs)
	BF = Vector{Float64}[ones(NumBF) 
					for n in 1:sg.NumGrdPts]			# BF[n][p] where n =1:NumGrdPts, p=1:NumBF
	
	if NumDeriv==0
		dBFdz = AA{Float64}[]
		dzdx  = Float64[]
		d2BFdz2 = AAA{Float64}[]
		d2zdx2 = Float64[]
	elseif NumDeriv==1
		dBFdz = AA{Float64}[[ones(NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:sg.NumGrdPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:sg.D
			push!(dzdx,2/(sg.ub[i] - sg.lb[i]))
		end
		d2BFdz2 = AAA{Float64}[]
		d2zdx2 = Float64[]
	else
		NumDerivArgs>2 ? println("****** MAX OF 2 DERIVATIVES ******") : nothing
		dBFdz = AA{Float64}[[ones(NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:sg.NumGrdPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:sg.D
			push!(dzdx,2/(sg.ub[i] - sg.lb[i]))
		end
		d2BFdz2 = AAA{Float64}[[[ones(NumBF) 
					for i in k:NumDerivArgs]
					for k in 1:NumDerivArgs]
					for n in 1:sg.NumGrdPts] 			#= d2BFdz2[n][i][k][p] where n =1:NumGrdPts,
														 i is position of 1st derivative, 
														 k = j - i + 1 where j in position of 2nd derivative, and p=1:NumBF =#
		d2zdx2 = zeros(Float64,NumDerivArgs) 
	end

	dBFdx = deepcopy(dBFdz)
	d2BFdx2 = deepcopy(d2BFdz2)	

	return SmolyakBasis(sg.D, sg.mu,  sg.lb, sg.ub, sg.Binds, 
		sg.NumGrdPts, NumBF, NumDeriv, NumDerivArgs, max_order,
		sg.xGrid, sg.zGrid, T, dT, d2T, 
		BF, dBFdz, d2BFdz2, dzdx, d2zdx2, dBFdx, d2BFdx2)
end

function SmolyakBasis(x::Vector{Float64},sg::SmolyakGrid;NumDeriv::Int64=2,NumDerivArgs::Int64=sg.D)
		
		z = zeros(sg.D)
		NumPts = 1  
		x2z!(x,z,sg.lb,sg.ub) 						#= x should be D x NumPts =#

		# Components for evaluation of Basis Functions
		NumBF = length(sg.Binds)
		max_order = 0
		for i in eachindex(sg.Binds), j in 1:sg.D
			max_order = max(max_order,sg.Binds[i][j]+1)
		end

		# Preallocate memory for basis functions
		T = zeros(sg.D, max_order) 	
		dT = zeros(sg.D, max_order)
		d2T = zeros(sg.D, max_order)

		# For Basis Functions and transformation back: Allocate memory for BF, pinvBF, and derivatives -> then makeBasis!(sb)
		NumDerivArgs = max(1,NumDerivArgs)
		BF = Vector{Float64}[ones(NumBF) 
						for n in 1:NumPts]			# BF[n][p] where n =1:NumGrdPts, p=1:NumBF
		
		if NumDeriv==0
			dBFdz = AA{Float64}[]
			dzdx  = Float64[]
			d2BFdz2 = AAA{Float64}[]
			d2zdx2 = Float64[]
		elseif NumDeriv==1
			dBFdz = AA{Float64}[[ones(NumBF) 
						for i in 1:NumDerivArgs] 
						for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
															 i is position of 1st derivative,and p=1:NumBF =#
			dzdx  = Float64[]
			for i in 1:sg.D
				push!(dzdx,2/(sg.ub[i] - sg.lb[i]))
			end
			d2BFdz2 = AAA{Float64}[]
			d2zdx2 = Float64[]
		else
			NumDerivArgs>2 ? println("****** MAX OF 2 DERIVATIVES ******") : nothing
			dBFdz = AA{Float64}[[ones(NumBF) 
						for i in 1:NumDerivArgs] 
						for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
															 i is position of 1st derivative,and p=1:NumBF =#
			dzdx  = Float64[]
			for i in 1:sg.D
				push!(dzdx,2/(sg.ub[i] - sg.lb[i]))
			end
			d2BFdz2 = AAA{Float64}[[[ones(NumBF) 
						for i in k:NumDerivArgs]
						for k in 1:NumDerivArgs]
						for n in 1:NumPts] 			#= d2BFdz2[n][i][k][p] where n =1:NumGrdPts,
															 i is position of 1st derivative, 
															 k = j - i + 1 where j in position of 2nd derivative, and p=1:NumBF =#
			d2zdx2 = zeros(Float64,NumDerivArgs) 
		end

		dBFdx = deepcopy(dBFdz)
		d2BFdx2 = deepcopy(d2BFdz2)	

		return SmolyakBasis(sg.D, sg.mu,  sg.lb, sg.ub, sg.Binds, 
			NumPts, NumBF, NumDeriv, NumDerivArgs, max_order,
			x, z, T, dT, d2T, 
			BF, dBFdz, d2BFdz2, dzdx, d2zdx2, dBFdx, d2BFdx2)
end


function SmolyakBasis(x::AA{Float64},sg::SmolyakGrid;NumDeriv::Int64=2,NumDerivArgs::Int64=sg.D)
	
	z = Vector{Float64}[zeros(sg.D) for r in eachindex(x)]
	NumPts = length(x)  
	x2z!(x,z,sg.lb,sg.ub) 						#= x should be D x NumPts =#

	# Components for evaluation of Basis Functions
	NumBF = length(sg.Binds)
	max_order = 0
	for i in eachindex(sg.Binds), j in 1:sg.D
		max_order = max(max_order,sg.Binds[i][j]+1)
	end

	# Preallocate memory for basis functions
	T = zeros(sg.D, max_order) 	
	dT = zeros(sg.D, max_order)
	d2T = zeros(sg.D, max_order)

	# For Basis Functions and transformation back: Allocate memory for BF, pinvBF, and derivatives -> then makeBasis!(sb)
	NumDerivArgs = max(1,NumDerivArgs)
	BF = Vector{Float64}[ones(NumBF) 
					for n in 1:NumPts]			# BF[n][p] where n =1:NumGrdPts, p=1:NumBF
	
	if NumDeriv==0
		dBFdz = AA{Float64}[]
		dzdx  = Float64[]
		d2BFdz2 = AAA{Float64}[]
		d2zdx2 = Float64[]
	elseif NumDeriv==1
		dBFdz = AA{Float64}[[ones(NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:sg.D
			push!(dzdx,2/(sg.ub[i] - sg.lb[i]))
		end
		d2BFdz2 = AAA{Float64}[]
		d2zdx2 = Float64[]
	else
		NumDerivArgs>2 ? println("****** MAX OF 2 DERIVATIVES ******") : nothing
		dBFdz = AA{Float64}[[ones(NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:sg.D
			push!(dzdx,2/(sg.ub[i] - sg.lb[i]))
		end
		d2BFdz2 = AAA{Float64}[[[ones(NumBF) 
					for i in k:NumDerivArgs]
					for k in 1:NumDerivArgs]
					for n in 1:NumPts] 			#= d2BFdz2[n][i][k][p] where n =1:NumGrdPts,
														 i is position of 1st derivative, 
														 k = j - i + 1 where j in position of 2nd derivative, and p=1:NumBF =#
		d2zdx2 = zeros(Float64,NumDerivArgs) 
	end

	dBFdx = deepcopy(dBFdz)
	d2BFdx2 = deepcopy(d2BFdz2)	

	return SmolyakBasis(sg.D, sg.mu,  sg.lb, sg.ub, sg.Binds, 
		NumPts, NumBF, NumDeriv, NumDerivArgs, max_order,
		x, z, T, dT, d2T, 
		BF, dBFdz, d2BFdz2, dzdx, d2zdx2, dBFdx, d2BFdx2)
end 

function SmolyakBasis(x::Vector{Float64},shd::SmolyakHD;NumDeriv::Int64=2,NumDerivArgs::Int64=shd.D)
	
	z = zeros(shd.D)
	NumPts = 1  
	x2z!(x,z,shd.lb,shd.ub) 						#= x should be D x NumPts =#

	# Components for evaluation of Basis Functions
	NumBF = length(shd.Binds)
	max_order = 0
	for i in eachindex(shd.Binds), j in 1:shd.D
		max_order = max(max_order,shd.Binds[i][j]+1)
	end

	# Preallocate memory for basis functions
	T = zeros(shd.D, max_order) 	
	dT = zeros(shd.D, max_order)
	d2T = zeros(shd.D, max_order)

	# For Basis Functions and transformation back: Allocate memory for BF, pinvBF, and derivatives -> then makeBasis!(sb)
	NumDerivArgs = max(1,NumDerivArgs)
	BF = Vector{Float64}[ones(NumBF) 
					for n in 1:NumPts]			# BF[n][p] where n =1:NumGrdPts, p=1:NumBF
	
	if NumDeriv==0
		dBFdz = AA{Float64}[]
		dzdx  = Float64[]
		d2BFdz2 = AAA{Float64}[]
		d2zdx2 = Float64[]
	elseif NumDeriv==1
		dBFdz = AA{Float64}[[ones(NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:shd.D
			push!(dzdx,2/(shd.ub[i] - shd.lb[i]))
		end
		d2BFdz2 = AAA{Float64}[]
		d2zdx2 = Float64[]
	else
		NumDerivArgs>2 ? println("****** MAX OF 2 DERIVATIVES ******") : nothing
		dBFdz = AA{Float64}[[ones(NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:shd.D
			push!(dzdx,2/(shd.ub[i] - shd.lb[i]))
		end
		d2BFdz2 = AAA{Float64}[[[ones(NumBF) 
					for i in k:NumDerivArgs]
					for k in 1:NumDerivArgs]
					for n in 1:NumPts] 			#= d2BFdz2[n][i][k][p] where n =1:NumGrdPts,
														 i is position of 1st derivative, 
														 k = j - i + 1 where j in position of 2nd derivative, and p=1:NumBF =#
		d2zdx2 = zeros(Float64,NumDerivArgs) 
	end

	dBFdx = deepcopy(dBFdz)
	d2BFdx2 = deepcopy(d2BFdz2)	

	return SmolyakBasis(shd.D, shd.mu,  shd.lb, shd.ub, shd.Binds, 
		NumPts, NumBF, NumDeriv, NumDerivArgs, max_order,
		x, z, T, dT, d2T, 
		BF, dBFdz, d2BFdz2, dzdx, d2zdx2, dBFdx, d2BFdx2)
end 

function SmolyakBasis(x::AA{Float64},shd::SmolyakHD;NumDeriv::Int64=2,NumDerivArgs::Int64=shd.D)
	
	z = Vector{Float64}[zeros(sg.D) for r in eachindex(x)]
	NumPts = length(x)  
	x2z!(x,z,sg.lb,sg.ub) 						#= x should be D x NumPts =#

	x2z!(x,z,shd.lb,shd.ub) 						#= x should be D x NumPts =#

	# Components for evaluation of Basis Functions
	NumBF = length(shd.Binds)
	max_order = 0
	for i in eachindex(shd.Binds), j in 1:shd.D
		max_order = max(max_order,shd.Binds[i][j]+1)
	end

	# Preallocate memory for basis functions
	T = zeros(sg.D, max_order) 	
	dT = zeros(sg.D, max_order)
	d2T = zeros(sg.D, max_order)

	# For Basis Functions and transformation back: Allocate memory for BF, pinvBF, and derivatives -> then makeBasis!(sb)
	NumDerivArgs = max(1,NumDerivArgs)
	BF = Vector{Float64}[ones(NumBF) 
					for n in 1:NumPts]			# BF[n][p] where n =1:NumGrdPts, p=1:NumBF
	
	if NumDeriv==0
		dBFdz = AA{Float64}[]
		dzdx  = Float64[]
		d2BFdz2 = AAA{Float64}[]
		d2zdx2 = Float64[]
	elseif NumDeriv==1
		dBFdz = AA{Float64}[[ones(NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:shd.D
			push!(dzdx,2/(shd.ub[i] - shd.lb[i]))
		end
		d2BFdz2 = AAA{Float64}[]
		d2zdx2 = Float64[]
	else
		NumDerivArgs>2 ? println("****** MAX OF 2 DERIVATIVES ******") : nothing
		dBFdz = AA{Float64}[[ones(NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:shd.D
			push!(dzdx,2/(shd.ub[i] - shd.lb[i]))
		end
		d2BFdz2 = AAA{Float64}[[[ones(NumBF) 
					for i in k:NumDerivArgs]
					for k in 1:NumDerivArgs]
					for n in 1:NumPts] 			#= d2BFdz2[n][i][k][p] where n =1:NumGrdPts,
														 i is position of 1st derivative, 
														 k = j - i + 1 where j in position of 2nd derivative, and p=1:NumBF =#
		d2zdx2 = zeros(Float64,NumDerivArgs) 
	end

	dBFdx = deepcopy(dBFdz)
	d2BFdx2 = deepcopy(d2BFdz2)	

	return SmolyakBasis(shd.D, shd.mu,  shd.lb, shd.ub, shd.Binds, 
		NumPts, NumBF, NumDeriv, NumDerivArgs, max_order,
		x, z, T, dT, d2T, 
		BF, dBFdz, d2BFdz2, dzdx, d2zdx2, dBFdx, d2BFdx2)
end 

# Constructor function without Smolyak Grid Call
function SmolyakBasis(x::Vector{Float64}, mu::ScalarOrVec{Int64},
						lb::Vector{Float64}=-1*ones(Float64,length(mu)), ub::Vector{Float64}=ones(Float64,length(mu));
						NumDeriv::Int64=2,NumDerivArgs::Int64=length(mu),D::Int64=length(mu))
	
	NumGrdPts, Ginds = Smolyak.SmolIdx(tuple(mu...))
	Binds = Vector{Int64}[zeros(Int64,D) for r in 1:NumGrdPts]
	Smolyak.makeBasisIdx!(Binds,Ginds,tuple(mu...)) # Basis Function Indices
	z = zeros(D) 
	x2z!(x,z,lb,ub) 						#= x should be D x NumPts =#
	NumPts = 1
	
	# Make Grid and Indices
	
	# Components for evaluation of Basis Functions
	NumBF = length(Binds)
	max_order = 0
	for i in eachindex(Binds), j in 1:D
		max_order = max(max_order,Binds[i][j]+1)
	end

	# Preallocate memory for basis functions
	T = zeros(D, max_order) 	
	dT = zeros(D, max_order)
	d2T = zeros(D, max_order)

	# For Basis Functions and transformation back: Allocate memory for BF, pinvBF, and derivatives -> then makeBasis!(sb)
	NumDerivArgs = max(1,NumDerivArgs)
	BF = Vector{Float64}[ones(NumBF) 
					for n in 1:NumPts]			# BF[n][p] where n =1:NumGrdPts, p=1:NumBF
	
	if NumDeriv==0
		dBFdz = AA{Float64}[]
		dzdx  = Float64[]
		d2BFdz2 = AAA{Float64}[]
		d2zdx2 = Float64[]
	elseif NumDeriv==1
		dBFdz = AA{Float64}[[ones(NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:D
			push!(dzdx,2/(ub[i] - lb[i]))
		end
		d2BFdz2 = AAA{Float64}[]
		d2zdx2 = Float64[]
	else
		NumDerivArgs>2 ? println("****** MAX OF 2 DERIVATIVES ******") : nothing
		dBFdz = AA{Float64}[[ones(NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:D
			push!(dzdx,2/(ub[i] - lb[i]))
		end
		d2BFdz2 = AAA{Float64}[[[ones(NumBF) 
					for i in k:NumDerivArgs]
					for k in 1:NumDerivArgs]
					for n in 1:NumPts] 			#= d2BFdz2[n][i][k][p] where n =1:NumGrdPts,
														 i is position of 1st derivative, 
														 k = j - i + 1 where j in position of 2nd derivative, and p=1:NumBF =#
		d2zdx2 = zeros(Float64,NumDerivArgs) 
	end

	dBFdx = deepcopy(dBFdz)
	d2BFdx2 = deepcopy(d2BFdz2)	

	return SmolyakBasis(D, mu,  lb, ub, Binds, 
		NumPts, NumBF, NumDeriv, NumDerivArgs, max_order,
		x, z, T, dT, d2T, 
		BF, dBFdz, d2BFdz2, dzdx, d2zdx2, dBFdx, d2BFdx2)
end 

# Constructor function without Smolyak Grid Call
function SmolyakBasis(x::AA{Float64}, mu::ScalarOrVec{Int64},
						lb::Vector{Float64}=-1*ones(Float64,length(mu)), ub::Vector{Float64}=ones(Float64,length(mu));
						NumDeriv::Int64=2,NumDerivArgs::Int64=length(mu),D::Int64=length(mu))
	
	NumGrdPts, Ginds = Smolyak.SmolIdx(tuple(mu...))
	Binds = Vector{Int64}[zeros(Int64,D) for r in 1:NumGrdPts]
	Smolyak.makeBasisIdx!(Binds,Ginds,tuple(mu...)) # Basis Function Indices
	z = Vector{Float64}[zeros(D) for r in eachindex(x)]
	x2z!(x,z,lb,ub) 						#= x should be D x NumPts =#
	NumPts = length(x) 
	
	# Make Grid and Indices
	
	# Components for evaluation of Basis Functions
	NumBF = length(Binds)
	max_order = 0
	for i in eachindex(Binds), j in 1:D
		max_order = max(max_order,Binds[i][j]+1)
	end

	# Preallocate memory for basis functions
	T = zeros(D, max_order) 	
	dT = zeros(D, max_order)
	d2T = zeros(D, max_order)

	# For Basis Functions and transformation back: Allocate memory for BF, pinvBF, and derivatives -> then makeBasis!(sb)
	NumDerivArgs = max(1,NumDerivArgs)
	BF = Vector{Float64}[ones(NumBF) 
					for n in 1:NumPts]			# BF[n][p] where n =1:NumGrdPts, p=1:NumBF
	
	if NumDeriv==0
		dBFdz = AA{Float64}[]
		dzdx  = Float64[]
		d2BFdz2 = AAA{Float64}[]
		d2zdx2 = Float64[]
	elseif NumDeriv==1
		dBFdz = AA{Float64}[[ones(NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:D
			push!(dzdx,2/(ub[i] - lb[i]))
		end
		d2BFdz2 = AAA{Float64}[]
		d2zdx2 = Float64[]
	else
		NumDerivArgs>2 ? println("****** MAX OF 2 DERIVATIVES ******") : nothing
		dBFdz = AA{Float64}[[ones(NumBF) 
					for i in 1:NumDerivArgs] 
					for n in 1:NumPts]			#= dBFdz[n][i][p] where n = 1:NumGrdPts, 
														 i is position of 1st derivative,and p=1:NumBF =#
		dzdx  = Float64[]
		for i in 1:D
			push!(dzdx,2/(ub[i] - lb[i]))
		end
		d2BFdz2 = AAA{Float64}[[[ones(NumBF) 
					for i in k:NumDerivArgs]
					for k in 1:NumDerivArgs]
					for n in 1:NumPts] 			#= d2BFdz2[n][i][k][p] where n =1:NumGrdPts,
														 i is position of 1st derivative, 
														 k = j - i + 1 where j in position of 2nd derivative, and p=1:NumBF =#
		d2zdx2 = zeros(Float64,NumDerivArgs) 
	end

	dBFdx = deepcopy(dBFdz)
	d2BFdx2 = deepcopy(d2BFdz2)	

	return SmolyakBasis(D, mu,  lb, ub, Binds, 
		NumPts, NumBF, NumDeriv, NumDerivArgs, max_order,
		x, z, T, dT, d2T, 
		BF, dBFdz, d2BFdz2, dzdx, d2zdx2, dBFdx, d2BFdx2)
end 



# -------- New state to evaluate Smolyak Interpolant f(D,mu) where sb.use ------

# New state vector -> applicable when sb.x is a vector
function new_x!(sb::SmolyakBasis,x::Vector{Float64})
	for d in eachindex(x)
		sb.x[d] = x[d]
	end
	x2z!(sb.x,sb.z,sb.lb,sb.ub)
end

# ----------- Chebyshev Polynomials & derivatives ----------- #

# Basis Function evaluated at a vector z
function Tn!(sb::SmolyakBasis)
	@inbounds for n in 1:sb.max_order, d in 1:sb.D
		if ==(n,1)
			sb.T[d,n] = 1.0
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 0.0 : nothing
			===(sb.NumDeriv,2) ? sb.d2T[d,n] = 0.0 : nothing				
		elseif ==(n,2)
			sb.T[d,n] = sb.z[d]
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 1.0 : nothing
			===(sb.NumDeriv,2) ? sb.d2T[d,n] = 0.0 : nothing				
		else
			sb.T[d,n] = 2.*sb.z[d]*sb.T[d,n-1] - sb.T[d,n-2]
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 2.*sb.T[d,n-1] + 2.*sb.z[d]*sb.dT[d,n-1] - sb.dT[d,n-2] : nothing
			===(sb.NumDeriv,2) ? sb.d2T[d,n] = 4.*sb.dT[d,n-1] + 2.*sb.z[d]*sb.d2T[d,n-1] - sb.d2T[d,n-2] : nothing
		end
	end
end

# Basis Function evaluated at a vector z at a grid point i = 1:NumPts
function Tn!(sb::SmolyakBasis,i::Int64)
	@inbounds for n in 1:sb.max_order, d in 1:sb.D
		if ==(n,1)
			sb.T[d,n] = 1.0
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 0.0 : nothing
			===(sb.NumDeriv,2) ? sb.d2T[d,n] = 0.0 : nothing				
		elseif ==(n,2)
			sb.T[d,n] = sb.z[i][d]
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 1.0 : nothing
			===(sb.NumDeriv,2) ? sb.d2T[d,n] = 0.0 : nothing				
		else
			sb.T[d,n] = 2.*sb.z[i][d]*sb.T[d,n-1] - sb.T[d,n-2]
			>=(sb.NumDeriv,1) ? sb.dT[d,n] = 2.*sb.T[d,n-1] + 2.*sb.z[i][d]*sb.dT[d,n-1] - sb.dT[d,n-2] : nothing
			===(sb.NumDeriv,2) ? sb.d2T[d,n] = 4.*sb.dT[d,n-1] + 2.*sb.z[i][d]*sb.d2T[d,n-1] - sb.d2T[d,n-2] : nothing
		end
	end
end


# ----------- BF & derivaitves ----------- # 

# Construct Basis Function
function BF!(sb::SmolyakBasis, BFIdx::Int64, DimIdx::Int64, GridIdx::Int64=1)
	sb.BF[GridIdx][BFIdx] *= sb.T[DimIdx,sb.Binds[BFIdx][DimIdx]+1]   
end

# Construct dBFdz! constructs 1st derivative of BF wrt z ∈ [-1,1], the transformed domain of state vector.   for first n arguments of state vector 
function dBFdz!(sb::SmolyakBasis, BFIdx::Int64, DimIdx::Int64, GridIdx::Int64=1, N::Int64=sb.NumDerivArgs)
	@inbounds for d in 1:N 		# n specifies 1:n derivatives to avoid unnecessary computations is derivatives of first n arguments required
		if ==(d,DimIdx)
			sb.dBFdz[GridIdx][d][BFIdx] *= sb.dT[DimIdx,sb.Binds[BFIdx][DimIdx]+1 ]
		else
			sb.dBFdz[GridIdx][d][BFIdx] *= sb.T[DimIdx,sb.Binds[BFIdx][DimIdx]+1 ]
		end
	end
end

# Construct d2BFdz2! constructs 1st derivative of BF wrt z ∈ [-1,1], the transformed domain of state vector for first n arguments of state vector 
function d2BFdz2!(sb::SmolyakBasis, BFIdx::Int64, DimIdx::Int64, GridIdx::Int64=1, N::Int64=sb.NumDerivArgs)
	@inbounds for i in 1:N, j in i:N
		k = j-i+1 				# Translate position of argument to index used to access the jth argument's derivative.
		if DimIdx==i==j
			sb.d2BFdz2[GridIdx][i][k][BFIdx] *= sb.d2T[DimIdx,sb.Binds[BFIdx][DimIdx]+1]
		elseif DimIdx==j!=i  
			sb.d2BFdz2[GridIdx][i][k][BFIdx] *= sb.dT[DimIdx,sb.Binds[BFIdx][DimIdx]+1]
		elseif DimIdx==i!=j
			sb.d2BFdz2[GridIdx][i][k][BFIdx] *= sb.dT[DimIdx,sb.Binds[BFIdx][DimIdx]+1]
		else
			sb.d2BFdz2[GridIdx][i][k][BFIdx] *= sb.T[DimIdx,sb.Binds[BFIdx][DimIdx]+1]
		end
	end
end

# Derivative of Transformations:dzdx, d2zdx2

#= Derivative constant over grid points under linear transform =#
function dzdx!(sb::SmolyakBasis,N::Int64=sb.NumDerivArgs)
	@inbounds for d in 1:N
		sb.dzdx[d] = 2/(sb.ub[d] - sb.lb[d])  
	end
end

# 1st derivatives for first n arguments of state vector
function dBFdx!(sb::SmolyakBasis,N::Int64=sb.NumDerivArgs)
	@inbounds for n in 1:sb.NumPts, d in 1:N, p in 1:sb.NumBF
		sb.dBFdx[n][d][p] = sb.dBFdz[n][d][p]*sb.dzdx[d] # matrix x scalar
	end
end

# 2nd derivatives for first n arguments of state vector
function d2BFdx2!(sb::SmolyakBasis,N::Int64=sb.NumDerivArgs)
	@inbounds for n in 1:sb.NumPts, i in 1:N, j in i:N, p in 1:sb.NumBF
		k = j-i+1 
		sb.d2BFdx2[n][i][k][p] = sb.d2BFdz2[n][i][k][p]*sb.dzdx[i]*sb.dzdx[j] # matrix x scalar x scalar
	end
end

#= --------------------------------------- =#
#= Construct Basis Functions & Derivatives =#
#= --------------------------------------- =#

# Initialise Basis Functions with 1's
function initBF!(sb::SmolyakBasis{Vector{Float64}},N::Int64=sb.NumDerivArgs)
	if ===(sb.NumDeriv,2)
		for n in eachindex(sb.x)
			fill!(sb.BF[n],1.)
			@inbounds for i in 1:N 
				fill!(sb.dBFdz[n][i],1.)
				for j in i:N
					k = j-i+1
					fill!(sb.d2BFdz2[n][i][k],1.)
				end
			end
		end
	elseif ===(sb.NumDeriv,1)
		for n in eachindex(sb.x)
			fill!(sb.BF[n],1.)
			@inbounds for i in 1:N 
				fill!(sb.dBFdz[n][i],1.)
			end
		end
	else
		for n in eachindex(sb.x)
			fill!(sb.BF[n],1.)
		end
	end
end

# Initialise Basis Functions with 1's
function initBF!(sb::SmolyakBasis{Float64},N::Int64=sb.NumDerivArgs)
	if ===(sb.NumDeriv,2)
		fill!(sb.BF[1],1.)
		@inbounds for i in 1:N 
			fill!(sb.dBFdz[1][i],1.)
			for j in i:N
				k = j-i+1
				fill!(sb.d2BFdz2[1][i][k],1.)
			end
		end
	elseif ===(sb.NumDeriv,1)
		fill!(sb.BF[1],1.)
		@inbounds for i in 1:N 
			fill!(sb.dBFdz[1][i],1.)
		end
	else
		fill!(sb.BF[1],1.)
	end
end

# Makes Basis Functions with sb.NumDeriv derivatives of the first n arguments of state vector
function makeBasis!(sb::SmolyakBasis{Vector{Float64}},N::Int64=sb.NumDerivArgs)
	sb.NumPts = length(sb.x)
	initBF!(sb,N) # Need to start with BF and derivatives as 1 because will take product over loops
	if ===(sb.NumDeriv,2)
		for i in eachindex(sb.x)
			Tn!(sb,i)
			for d in eachindex(sb.x[i]), p in eachindex(sb.BF[i])
				BF!(sb, p, d, i)
				dBFdz!(sb, p, d, i, N)
				d2BFdz2!(sb, p, d, i, N) 			# Hess
			end
		end
		dBFdx!(sb, N)
		d2BFdx2!(sb, N)
	elseif ===(sb.NumDeriv,1)
		for i in eachindex(sb.x)
			Tn!(sb,i)
			for d in eachindex(sb.x[i]), p in eachindex(sb.BF[i])
				BF!(sb, p, d, i)
				dBFdz!(sb, p, d, i, N) 				# Jac
			end
		end
		dBFdx!(sb, N)
	elseif ===(sb.NumDeriv,0)
		for i in eachindex(sb.x)
			Tn!(sb,i)
			for d in eachindex(sb.x[i]), p in eachindex(sb.BF[i])
				BF!(sb, p, d, i) 					# Now need to multiply over dimension to create a NBF-vector of Basis Function at grid point i. END OF LOOP → NumBF x NumPts for all grid points. 
			end 									
		end
	else
		print("Warning: sb.NumDeriv∈{0,1,2}")
	end
end

# Makes Basis Functions with sb.NumDeriv derivatives of the first n arguments of state vector
function makeBasis!(sb::SmolyakBasis{Float64},N::Int64=sb.NumDerivArgs)
	initBF!(sb,N) # Need to start with BF and derivatives as 1 because will take product over loops
	if ===(sb.NumDeriv,2)
		Tn!(sb)
		for d in eachindex(sb.x), p in eachindex(sb.BF[1])
			BF!(sb, p, d, 1)
			dBFdz!(sb, p, d, 1, N)
			d2BFdz2!(sb, p, d, 1, N) 			# Hess
		end
		dBFdx!(sb, N)
		d2BFdx2!(sb, N)
	elseif ===(sb.NumDeriv,1)
		Tn!(sb)
		for d in eachindex(sb.x), p in eachindex(sb.BF[1])
			BF!(sb, p, d, 1)
			dBFdz!(sb, p, d, 1, N) 				# Jac
		end
		dBFdx!(sb, N)
	elseif ===(sb.NumDeriv,0)
		Tn!(sb)
		for d in eachindex(sb.x), p in eachindex(sb.BF[1])
			BF!(sb, p, d, 1) 					# Now need to multiply over dimension to create a NBF-vector of Basis Function at grid point i. END OF LOOP → NumBF x NumPts for all grid points. 
		end 									
	else
		print("Warning: sb.NumDeriv∈{0,1,2}")
	end
end

# New state vector -> applicable when sb.x is a vector
function updateBasis!(sb::SmolyakBasis,x::Vector{Float64},N::Int64=sb.NumDerivArgs)
	new_x!(sb,x)
	makeBasis!(sb,N)
end

function show(io::IO, sb::SmolyakBasis)
	msg = "\n\tCreated Smolyak Basis:\n"
	msg *= "\t- Dim: $(sb.D), mu: $(sb.mu)\n"
	msg *= "\t- NumPts: $(sb.NumPts)\n"
	msg *= "\t- Number of Basis Functions: $(sb.NumBF)\n"
	if ==(sb.NumDeriv,0) msg *= "\t- No Derivative supplied. Do not call dBFdx or d2BFdx2.\n" end
	if ==(sb.NumDeriv,1) msg *= "\t- with dBFdx up to first $(sb.NumDerivArgs) arguments. Do not call d2BFdx2.\n" end
	if ==(sb.NumDeriv,2) msg *= "\t- with dBFdx & d2BFdx2 up to first $(sb.NumDerivArgs) arguments.\n" end
	print(io, msg)
end
