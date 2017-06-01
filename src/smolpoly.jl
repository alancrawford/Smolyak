#= ************************ =#
#= Smolyak Polynomial type  =#
#= ************************ =#

"""
## Description

Smolyak Polynomial type. Both Anisotrophic and Isotrophic Grids are supported 
and they are constructed efficiently following the methodology outlined in
Judd, Maliar, Maliar, Valero (2014). The code is designed for Julia v0.4.

#### Fields

- `NumPts :: Int64` : Number of points
- `NumCoef :: Int64` : Number of elements in Coefficient Vector
- `Coef :: Vector{Float64}`	: Coefficients of Smplyak Interpolant
- `Value :: ScalarOrVec{Float64}` : Value of Interpolant at each of the sb.NumPts
- `Grad :: AA{Float64}`	: Gradient of Interpolant up to first sp.NumDerivArgs  at grid point n = 1:NumPts
- `Hess :: AM{Float64}` : Hessian of Interpolant up to first sp.NumDerivArgs at grid point n = 1:NumPts
- `NumDeriv :: Int64` : Number of Derivative
- `NumDerivArgs :: Int64` : 1 to NumDerivArgs 
- `pinvBFt :: Matrix{Float64}` : Transpose of Moore-Penrose Pseudo Inverse of sb.BF (for case where NumPts ≥ NumCoef)

## Constructor functions

The constructor function creates the fields to contain the Smolyak Polynomial.

`sp = SmolyakPoly(sb, Coef, NumDeriv, NumDerivArgs, NumPts)`

where:

- `sb :: SmolyakBasis`
- `Coef :: Vector=rand(sb.NumBF)`
- `NumDeriv :: Int64=sb.NumDeriv`
- `NumDerivArgs :: Int64=sb.NumDerivArgs`
- `NumPts :: Int64=sb.NumPts`

After creating the fields for the Smolyak Polynomial, for a given coefficient vector, sp.Coed, 
fill-in the value, gradient and hessian of fields of the Smolyak Polynomial 

- For polynomial value(s): `makeValue!(sp)`
- For gradient: `makeGrad!(sp)`
- For hessian: `makeHess!(sp)`

Alternatively, to find a new coefficient vector using least squares given function values, sp.Values, 
and Basis Functions, sb.BF, then:

**Step 1**: Calculate Moore-Penrose Pseudo-Inverse for sb.BF: `make_pinvBFt!(sp)`

**Step 2**: Solve for new coefficient vector: `makeCoef!(sp)`

## Examples

```julia  
using Smolyak
mu = [2,2,2]
lb = -2.*ones(length(mu))
ub = 3.*ones(length(mu))
sg = SmolyakGrid(mu,lb,ub)
sb = SmolyakBasis(sg)
makeBasis!(sb)
sp = SmolyakPoly(sb)
makeValue!(sp,sb)
makeGrad!(sp,sb)
makeHess!(sp,sb)
```

For more detailed example see [Interpolation_Example.jl](./test/Interpolation_Example.jl).
	
"""
type SmolyakPoly{T}
	NumPts 		:: Int64 					# Number of points
	NumCoef 	:: Int64					# Number of elements in Coefficient Vector
	Coef 	 	:: Vector{Float64} 			# Coefficients of Smplyak Interpolant
	Value 	 	:: T						# Value of Interpolant at each of the sb.NumPts
	Grad		:: AA{Float64} 				# Gradient of Interpolant up to first sp.NumDerivArgs  at grid point n = 1:NumPts
	Hess	 	:: AM{Float64} 				# Hessian of Interpolant up to first sp.NumDerivArgs at grid point n = 1:NumPts
	NumDeriv 	:: Int64 					# Number of Derivative
	NumDerivArgs:: Int64 					# 1 to NumDerivArgs 
	pinvBFt		:: Matrix{Float64} 			# Transpose of Moore-Penrose Pseudo Inverse of sb.BF (for case where NumPts ≥ NumCoef)

end

function SmolyakPoly(sb::SmolyakBasis{Float64}; Coef::Vector=rand(sb.NumBF), NumDeriv::Int64=sb.NumDeriv, NumDerivArgs::Int64=sb.NumDerivArgs, NumPts::Int64=sb.NumPts)		
	
	Value = zeros(1)
	Grad = Vector{Float64}[zeros(NumDerivArgs) for n in 1:NumPts]
	Hess = Matrix{Float64}[zeros(NumDerivArgs,NumDerivArgs) for n in 1:NumPts]
	pinvBF = Array{Float64}(sb.NumBF,NumPts)

	return SmolyakPoly(NumPts, length(Coef), Coef, Value, Grad, Hess, NumDeriv, NumDerivArgs,pinvBF)
end

function SmolyakPoly(sb::SmolyakBasis{Vector{Float64}}; Coef::Vector=rand(sb.NumBF), NumDeriv::Int64=sb.NumDeriv, NumDerivArgs::Int64=sb.NumDerivArgs, NumPts::Int64=sb.NumPts)		
	
	Value = zeros(NumPts)
	Grad = Vector{Float64}[zeros(NumDerivArgs) for n in 1:NumPts]
	Hess = Matrix{Float64}[zeros(NumDerivArgs,NumDerivArgs) for n in 1:NumPts]
	pinvBF = Array{Float64}(sb.NumBF,NumPts)

	return SmolyakPoly(NumPts, length(Coef), Coef, Value, Grad, Hess, NumDeriv, NumDerivArgs, pinvBF)
end

# In place fn value update
function makeValue!(sp::SmolyakPoly{Float64},sb::SmolyakBasis;Coef::Vector{Float64}=sp.Coef)
	return sp.Value = dot(sb.BF[1],Coef)
end

# In place fn value update
function makeValue!(sp::SmolyakPoly{Vector{Float64}},sb::SmolyakBasis;Coef::Vector{Float64}=sp.Coef)
	for n in eachindex(sp.Value)
		sp.Value[n] = dot(sb.BF[n],Coef)
	end
	return sp.Value
end

# In place 1st Derivative Update
function makeGrad!(sp::SmolyakPoly,sb::SmolyakBasis;Coef::Vector{Float64}=sp.Coef,N::Int64=sp.NumDerivArgs)
	@inbounds for i in 1:N, n in eachindex(sp.Value)
		sp.Grad[n][i] = dot(sb.dBFdx[n][i],Coef)
	end
	return sp.Grad
end

# In place 2nd Derivative Update
function makeHess!(sp::SmolyakPoly,sb::SmolyakBasis;Coef::Vector{Float64}=sp.Coef,N::Int64=sp.NumDerivArgs)
	@inbounds for i in 1:N, j in i:N, n in eachindex(sp.Value)
		k = j-i+1
		sp.Hess[n][i,j] = dot(sb.d2BFdx2[n][i][k],Coef)
		!=(i,j) ? sp.Hess[n][j,i] = sp.Hess[n][i,j] : nothing
	end
	return sp.Hess
end

# Get Inverse of sb.BF to calculate sp.Coef by least squares
function make_pinvBFt!(sp::SmolyakPoly,sb::SmolyakBasis)
	if >=(sp.NumPts,sp.NumCoef)
		BF = Array{Float64}(sp.NumPts,sp.NumCoef)
		@inbounds for n in eachindex(sp.Value), p in eachindex(sp.Coef)
			BF[n,p] = sb.BF[n][p]
		end
		return sp.pinvBFt = pinv(BF)' # Transpose because matrix multiplication goes down cols
	else
		println("No single solution to least squares: NumPts < NumCoef. Not Created sp.pinvBFt")
	end
end

# Calculate New Coefficient using Least Squares with precalculate Moore-Penrose Pseudo-Inverse - may be more efficient than backslash if pinvBF not changing
function makeCoef!(sp::SmolyakPoly;f::Vector{Float64}=sp.Value)
	@inbounds for k in 1:sp.NumCoef
		s = 0.0
		for n in 1:sp.NumPts
		 	s += sp.pinvBFt[n,k]*f[n]
		end
		sp.Coef[k] = s
	end	
end

function show(io::IO, sp::SmolyakPoly)
	msg = "\n\tCreated Smolyak Interpolant on $(sp.NumPts) grid points:\n"
	if ==(sp.NumDeriv,0) msg *= "\t- No derivatives supplied. Do not call Gradient or Hessian.\n" end
	if ==(sp.NumDeriv,1) msg *= "\t- with gradient up to first $(sp.NumDerivArgs) arguments of Interpolant. Do not call Hessian.\n" end
	if ==(sp.NumDeriv,2) msg *= "\t- with gradient & hessian of first $(sp.NumDerivArgs) arguments of Interpolant.\n" end
	print(io, msg)
end





