#= Allocations!
sp.f = Vector{Float64}(sb.NumPts)
sp.Grad=Vector{Float64}[ones(Float64,sb.D) for i in 1:sb.NumPts]
sp.Hess=Matrix{Float64}[ones(Float64,sb.D,sb.D) for i in 1:sb.NumPts]

# Then do looped matrix multiplication

for n in ... and so on... 
s = 0.0
for p = 1:sb.NumBF
	s += sb.d2BFdx[n][i][j][p]*θ[p] 
end
sp.Hess[n][i,j] = s

=#

#= Smolyak Interpolant =#

#= ************************ =#
#= Smolyak Polynomial type  =#
#= ************************ =#


type SmolyakPoly
	NumPts 		:: Int64 					# Number of points
	NumCoef 	:: Int64					# Number of elements in Coefficient Vector
	Coef 	 	:: Vector{Float64} 			# Coefficients of Smplyak Interpolant
	Value 	 	:: ScalarOrVec{Float64} 	# Value of Interpolant at each of the sb.NumPts
	Grad		:: AA{Float64} 				# Gradient of Interpolant up to first sp.NumDerivArgs  at grid point n = 1:NumPts
	Hess	 	:: AM{Float64} 				# Hessian of Interpolant up to first sp.NumDerivArgs at grid point n = 1:NumPts
	NumDeriv 	:: Int64 					# Number of Derivative
	NumDerivArgs:: Int64 					# 1 to NumDerivArgs 
	pinvBFt		:: Matrix{Float64} 			# Transpose of Moore-Penrose Pseudo Inverse of sb.BF (for case where NumPts ≥ NumCoef)

	function SmolyakPoly(sb::SmolyakBasis, Coef::Vector=rand(sb.NumBF), NumDeriv::Int64=sb.NumDeriv, NumDerivArgs::Int64=sb.NumDerivArgs, NumPts::Int64=sb.NumPts)		
		
		Value = Vector{Float64}(NumPts)
		Grad = Vector{Float64}[Array{Float64}(NumDerivArgs) for n in 1:NumPts]
		Hess = Matrix{Float64}[zeros(NumDerivArgs,NumDerivArgs) for n in 1:NumPts]
		pinvBF = Array{Float64}(sb.NumBF,NumPts)

		new(NumPts, length(Coef), Coef, Value, Grad, Hess, NumDeriv, NumDerivArgs)
	end

end

# Vector to multiply Basis Function vectors with coefficients
function VecMult(X::Vector{Float64},Coef::Vector{Float64},s::Float64=0.0)
	for i in eachindex(X)
		s += X[i]*Coef[i]
	end
	return s
end

# In place fn value update
function getValue!(sp::SmolyakPoly,sb::SmolyakBasis,Coef::Vector{Float64}=sp.Coef)
	for n in 1:sb.NumPts
		sp.Value[n] = VecMult(sb.BF[n],Coef)
	end
end

# In place 1st Derivative Update
function getGrad!(sp::SmolyakPoly,sb::SmolyakBasis,Coef::Vector{Float64}=sp.Coef,N::Int64=sp.NumDerivArgs)
	for i in 1:N, n in 1:sb.NumPts,
		sp.Grad[n][i] = VecMult(sb.dBFdx[n][i],Coef)
	end
end

# In place 2nd Derivative Update
function getHess!(sp::SmolyakPoly,sb::SmolyakBasis,Coef::Vector{Float64}=sp.Coef,N::Int64=sp.NumDerivArgs)
	for i in 1:N, j in i:N, n in 1:sb.NumPts,
		k = j-i+1
		sp.Hess[n][i,j] = VecMult(sb.d2BFdx2[n][i][k],Coef)
		!=(i,j) ? sp.Hess[n][j,i] = sp.Hess[n][i,j] : nothing
	end
end

# Get Inverse of sb.BF to calculate sp.Coef by least squares
function get_pinvBFt!(sp::SmolyakPoly,sb::SmolyakBasis)
	if >=(sp.NumPts,sp.NumCoef)
		BF = Array{Float64}(sp.NumPts,sp.NumCoef)
		for n in eachindex(sp.Value), p in eachindex(sp.Coef)
			BF[n,p] = sb.BF[n][p]
		end
		return sp.pinvBFt = pinv(BF)' # Transpose because matrix multiplication goes down cols
	else
		println("No single solution to least squares: NumPts < NumCoef. Not Created sp.pinvBFt")
	end
end

# Calculate New Coefficient using Least Squares with precalculate Moore-Penrose Pseudo-Inverse - may be more efficient than backslash if pinvBF not changing
function getCoef!(sp::SmolyakPoly,f::Vector{Float64}=sp.Value)
	for k in 1:sp.NumCoef
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





