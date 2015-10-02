
#= ******************************* =#
#= Function approximation examples =#
#= ******************************* =#

#=
	**********************************
	Collocated function approximation
	**********************************

Goal: 

	Find Smolyak Interpolant Approximating truefun(x) on 
	hypercube: x[1]∈[0,3] x[2]∈[0,2], x[3]∈[0,1] using
	Collocation on Smolyak Grid Points.

Process:

1. Setup a Smolyak Grid 
2. Create Basis Functions
3. Evaluate truefun on Smolyak Grid
4. Find Interpolating Coefficients.
5. Check the Answer!
6. How accurate is the answer?

=# 

#= 		**********************
	 	1. Create Smolyak Grid 
		**********************
	
This automatically translates X → Z = [-1,1] using a linear transformation using chosen bounds

Call:
	SmolyakGrid(D::Int64, mu::ScalarOrVec{Int64},lb::Vector{Float64}=-1*ones(Float64,D), ub::Vector{Float64}=ones(Float64,D))

  Summary:

  type Smolyak.SmolyakGrid <: Any

  Fields:

  D         :: Int64
  mu        :: Union{Array{Int64,1},Int64}
  NumGrdPts :: Int64
  lb        :: Array{Float64,1}
  ub        :: Array{Float64,1}
  zGrid     :: Array{Array{Float64,1},1}
  xGrid     :: Array{Array{Float64,1},1}
  Binds     :: Array{Array{Int64,1},1}

=#

# *************** Define True FunctionL truefun(x) ********************* 
function truefun(x)
	return f = 1 + x[1] + x[2] + x[3] + x[1]*x[2] + x[2]*x[3] + x[1]^4 + x[2]^5 + x[3]^5
end

# Setup Corresponding SmolyakGrid =#

using Smolyak, PyPlot
D = 3 								# Dimensionality of the problem 
mu = [2,2,2]						# Vector controlling polynomial accuracy in each dimension
Xlow = zeros(D) 					# Lower bounds over approximation domain
Xhigh = collect(3.:-1:1.) 			# Upper bounds over approximation domain
sg = SmolyakGrid(D,mu,Xlow,Xhigh) 	# Set up Smolyak Grid


#= 	***********************************************************
	2. Create Smolyak Basis Function of Interpolating Polynomial 
	***********************************************************
	
Create Smolyak Basis Function of Interpolating Polynomial 

When formulating Basis Functions on the Smolyak Grid points call:

	sb = SmolyakBasis(sg::SmolyakGrid,NumDeriv::Int64=2,NumDerivArgs::Int64=sg.D)

Alternatively, if you wish to calcuate Smolyak basis functions on D dimesnsions
with mu controlling the order of basis functions in each dimension on a vector 
(or array of vectors), x, with n=1:length(x) and D=length(x[n]), call:

	sb = SmolyakBasis(x::VecOrAA{Float64},sg::SmolyakGrid,NumDeriv::Int64=2,NumDerivArgs::Int64=sg.D)

Once, setup the parameters and allocated memory for Smolyak Basis, and then call

	makeBasis!(sb)

to create it.


  Summary:

  type Smolyak.SmolyakBasis <: Any

  Fields:

  D            :: Int64
  mu           :: Union{Array{Int64,1},Int64}
  lb           :: Array{Float64,1}
  ub           :: Array{Float64,1}
  Binds        :: Array{Array{Int64,1},1}
  NumPts       :: Int64
  NumBF        :: Int64
  NumDeriv     :: Int64
  NumDerivArgs :: Int64
  max_order    :: Int64
  x            :: Union{Array{Array{Float64,1},1},Array{Float64,1}}
  z            :: Union{Array{Array{Float64,1},1},Array{Float64,1}}
  T            :: Array{Float64,2}
  dT           :: Array{Float64,2}
  d2T          :: Array{Float64,2}
  BF           :: Array{Array{Float64,1},1}
  dBFdz        :: Array{Array{Array{Float64,1},1},1}
  d2BFdz2      :: Array{Array{Array{Array{Float64,1},1},1},1}
  dzdx         :: Array{Float64,1}
  d2zdx2       :: Array{Float64,1}
  dBFdx        :: Array{Array{Array{Float64,1},1},1}
  d2BFdx2      :: Array{Array{Array{Array{Float64,1},1},1},1}
=#

NumDeriv=0
NumDerivArgs=sg.D
sb = SmolyakBasis(sg,NumDeriv,NumDerivArgs) 	# Create Smolyak Basis Type
makeBasis!(sb) 									# This sets up basis matrix: sb.BF



#= 	***********************************
	3. Evaluate truefun on Smolyak Grid
	***********************************

Two steps:
3.1. Create Smolyak Polynomial type to read true function values into
3.2. Insert true values of function in SmolyakPoly.Value

To setup SmolyakPoly call: 

	SmolyakPoly(sb::SmolyakBasis, Coef::Vector=rand(sb.NumBF), NumDeriv::Int64=sb.NumDeriv, NumDerivArgs::Int64=sb.NumDerivArgs, NumPts::Int64=sb.NumPts)		

  Summary:

  type Smolyak.SmolyakPoly <: Any

  Fields:

  NumPts       :: Int64
  NumCoef      :: Int64
  Coef         :: Array{Float64,1}
  Value        :: Union{Array{Float64,1},Float64}
  Grad         :: Array{Array{Float64,1},1}
  Hess         :: Array{Array{Float64,2},1}
  NumDeriv     :: Int64
  NumDerivArgs :: Int64
  pinvBFt      :: Array{Float64,2}

If do not have coefficient, constructor function allocates conformable memory.

=#

sp = SmolyakPoly(sb) 				# Create SmolyakPoly Type  
for i in 1:sb.NumPts
	sp.Value[i] = truefun(sb.x[i]) 	# Use truefun to input correct values of function into
end

#= 	***********************************
	4. Get Interpolating Coefficients
	***********************************
=#
get_pinvBFt!(sp,sb)		# Get Moore-Penrose Pseudo Inverse Matrix of Smolyak Basis
getCoef!(sp) 			# Get Collocating Coefficients.

#=  **********************************
	5. Check Answer from Smolyak Package
	**************************************

1. Create Matrix BF
2. Solve system CoefTrue = BF\sp.Value
3. Compare CoefTrue and sp.Coef

=#

# 1. Create Matrix BF
BF = Array{Float64}(sb.NumPts,sb.NumBF)
for n in 1:sb.NumPts, p in 1:sb.NumBF
	BF[n,p] = sb.BF[n][p]
end

# 2. Calculate Coefficient Vector for Comparison
CoefTrue = BF\sp.Value

# 3. Do the Test!
isapprox(CoefTrue,sp.Coef) ? println("Correct Answer: CoefTrue == sp.Coef") : println("Oh no, something went wrong!")


#=  **********************************
	5. Check Answer from Smolyak Package
	**************************************

1. Create Matrix BF
2. Solve system CoefTrue = BF\sp.Value
3. Compare CoefTrue and sp.Coef

=#

NumObs = 1000
X = Vector{Float64}[Array{Float64}(D) for i = 1:NumObs]
f = Array{Float64}(NumObs)
for i in 1:NumObs
	X[i] = [ Xlow[d]+( Xhigh[d]- Xlow[d])*rand() for d in 1:D]
	f[i] = truefun(X[i])
end

sbX = SmolyakBasis(X,D,mu,Xlow,Xhigh,NumDeriv,NumDerivArgs)
makeBasis!(sbX)
spX = SmolyakPoly(sbX)
copy!(spX.Coef,sp.Coef) 		# Copy over Coefficient from Collocation
getValue!(spX,sbX) 				# Interpolated Values

# Error Checking
max_abs_error = 0.0
mean_abs_error = 0.0
for i in 1:NumObs
	max_abs_error = max(max_abs_error, abs(spX.Value[i]/f[i]-1))
	mean_abs_error += abs(spX.Value[i]/f[i]-1)
end
println("\nWhen D = $D, and µ = $mu then:")
println("Max Abs. Error is log 10 units = $(round(log10(max_abs_error),4))")
println("Mean Abs. Error is log 10 units = $(round(log10(mean_abs_error/NumObs),4))")
println("Now change mu and see how accuracy is affected. Also change the truefun(x) formula.... ")