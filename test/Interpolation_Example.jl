
#= Function approximation example =#

# True Function to generate the data 
function truefun(x)
	return f = 1 + x[1] + x[2] + x[3] + x[1]*x[2] + x[2]*x[3] + x[1]^4 + x[2]^9 + x[3]^5
end

#= Function Approximation Code =#

# Set up value of function to interpolate
NumObs = 100
X = hcat(rand(NumObs), 0.5*rand(NumObs), 0.25*rand(NumObs))';
f = Float64[]
for i in 1:NumObs
	push!(f,truefun(X[:,i]))
end

# Setup function approximation

using Smolyak, PyPlot

# Setup Grid Nodes
D = 3 							# Dimensionality of the problem 
mu = vcat(2,2,2)				# Vector controlling polynomial accuracy in each dimension
Xlow = zeros(D) 				# Lower bounds over approximation domain
Xhigh = vcat(3.,2.,1.) 			# Upper bounds over approximation domain

# ------- Collocated function approximation ------------- #

# Grid of Smolyak grid: notice this automatically translates X → Z = [-1,1] using a linear transformation using chosen bounds
sg = SmolyakGrid(D,mu,Xlow,Xhigh)

#  Call matrix of Basis function over SmolyakGrid where call is SmolyakBasis(sg::SmolyakGrid, NumDeriv::Int64=1)
sb = SmolyakBasis(sg,0,true) 	#= Create Smolyak Basis Type =#
makeBF!(sb) 					#= This sets up basis matrix sb.BF =#

# Get true value at on Smolyak Grid points sg.Xgrid which is D x NumBasisFun given mu
f = Float64[]
for i in 1:sb.NumPts
	push!(f,truefun(sb.x[:,i]))
end

# Get Coefficients in collocated case
coef = sb.pinvBF'*f

# Check the answer
sp = SmolyakPoly(sb,coef)
f-sp.f 
figure(1)
scatter(f,sp.f)

# ---------------------------
# Check off grid
NumObs = 100
X = hcat(rand(NumObs), 0.5*rand(NumObs), 0.25*rand(NumObs))';
fX = Float64[]
for i in 1:NumObs
	push!(fX,truefun(X[:,i]))
end
sbX = SmolyakBasis(X,sg,0,true) # Smolyak Basis at D x N matrix of N grid function evaluation points
 
spX = SmolyakPoly(sbX,coef) 				# Approximated function
fX-spX.f
figure(2)
scatter(fX,spX.f)


#--- Now with more accuracy
mu = vcat(2,3,2)	
sg2 = SmolyakGrid(D,mu,Xlow,Xhigh)
sb2 = SmolyakBasis(sg2,0,true)
makeBF!(sb2) 
f2 = Float64[]
for i in 1:sb2.NumPts
	push!(f2,truefun(sb2.x[:,i]))
end 
# Get Coefficients in collocated case
coef2 = sb2.pinvBF'*f2

# Check the answer
f2X = Float64[]
for i in 1:NumObs
	push!(f2X,truefun(X[:,i]))
end
sb2X = SmolyakBasis(X,sg2,0,true) # Smolyak Basis at D x N matrix of N grid function evaluation points
makeBF!(sb2X)
sp2X = SmolyakPoly(sb2X,coef2) 				# Approximated function
f2X-sp2X.f
figure(3)
scatter(f2,sp2X.f)
