
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
µ = vcat(2,3,2)					# Vector controlling polynomial accuracy in each dimension
Xlow = zeros(D) 				# Lower bounds over approximation domain
Xhigh = vcat(3.,2.,1.) 		# Upper bounds over approximation domain

# ------- Collocated function approximation ------------- #

# Grid of Smolyak grid: notice this automatically translates X → Z = [-1,1] using a linear transformation using chosen bounds
sg = SmolyakGrid(D,µ,Xlow,Xhigh)

#  Call matrix of Basis function over SmolyakGrid where call is SmolyakBasis(sg::SmolyakGrid, CalcInv::Bool=false, NumDeriv::Int64=1, SpOut::Bool=false)
sb = SmolyakBasis(sg,true,0) 

#= Note: This sets up basis matrix Ψ so that fhat = Ψ'*θ (i.e. each column i) - (i.e. each column of sb.Ψ is a vector of basis function evaluated at X[:,i])=#

# Get true value at on Smolyak Grid points sg.Xgrid which is D x NumBasisFun given µ
NBF = size(sg.xGrid,2)
f = Float64[]
for i in 1:NBF
	push!(f,truefun(sg.xGrid[:,i]))
end

# Get Coefficients in collocated case
θ = inv(sb.Ψ')*f

# Check the answer
fhat = sb.Ψ'*θ
figure(1)
scatter(f,fhat)

# ---------------------------
# Check off grid
NumObs = 100
X = hcat(rand(NumObs), 0.5*rand(NumObs), 0.25*rand(NumObs))';
f1 = Float64[]
for i in 1:NumObs
	push!(f1,truefun(X[:,i]))
end
sbX = SmolyakBasis(X,sg,true,0) # Smolyak Basis at D x N matrix of N grid function evaluation points
f1hat = sbX.Ψ'*θ 				# Approximated function
figure(2)
scatter(f1,f1hat)


#--- Now with more accuracy
µ = vcat(2,3,2)	
sg2 = SmolyakGrid(D,µ,Xlow,Xhigh)
sb2 = SmolyakBasis(sg,true,0)
NBF = size(sg2.xGrid,2)
f = Float64[]
for i in 1:NBF
	push!(f,truefun(sg2.xGrid[:,i]))
end 
# Get Coefficients in collocated case
θ = inv(sb2.Ψ')*f

# Check the answer
sb2X = SmolyakBasis(X,sg,true,0) # Smolyak Basis at D x N matrix of N grid function evaluation points
f1hat = sb2X.Ψ'*θ 				# Approximated function
figure(4)
scatter(f1,f1hat)
