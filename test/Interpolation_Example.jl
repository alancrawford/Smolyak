
#= Function approximation example 

File contains two example of ways to find interpolating coefficients 
of a function using the Smolyak Package. In both cases I do not use closed form formulas, instead i 
use lagrange interpolation or least squares (though could obviously do other methods in case 2).

1. The function is evaluated on the nodes of an Iso/Anisotrophic grid of D-dimensions and 
use lagrange interpolation to find the coefficients.

2. The points at which the function is evaluated is not necessarily the D-dimensional 
Smolyak Grid (i.e. perhaps a from a Simulation or picked in some other way), but you 
would still like to intepolate using the Smolyak Interpolant, f(D,µ). In this case, 
build Smolyak Grid, H(D,µ), and use it only to make the correpsonding Smolyak Basis 
functions on the grid points [Note 1]. Then I use least squares to find interpolating 
coefficients. Though, of course lots of options available in this last step - I know 
Maliar and Maliar (2014) have a nice section on this in their handbook chapter. 

[Note 1]: (Technically you don't need the grid, you only need bounds and accuracy index, µ, 
but i create grid anyway do decided to call Smolyak Basis using Smolyak Grid - lazy, I know... ).

I will do exercise 1 and 2 with two functions: 

1. fun1: can easily be handled with polynomial interpolation accuracy, i.e. µ=2. 
2. fun2: Make use of Anisotrophic accuracy index to efficiently get accurate interpolant.

In both cases i restrict domains to be finite (not in any special way, just randomly picked domain.

=#

#= --------- True Functions: fun1 and fun2 ----------- =# 

function fun1(x)
	return f = 1 + x[1] + x[2] + x[3] + x[1]*x[2] + x[2]*x[3] + x[1]^4 + x[2]^3 + x[3]^3
end

function fun2(x)
	return f = 1 + x[1] + x[2] + x[3] + x[1]*x[2] + x[2]*x[3] + x[1]^4 + x[2]^7 + x[3]^3
end

# Load Packages

using Smolyak, PyPlot

# ---- 1. Collocation: Lagrange Interpolation ------- #

# Setup Smolyak Grid 
D = 3; 								# Dimensionality of the problem 
µ = vcat(2,2,2);					# Vector controlling polynomial accuracy in each dimension
Xlow = zeros(D);					# Lower bounds over approximation domain
Xhigh = vcat(3.,2.,1.); 			# Upper bounds over approximation domain

# Smolyak grid
sg = SmolyakGrid(D,µ,Xlow,Xhigh)

# Basis Functions
sb = SmolyakBasis(sg,true,0) 

# Note: SmolyakBasis(sg::SmolyakGrid, CalcInv::Bool=false, NumDeriv::Int64=1, SpOut::Bool=false)

# Get values of fun1 and fun2 on the Smolyak Grid 
NBF = size(sg.xGrid,2)
f1 = Float64[]
f2 = Float64[]
for i in 1:NBF
	push!(f1,fun1(sg.xGrid[:,i]))
	push!(f2,fun2(sg.xGrid[:,i]))
end

# Get Coefficients in collocated case
θ1 = inv(sb.Ψ')*f1
θ2 = inv(sb.Ψ')*f2

#= ... or could use ready made inverse in SmolyakBasis 

α1 = sb.pinvΨ*f1
α2 = sb.pinvΨ*f1

Notes: 
	1. 	Called it pinvΨ not invΨ because will be Moore-Penrose Pseudo Inverse 
		when do least squares.
	2. 	Also, note that it is transposed - admittedly following corrected bug
		fixed yesterday (in case you already tried it!)

=#

# Check the answer
f1hat = sb.Ψ'*θ1
f2hat = sb.Ψ'*θ2

# Plot actual against fitted to check it worked - Excuse my rubbish PyPlotting... !
close()
figcount = 0

figcount += 1
figure(figcount)
title("Actual vs Interpolated Values on grid: fun1, µ=[2,2,2]")
scatter(f1,f1hat,marker= "o",color="red")
xlabel("Actual")
ylabel("Interpolated")

figcount += 1
figure(figcount)
title("Actual vs Interpolated Values on grid: fun2, µ=[2,2,2]")
scatter(f2,f2hat,marker="x",color="blue")
xlabel("Actual")
ylabel("Interpolated")


# To check interpolation use simulated points in domain and plot interpolated vs actual
NumObs = 100
X = hcat(3*rand(NumObs), 2*rand(NumObs), rand(NumObs))'; # Random points in domain
sim1 = Float64[]
sim2 = Float64[]
for i in 1:NumObs
	push!(sim1,fun1(X[:,i]))
	push!(sim2,fun2(X[:,i]))
end

# Setup Smolyak Basis function on X
sbX = SmolyakBasis(X,sg,true,0)
sim1hat = sbX.Ψ'*θ1
sim2hat = sbX.Ψ'*θ2

figcount += 1
figure(figcount)
title("Simulated Actual vs Interpolated Values: fun1, µ=[2,2,2]")
scatter(sim1,sim1hat,marker= "o",color="red")
xlabel("Actual")
ylabel("Interpolated")

figcount += 1
figure(figcount)
title("Simulated Actual vs Interpolated Values: fun2, µ=[2,2,2]")
scatter(sim2,sim2hat,marker= "x",color="blue")
xlabel("Actual")
ylabel("Interpolated")

# Obviously the interpolation on fun2 could be better, up accuracy for higher powered term, x[2]
µ_aniso = [2,3,2]
sg_aniso = SmolyakGrid(D,µ_aniso,Xlow,Xhigh)		# Anisotrophic Smolyak Grid
sb_aniso = SmolyakBasis(sg_aniso,true,0)			# Anisotrophic Basis Functions

# Evaluate fun2 on new grid
NBF_aniso = size(sg_aniso.xGrid,2)
f2_aniso = Float64[]
for i in 1:NBF_aniso
	push!(f2_aniso,fun2(sg_aniso.xGrid[:,i]))
end

# Get new interpolating coefficient on new grid using Lagrange Interpolation
θ2_aniso = inv(sb_aniso.Ψ')*f2_aniso

# Equivalently ... α2_aniso = sb_aniso.pinvΨ*f2_aniso
sbX_aniso = SmolyakBasis(X,sg_aniso,true,0)
sim2hat_aniso = sbX_aniso.Ψ'*θ2_aniso

figcount += 1
figure(figcount)
title("Simulated Actual vs Interpolated Values: fun2, µ=[2,3,2]")
scatter(sim2,sim2hat_aniso,marker= "x",color="blue")
xlabel("Actual")
ylabel("Interpolated")

# Plot show a lot better accuracy

# ---- 2. Regression: Least Squares Interpolation ------- #

# Simulate grid points & function values 
NumObs = 100;
X = hcat(3*rand(NumObs), 2*rand(NumObs), rand(NumObs))'; # Random points in domain
y1 = Float64[];
y2 = Float64[];
for i in 1:NumObs
	push!(y1,fun1(X[:,i]))
	push!(y2,fun2(X[:,i]))
end

# Smolyak Grid and Basis Matrix 
D = 3; 								# Dimensionality of the problem 
µ = vcat(2,2,2);					# Vector controlling polynomial accuracy in each dimension
Xlow = zeros(D); 					# Lower bounds over approximation domain
Xhigh = vcat(3.,2.,1.); 			# Upper bounds over approximation domain

sg = SmolyakGrid(D,µ,Xlow,Xhigh)	# Smolyak grid
sb_reg = SmolyakBasis(X,sg,true,0)  # Smolyak Basis on simulated Domain

# GET VALUES OF FUNCTIONS

α1 = sb_reg.pinvΨ*y1
α2 = sb_reg.pinvΨ*y2

# Check the answer
y1hat = sb_reg.Ψ'*α1
y2hat = sb_reg.Ψ'*α2

# Plot actual against fitted 

figcount += 1
figure(figcount)
title("Actual vs Interpolated Values on Simulated Points: fun1, µ=[2,2,2]")
scatter(y1,y1hat,marker= "o",color="red")
xlabel("Actual")
ylabel("Interpolated")

figcount += 1
figure(figcount)
title("Actual vs Interpolated Values on Simulated Points: fun2, µ=[2,2,2]")
scatter(y2,y2hat,marker="x",color="blue")
xlabel("Actual")
ylabel("Interpolated")

#= 	*********
	Again, zooming in you can see interpolation could be more accurate for case 2
	*********
=#

# ... as before can modify accuracy for fun2
µ_aniso = [2,3,2]
sg_aniso = SmolyakGrid(D,µ_aniso,Xlow,Xhigh)		# Anisotrophic Smolyak Grid
sb_reg_aniso = SmolyakBasis(X,sg_aniso,true,0)		# Anisotrophic Basis Functions
α2_aniso = sb_reg_aniso.pinvΨ*y2 					# Get new coefficients

# Get new interpolated points using µ = [2,3,2]
y2hat_aniso = sb_reg_aniso.Ψ'*α2_aniso

figcount += 1
figure(figcount)
title("Actual vs Interpolated Values on Simulated Points: fun2, µ=[2,3,2]")
scatter(y2,y2hat_aniso,marker="x",color="blue")
xlabel("Actual")
ylabel("Interpolated")