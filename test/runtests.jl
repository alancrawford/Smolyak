
using Smolyak
using Test

# make a linear function to predict
global slopes = rand(4)

for basis_fun_type in [:chebyshev, :spread], mu_ in 1:4
	@testset "testing interpolation on grid with mu=$mu_ with basis fun = $(basis_fun_type)" begin

		println("1D") 

		truefun1(x) = 1.1 + x[1]^3

		mu = [mu_]
		xbnds = [[-2., rand()] for i in 1:length(mu)]

		# Smolyak Components
		sk = SmolyakKernel(mu, xbnds)
		sg = SmolyakGrid(sk)
		sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
		sp = SmolyakPoly(sb; NumDeriv=0)

		# Get true values of function at grid points
		W = truefun1.(xgrid(sg))

		# Generate corresponding Smolyak Basis Functions
		BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));

		# Solve for the coefficients
		θ = BF\W

		# Update the coefficients
		coef!(θ, sp)  # Update coefficient in Smolyak Polynomial

		# Evaluate Smolayk poly on Smolyak Grid
		What = value(xgrid(sg) , sp)

		# Check maximum difference
		@test maximum(abs, W-What) < 1e-12  
	
		println("2D") 
		truefun2(x) = 1.1 + (x[1]-x[2]^2)^2 

		mu = [mu_,mu_]
		xbnds = [[-2., rand()] for i in 1:length(mu)]

		# Smolyak Components
		sk = SmolyakKernel(mu, xbnds)
		sg = SmolyakGrid(sk)
		sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
		sp = SmolyakPoly(sb; NumDeriv=0)

		# Get true values of function at grid points
		W = truefun2.(xgrid(sg))

		# Generate corresponding Smolyak Basis Functions
		BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));

		# Solve for the coefficients
		θ = BF\W

		# Update the coefficients
		coef!(θ, sp)  # Update coefficient in Smolyak Polynomial

		# Evaluate Smolayk poly on Smolyak Grid
		What = value(xgrid(sg) , sp)

		# Check maximum difference
		@test maximum(abs, W-What) <1e-12  

		println("3D") 
		truefun3(x) = 1.1 + (x[1]-x[2]^2)^2 + x[3]^2

		mu = [mu_,mu_,mu_]
		xbnds = [[-2., rand()] for i in 1:length(mu)]

		# Smolyak Components
		sk = SmolyakKernel(mu, xbnds)
		sg = SmolyakGrid(sk)
		sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
		sp = SmolyakPoly(sb; NumDeriv=0)

		# Get true values of function at grid points
		W = truefun3.(xgrid(sg))

		# Generate corresponding Smolyak Basis Functions
		BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));

		# Solve for the coefficients
		θ = BF\W

		# Update the coefficients
		coef!(θ, sp)  # Update coefficient in Smolyak Polynomial

		# Evaluate Smolayk poly on Smolyak Grid
		What = value(xgrid(sg) , sp)

		# Check maximum difference
		@test maximum(abs, W-What) <1e-12

		println("4D")
		truefun4(x) = 1.1 + (x[1]-x[2]^2)^2 + x[3]^2 - (x[3]+x[4])^2

		mu = [mu_,mu_,mu_,mu_]
		xbnds = [[-2., rand()] for i in 1:length(mu)]

		# Smolyak Components
		sk = SmolyakKernel(mu, xbnds)
		sg = SmolyakGrid(sk)
		sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
		sp = SmolyakPoly(sb; NumDeriv=0)

		# Get true values of function at grid points
		W = truefun4.(xgrid(sg))

		# Generate corresponding Smolyak Basis Functions
		BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));

		# Solve for the coefficients
		θ = BF\W

		# Update the coefficients
		coef!(θ, sp)  # Update coefficient in Smolyak Polynomial

		# Evaluate Smolayk poly on Smolyak Grid
		What = value(xgrid(sg) , sp)

		# Check maximum difference
		@test maximum(abs, W-What) <1e-12 
	end
end
for basis_fun_type in [:chebyshev, :spread], mu_ in 1:4

	@testset "testing interpolation off grid with mu=$(mu_)" begin

		# random point picker
		rpoint(lb,ub) = (ub - lb)*rand() + lb
		rpoint(bounds) = rpoint(bounds...)

		println("one dimension")

		function truefun1(x)
			return 1.1 + slopes[1]*x[1]
		end 

		mu = [mu_]
		xbnds = [[-2., 12.] for i in 1:length(mu)]

		# Smolyak Components
		sk = SmolyakKernel(mu, xbnds)
		sg = SmolyakGrid(sk)
		sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
		sp = SmolyakPoly(sb; NumDeriv=0)

		# Get true values of function at grid points
		W = truefun1.(xgrid(sg))

		# Solve for the coefficients
		BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));
		θ = BF\W

		# Update coefficient in Smolyak Polynomial
		coef!(θ, sp)  

		# make basis on random point
		NumObs = 10
		xx = [rpoint.(xbnds) for n in 1:NumObs]

		What = value(xx , sp)

		@test isapprox(What, truefun1.(xx), atol=1e-6)

		println("2D")
		function truefun2(x)
			return 1.1 + slopes[1]*x[1] - slopes[2]*x[2]
		end
		mu = [mu_,mu_]
		xbnds = [[-2., 12.] for i in 1:length(mu)]

		# Smolyak Components
		sk = SmolyakKernel(mu, xbnds)
		sg = SmolyakGrid(sk)
		sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
		sp = SmolyakPoly(sb; NumDeriv=0)

		# Solve for the coefficients
		BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));
		θ = BF\truefun2.(xgrid(sg))

		# Update coefficient in Smolyak Polynomial
		coef!(θ, sp)  

		# make basis on random point
		NumObs = 10
		xx = [rpoint.(xbnds) for n in 1:NumObs]

		What = value(xx , sp)

		@test isapprox(What, truefun2.(xx), atol=1e-6)	

		println("3D") 
		function truefun3(x)
			return 1.1 + slopes[1]*x[1] - slopes[2]*x[2] + slopes[3]*x[3]
		end

		mu = [mu_,mu_,mu_]
		xbnds = [[-2., 12.] for i in 1:length(mu)]

		# Smolyak Components
		sk = SmolyakKernel(mu, xbnds)
		sg = SmolyakGrid(sk)
		sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
		sp = SmolyakPoly(sb; NumDeriv=0)

		# Solve for the coefficients
		BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));
		θ = BF\truefun3.(xgrid(sg))

		# Update coefficient in Smolyak Polynomial
		coef!(θ, sp)  

		# make basis on random point
		NumObs = 10
		xx = [rpoint.(xbnds) for n in 1:NumObs]

		What = value(xx , sp)

		@test isapprox(What, truefun3.(xx), atol=1e-6)

		println("4D")
		function truefun4(x)
			return  1.1 + slopes[1]*x[1] - slopes[2]*x[2] + slopes[3]*x[3] * slopes[4] * x[4]
		end
		mu = [mu_,mu_,mu_,mu_]
		xbnds = [[-2., 12.] for i in 1:length(mu)]

		# Smolyak Components
		sk = SmolyakKernel(mu, xbnds)
		sg = SmolyakGrid(sk)
		sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
		sp = SmolyakPoly(sb; NumDeriv=0)

		# Solve for the coefficients
		BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));
		θ = BF\truefun4.(xgrid(sg))

		# Update coefficient in Smolyak Polynomial
		coef!(θ, sp)  

		# make basis on random point
		NumObs = 10
		xx = [rpoint.(xbnds) for n in 1:NumObs]

		What = value(xx , sp)

		if mu_ == 1
			println("approximation level mu=1 is too low in 4D with multiplicative component.")
			for i in 1:NumObs
				#@test !isapprox(What, truefun.(xx), atol=1e-6)
			end

		else
			for i in 1:NumObs
				@test isapprox(What, truefun4.(xx), atol=1e-6)
			end
		end
	end
end
