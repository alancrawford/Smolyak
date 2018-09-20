


module SmolyakTest
	using Smolyak, FactCheck

	for basis_fun_type in [:chebyshev, :spread], mu_ in 1:4

		facts("testing interpolation on grid with mu=$mu_ with basis fun = $(basis_fun_type)") do

			context("one dimension") do

				truefun(x) = 1.1 + x[1]^3

				mu = [mu_]
				xbnds = [[-2., rand()] for i in 1:length(mu)]

				# Smolyak Components
				sk = SmolyakKernel(mu, xbnds)
				sg = SmolyakGrid(sk)
				sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
				sp = SmolyakPoly(sb; NumDeriv=0)

				# Get true values of function at grid points
				W = truefun.(xgrid(sg))

				# Generate corresponding Smolyak Basis Functions
				BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));

				# Solve for the coefficients
				θ = BF\W

				# Update the coefficients
				coef!(θ, sp)  # Update coefficient in Smolyak Polynomial

				# Evaluate Smolayk poly on Smolyak Grid
				What = value(xgrid(sg) , sp)

				# Check maximum difference
				@fact maximum(abs, W-What) <1e-12 --> true 

			end

			context("2D") do
				truefun(x) = 1.1 + (x[1]-x[2]^2)^2 

				mu = [mu_,mu_]
				xbnds = [[-2., rand()] for i in 1:length(mu)]

				# Smolyak Components
				sk = SmolyakKernel(mu, xbnds)
				sg = SmolyakGrid(sk)
				sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
				sp = SmolyakPoly(sb; NumDeriv=0)

				# Get true values of function at grid points
				W = truefun.(xgrid(sg))

				# Generate corresponding Smolyak Basis Functions
				BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));

				# Solve for the coefficients
				θ = BF\W

				# Update the coefficients
				coef!(θ, sp)  # Update coefficient in Smolyak Polynomial

				# Evaluate Smolayk poly on Smolyak Grid
				What = value(xgrid(sg) , sp)

				# Check maximum difference
				@fact maximum(abs, W-What) <1e-12 --> true 

			end

			context("3D") do
				truefun(x) = 1.1 + (x[1]-x[2]^2)^2 + x[3]^2

				mu = [mu_,mu_,mu_]
				xbnds = [[-2., rand()] for i in 1:length(mu)]

				# Smolyak Components
				sk = SmolyakKernel(mu, xbnds)
				sg = SmolyakGrid(sk)
				sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
				sp = SmolyakPoly(sb; NumDeriv=0)

				# Get true values of function at grid points
				W = truefun.(xgrid(sg))

				# Generate corresponding Smolyak Basis Functions
				BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));

				# Solve for the coefficients
				θ = BF\W

				# Update the coefficients
				coef!(θ, sp)  # Update coefficient in Smolyak Polynomial

				# Evaluate Smolayk poly on Smolyak Grid
				What = value(xgrid(sg) , sp)

				# Check maximum difference
				@fact maximum(abs, W-What) <1e-12 --> true 
			end

			context("4D") do
				truefun(x) = 1.1 + (x[1]-x[2]^2)^2 + x[3]^2 - (x[3]+x[4])^2

				mu = [mu_,mu_,mu_,mu_]
				xbnds = [[-2., rand()] for i in 1:length(mu)]

				# Smolyak Components
				sk = SmolyakKernel(mu, xbnds)
				sg = SmolyakGrid(sk)
				sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
				sp = SmolyakPoly(sb; NumDeriv=0)

				# Get true values of function at grid points
				W = truefun.(xgrid(sg))

				# Generate corresponding Smolyak Basis Functions
				BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));

				# Solve for the coefficients
				θ = BF\W

				# Update the coefficients
				coef!(θ, sp)  # Update coefficient in Smolyak Polynomial

				# Evaluate Smolayk poly on Smolyak Grid
				What = value(xgrid(sg) , sp)

				# Check maximum difference
				@fact maximum(abs, W-What) <1e-12 --> true 
			end
		end
		

		facts("testing interpolation off grid with mu=$mu_") do

			# make a linear function to predict
			slopes = rand(4)

			# random point picker
			rpoint(lb,ub) = (ub - lb)*rand() + lb
			rpoint(bounds) = rpoint(bounds...)

			context("one dimension") do

				truefun(x) = 1.1 + slopes[1]*x[1]

				mu = [mu_]
				xbnds = [[-2., 12.] for i in 1:length(mu)]
				
				# Smolyak Components
				sk = SmolyakKernel(mu, xbnds)
				sg = SmolyakGrid(sk)
				sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
				sp = SmolyakPoly(sb; NumDeriv=0)

				# Get true values of function at grid points
				W = truefun.(xgrid(sg))

				# Solve for the coefficients
				BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));
				θ = BF\W

				# Update coefficient in Smolyak Polynomial
				coef!(θ, sp)  

				# make basis on random point
				NumObs = 10
				xx = [rpoint.(xbnds) for n in 1:NumObs]

				What = value(xx , sp)

				@fact What --> roughly(truefun.(xx))

			end

			context("2D") do
				truefun(x) = 1.1 + slopes[1]*x[1] - slopes[2]*x[2]

				mu = [mu_,mu_]
				xbnds = [[-2., 12.] for i in 1:length(mu)]
				
				# Smolyak Components
				sk = SmolyakKernel(mu, xbnds)
				sg = SmolyakGrid(sk)
				sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
				sp = SmolyakPoly(sb; NumDeriv=0)

				# Solve for the coefficients
				BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));
				θ = BF\truefun.(xgrid(sg))

				# Update coefficient in Smolyak Polynomial
				coef!(θ, sp)  

				# make basis on random point
				NumObs = 10
				xx = [rpoint.(xbnds) for n in 1:NumObs]

				What = value(xx , sp)

				@fact What --> roughly(truefun.(xx))
			end

			context("3D") do
				truefun(x) = 1.1 + slopes[1]*x[1] - slopes[2]*x[2] + slopes[3]*x[3]

				mu = [mu_,mu_,mu_]
				xbnds = [[-2., 12.] for i in 1:length(mu)]
				
				# Smolyak Components
				sk = SmolyakKernel(mu, xbnds)
				sg = SmolyakGrid(sk)
				sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
				sp = SmolyakPoly(sb; NumDeriv=0)

				# Solve for the coefficients
				BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));
				θ = BF\truefun.(xgrid(sg))

				# Update coefficient in Smolyak Polynomial
				coef!(θ, sp)  

				# make basis on random point
				NumObs = 10
				xx = [rpoint.(xbnds) for n in 1:NumObs]

				What = value(xx , sp)

				@fact What --> roughly(truefun.(xx))
			end

			context("4D") do
				truefun(x) = 1.1 + slopes[1]*x[1] - slopes[2]*x[2] + slopes[3]*x[3] * slopes[4] * x[4]

				mu = [mu_,mu_,mu_,mu_]
				xbnds = [[-2., 12.] for i in 1:length(mu)]
				
				# Smolyak Components
				sk = SmolyakKernel(mu, xbnds)
				sg = SmolyakGrid(sk)
				sb = SmolyakBasis(basis_fun_type, sk; NumDeriv=0)
				sp = SmolyakPoly(sb; NumDeriv=0)

				# Solve for the coefficients
				BF = VVtoMatrix(BasisFunctions(xgrid(sg), sb));
				θ = BF\truefun.(xgrid(sg))

				# Update coefficient in Smolyak Polynomial
				coef!(θ, sp)  

				# make basis on random point
				NumObs = 10
				xx = [rpoint.(xbnds) for n in 1:NumObs]

				What = value(xx , sp)

				if mu_ == 1
					println("approximation level mu=1 is too low in 4D with multiplicative component.")
					for i in 1:NumObs
						@fact What --> !roughly(truefun.(xx))
					end

				else
					for i in 1:NumObs
						@fact What --> roughly(truefun.(xx))
					end
				end

			end
		end
	end
end