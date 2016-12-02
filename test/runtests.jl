


module SmolyakTest
	using Smolyak, FactCheck

	for mu_ in 1:4

		facts("testing interpolation on grid with mu=$mu_") do

			context("one dimension") do

				truefun(x) = 1.1 + x[1]^3

				mu = [mu_]
				lb = -2.*ones(length(mu))
				ub = rand(length(mu))
				sg = SmolyakGrid(mu,lb,ub)
				sb = SmolyakBasis(sg)
				makeBasis!(sb)
				sp = SmolyakPoly(sb)

				for i in 1:sb.NumPts
					sp.Value[i] = truefun(sg.xGrid[i]) 	# Use truefun to input correct values of function into
				end

				@fact maxabs(makeValue!(sp,sb) - sp.Value) < 1e-16 --> true
			end

			context("2D") do
				truefun(x) = 1.1 + (x[1]-x[2]^2)^2 

				mu = [mu_,mu_]
				lb = -2.*ones(length(mu))
				ub = rand(length(mu))
				sg = SmolyakGrid(mu,lb,ub)
				sb = SmolyakBasis(sg)
				makeBasis!(sb)
				sp = SmolyakPoly(sb)

				for i in 1:sb.NumPts
					sp.Value[i] = truefun(sg.xGrid[i]) 	# Use truefun to input correct values of function into
				end

				@fact maxabs(makeValue!(sp,sb) - sp.Value) < 1e-16 --> true

			end

			context("3D") do
				truefun(x) = 1.1 + (x[1]-x[2]^2)^2 + x[3]^2

				mu = [mu_,mu_,mu_]
				lb = -2.*ones(length(mu))
				ub = rand(length(mu))
				sg = SmolyakGrid(mu,lb,ub)
				sb = SmolyakBasis(sg)
				makeBasis!(sb)
				sp = SmolyakPoly(sb)

				for i in 1:sb.NumPts
					sp.Value[i] = truefun(sg.xGrid[i]) 	# Use truefun to input correct values of function into
				end

				@fact maxabs(makeValue!(sp,sb) - sp.Value) < 1e-16 --> true

			end

			context("4D") do
				truefun(x) = 1.1 + (x[1]-x[2]^2)^2 + x[3]^2 - (x[3]+x[4])^2

				mu = [mu_,mu_,mu_,mu_]
				lb = -2.*ones(length(mu))
				ub = rand(length(mu))
				sg = SmolyakGrid(mu,lb,ub)
				sb = SmolyakBasis(sg)
				makeBasis!(sb)
				sp = SmolyakPoly(sb)

				for i in 1:sb.NumPts
					sp.Value[i] = truefun(sg.xGrid[i]) 	# Use truefun to input correct values of function into
				end

				@fact maxabs(makeValue!(sp,sb) - sp.Value) < 1e-16 --> true

			end
		end
		

		facts("testing interpolation off grid with mu=$mu_") do

			# make a linear function to predict
			slopes = rand(4)

			# random point picker
			rpoint(lb,ub) = (ub - lb)*rand() + lb

			context("one dimension") do

				truefun(x) = 1.1 + slopes[1]*x[1]

				mu = [mu_]
				lb = -2.*ones(length(mu))
				ub = 12.*ones(length(mu))
				sg = SmolyakGrid(mu,lb,ub)
				sb = SmolyakBasis(sg)
				makeBasis!(sb)
				sp = SmolyakPoly(sb)
				for i in 1:sb.NumPts
					sp.Value[i] = truefun(sg.xGrid[i]) 	# Use truefun to input correct values of function into
				end
				make_pinvBFt!(sp,sb)		
				makeCoef!(sp) 		

				# make basis on random point
				NumObs = 10
				X = Vector{Float64}[Array{Float64}(1) for i = 1:NumObs]
				for i in 1:NumObs
					X[i] = collect(rpoint(lb,ub))
				end
				sbX = SmolyakBasis(X,mu,lb,ub;NumDeriv=0)
				makeBasis!(sbX)
				spX = SmolyakPoly(sbX)
				copy!(spX.Coef,sp.Coef) 		
				makeValue!(spX,sbX) # Interpolated Values

				for i in 1:NumObs
					@fact spX.Value[i] --> roughly(truefun(X[i]))
				end
			end

			context("2D") do
				truefun(x) = 1.1 + slopes[1]*x[1] - slopes[2]*x[2]

				mu = [mu_,mu_]
				D = length(mu)
				lb = -2.*ones(length(mu))
				ub = 12.*ones(length(mu))
				sg = SmolyakGrid(mu,lb,ub)
				sb = SmolyakBasis(sg)
				makeBasis!(sb)
				sp = SmolyakPoly(sb)
				for i in 1:sb.NumPts
					sp.Value[i] = truefun(sg.xGrid[i]) 	# Use truefun to input correct values of function into
				end
				make_pinvBFt!(sp,sb)		
				makeCoef!(sp) 		

				# make basis on random point
				NumObs = 10
				X = Vector{Float64}[ Float64[lb[d]+( ub[d]- lb[d])*rand() for d in 1:D] for i in 1:NumObs]
				sbX = SmolyakBasis(X,mu,lb,ub;NumDeriv=0)
				makeBasis!(sbX)
				spX = SmolyakPoly(sbX)
				copy!(spX.Coef,sp.Coef) 		
				makeValue!(spX,sbX) # Interpolated Values

				for i in 1:NumObs
					@fact spX.Value[i] --> roughly(truefun(X[i]))
				end

			end

			context("3D") do
				truefun(x) = 1.1 + slopes[1]*x[1] - slopes[2]*x[2] + slopes[3]*x[3]

				mu = [mu_,mu_,mu_]
				D = length(mu)
				lb = -2.*ones(length(mu))
				ub = 12.*ones(length(mu))
				sg = SmolyakGrid(mu,lb,ub)
				sb = SmolyakBasis(sg)
				makeBasis!(sb)
				sp = SmolyakPoly(sb)
				for i in 1:sb.NumPts
					sp.Value[i] = truefun(sg.xGrid[i]) 	# Use truefun to input correct values of function into
				end
				make_pinvBFt!(sp,sb)		
				makeCoef!(sp) 		

				# make basis on random point
				NumObs = 10
				X = Vector{Float64}[ Float64[lb[d]+( ub[d]- lb[d])*rand() for d in 1:D] for i in 1:NumObs]
				sbX = SmolyakBasis(X,mu,lb,ub;NumDeriv=0)
				makeBasis!(sbX)
				spX = SmolyakPoly(sbX)
				copy!(spX.Coef,sp.Coef) 		
				makeValue!(spX,sbX) # Interpolated Values

				for i in 1:NumObs
					@fact spX.Value[i] --> roughly(truefun(X[i]))
				end
			end

			context("4D") do
				truefun(x) = 1.1 + slopes[1]*x[1] - slopes[2]*x[2] + slopes[3]*x[3] * slopes[4] * x[4]

				mu = [mu_,mu_,mu_,mu_]
				D = length(mu)
				lb = -2.*ones(length(mu))
				ub = 12.*ones(length(mu))
				sg = SmolyakGrid(mu,lb,ub)
				sb = SmolyakBasis(sg)
				makeBasis!(sb)
				sp = SmolyakPoly(sb)
				for i in 1:sb.NumPts
					sp.Value[i] = truefun(sg.xGrid[i]) 	# Use truefun to input correct values of function into
				end
				make_pinvBFt!(sp,sb)		
				makeCoef!(sp) 		

				# make basis on random point
				NumObs = 10
				X = Vector{Float64}[ Float64[lb[d]+( ub[d]- lb[d])*rand() for d in 1:D] for i in 1:NumObs]
				sbX = SmolyakBasis(X,mu,lb,ub;NumDeriv=0)
				makeBasis!(sbX)
				spX = SmolyakPoly(sbX)
				copy!(spX.Coef,sp.Coef) 		
				makeValue!(spX,sbX) # Interpolated Values

				if mu_ == 1
					println("approximation level mu=1 is too low in 4D with multiplicative component.")
					for i in 1:NumObs
						@fact spX.Value[i] --> not(roughly(truefun(X[i])))
					end

				else
					for i in 1:NumObs
						@fact spX.Value[i] --> roughly(truefun(X[i]))
					end
				end

			end
		end
	end
end