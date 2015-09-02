using Smolyak
D = 3
µ = 2
lb = -2ones(D)
ub = 2ones(D)
sg = SmolyakGrid(D,µ,lb,ub)
zgrid = round(cos(sg.thetaGrid),14)

@time T, ∂T, ∂2T = Tn(sg,zgrid);
Psi = Psifun(sg,T)
∂Psi∂x = ∂Psi∂x_fun(sg,T,∂T)
∂2Psi∂x2 = ∂2Psi∂x2_fun(sg,T,∂T,∂2T)

function Tn(sg::SmolyakGrid,zgrid::Array{Float64,2})
	M = maximum(sg.Binds)+1
	T = Array(Float64,sg.NumGrdPts,M,sg.D)
	∂T = Array(Float64,sg.NumGrdPts,M,sg.D)
	∂2T = Array(Float64,sg.NumGrdPts,M,sg.D)
	for n in 1:sg.NumGrdPts
		for m in 1:M
			for d in 1:sg.D
				if m==1
					T[n,m,d] = 1.0
					∂T[n,m,d] = 0.0
					∂2T[n,m,d] = 0.0				
				elseif m==2
					T[n,m,d] = zgrid[n,d]
					∂T[n,m,d] = 1.0
					∂2T[n,m,d] = 0.0				
				else
					T[n,m,d] = 2*zgrid[n,d]*T[n,m-1,d] - T[n,m-2,d]
					∂T[n,m,d] = 2*T[n,m-1,d] + 2*zgrid[n,d]*∂T[n,m-1,d] - ∂T[n,m-2,d]
					∂2T[n,m,d] = 4*∂T[n,m-1,d] + 2*zgrid[n,d]*∂2T[n,m-1,d] - ∂2T[n,m-2,d]
				end
			end
		end
	end
	return T, ∂T, ∂2T
end

function Psi_fun(sg::SmolyakGrid,T::Array{Float64,3})
	NBF = size(sg.Binds,1)
	Psi = ones(Float64,sg.NumGrdPts,NBF)
	for n in 1:sg.NumGrdPts
		for m in 1:NBF
			for d in 1:sg.D
				Psi[n,m] *= T[n,sg.Binds[m,d]+1,d]
			end
		end
	end
	return Psi
end

function ∂Psi∂x_fun(sg::SmolyakGrid,T::Array{Float64,3},∂T::Array{Float64,3})
	NBF = size(sg.Binds,1)
	∂Psi∂x = ones(Float64,sg.NumGrdPts,NBF,sg.D)
	for n in 1:sg.NumGrdPts
		for m in 1:NBF
			for d in 1:sg.D
				for i in 1:sg.D
					i==d 	?	∂Psi∂x[n,m,i] *= ∂T[n,sg.Binds[m,d]+1,d]
						 	:	∂Psi∂x[n,m,i] *= T[n,sg.Binds[m,d]+1,d]
				end
			end
		end
	end
	return ∂Psi∂x
end

function ∂2Psi∂x2_fun(sg::SmolyakGrid,T::Array{Float64,3},∂T::Array{Float64,3},∂2T::Array{Float64,3})
	NBF = size(sg.Binds,1)
	∂2Psi∂x2 = ones(Float64,sg.NumGrdPts,NBF,sg.D,sg.D)
	for n in 1:sg.NumGrdPts
		for m in 1:NBF
			for d in 1:sg.D
				for i in 1:sg.D
					for j in 1:sg.D
						if d==i==j
							∂2Psi∂x2[n,m,i,j] *= ∂2T[n,sg.Binds[m,d]+1,d]
						elseif d==j!=i  
							∂2Psi∂x2[n,m,i,j] *= ∂T[n,sg.Binds[m,d]+1,d]
						elseif d==i!=j
							∂2Psi∂x2[n,m,i,j] *= ∂T[n,sg.Binds[m,d]+1,d]																												
						else
							∂2Psi∂x2[n,m,i,j] *= T[n,sg.Binds[m,d]+1,d]
						end
					end
				end
			end
		end
	end
	return ∂2Psi∂x2
end






