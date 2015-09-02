	for m in 1:NBF
		for n in 1:D
			@show sg.Binds[m,n]
			#dB[m,n] = -sg.Binds[m,n]*sin(sg.Binds[m,n]*sg.thetaGrid[m,n])
		end
	end


for i in 1:sg.NumGrdPts

end

(-(-(2x) * (0.5 / sqrt(1 - x^2))) * n) / sqrt(1 - x^2)^2) * sin(n * acos(x))
		 + (n / sqrt(1 - x^2)) * ((n * (-1 / sqrt(1 - x^2))) * cos(n * acos(x))


n = sg.Binds[m,n]
theta = sg.thetaGrid[m,n]
sin_theta = sin(theta)
sin_ntheta = sin(n*theta)
cos_theta = cos(theta)
cos_ntheta = cos(n*theta)

x = 0.5;
theta = acos(x);
R = {};
dR = {};
d2R = {};
for n in 0:4
	push!(R,cos(n*theta))
	push!(dR,n*sin(n*theta)./sin(theta))
	push!(d2R, (dR[n+1]*cos(theta) - n*n*R[n+1])/(sin(theta)^2))
end
[R dR d2R]

x = 0.5
T = {}
dT = {}
d2T = {}
push!(T,1.0)
push!(dT,0.0)
push!(d2T,0.0)
push!(T,x)
push!(dT,1.0)
push!(d2T,0.0)
for n in 2:4
	@show n 
	push!(T,2x*T[n] - T[n-1])
	push!(dT,2*T[n] + 2x*dT[n] - dT[n-1])
	push!(d2T,4*dT[n] + 2x*d2T[n] - d2T[n-1])
end
[T dT d2T]



d2R[n+1] 	=	 -n*n*cos(n*Z) / sin(Z)^2 + (n * sin(n*Z)/ sin(Z)) * cos(Z))) / sin(Z)^2 
			= ( dR[n+1]*cos(Z)-(n^2)*R[n+1] ) / sin(Z)^2

T = Array(Float64,NBF,D,NGP)
for i in 1:NGP
	for m in 1:NBF
		for n in 1:D
			T[m,n,i] = cos((m-1)*Grid[i,n])
		end
	end
end

Psi = Array(Float64,NBF,NGP)
for i in 1:NGP
	Psi[i,:] = prod(T[:,:,i],2)
end
round(Psi,10)
#= Checked ! Yes =#

dT = Array(Float64,NBF,D,NGP)
for i in 1:NGP
	for m in 1:NBF
		for n in 1:D
			>( 1.0 - abs(cos(Grid[i,n])) , 1e-14) ?
				dT[m,n,i] = (m-1)*sin((m-1)*Grid[i,n])/sin(Grid[i,n]) :
				dT[m,n,i] = (m-1)*(m-1) * (cos(Grid[i,n])^m)
		end
	end
end

∂Psi∂x = Array(Float64,NBF,NGP,D)
for i in 1:NGP
	for n in 1:D
		B = copy(T)
		B[:,n] = dT[:,n]
		∂Psi∂x[:,i,n] = prod(B,2)[:,:,i]'
	end
end
round(∂Psi∂x[:,:,1]',10)
round(∂Psi∂x[:,:,2]',10)



