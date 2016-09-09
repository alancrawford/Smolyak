# TEST OF DERIVATIVES USING D = 2,  µ = 1 if using [-1,1]^2 

using Smolyak, FactCheck
D = 2
mu = ones(Int64,D)
sg = SmolyakGrid(mu)
sb = SmolyakBasis(sg;NumDeriv=2)
makeBasis!(sb)

# Polynomial
coef = rand(sb.NumBF)
sp = SmolyakPoly(sb;Coef=coef)
makeValue!(sp,sb)
makeGrad!(sp,sb)
makeHess!(sp,sb)

# ---------------------------------------------------------------------------
#= Analytical answers =#

z1 = Float64[]
z2 = Float64[]
for i in 1:sg.NumGrdPts
	push!(z1,sg.zGrid[i][1])
	push!(z2,sg.zGrid[i][2])
end

InputBF = hcat(ones(sg.NumGrdPts), z1, 2z1.^2-1, z2, 2z2.^2-1);
BF = []
for i in 1:5
	push!(BF,view(InputBF,i,:))
end

InputdBFdz1  = hcat(zeros(sg.NumGrdPts), ones(sg.NumGrdPts), 4z1, zeros(sg.NumGrdPts), zeros(sg.NumGrdPts));
dBFdz1 = []
for i in 1:5
	push!( dBFdz1,view(InputdBFdz1,i,:) )
end

InputdBFdz2  = hcat(zeros(sg.NumGrdPts), zeros(sg.NumGrdPts), zeros(sg.NumGrdPts), ones(sg.NumGrdPts), 4z2);
dBFdz2 = []
for i in 1:5
	push!( dBFdz2,view(InputdBFdz2,i,:) )
end

Inputd2BFdz11 = hcat(zeros(sg.NumGrdPts), zeros(sg.NumGrdPts), 4ones(sg.NumGrdPts), zeros(sg.NumGrdPts), zeros(sg.NumGrdPts));
d2BFdz11 = []
for i in 1:5
	push!(d2BFdz11,view(Inputd2BFdz11,i,:) )
end

Inputd2BFdz22 = hcat(zeros(sg.NumGrdPts), zeros(sg.NumGrdPts), zeros(sg.NumGrdPts), zeros(sg.NumGrdPts), 4ones(sg.NumGrdPts));
d2BFdz22 = []
for i in 1:5
	push!(d2BFdz22,view(Inputd2BFdz22,i,:) )
end

d2BFdz12 = [[0.0 for i in 1:5] for j in 1:5]

# ---------------------------------------------------------------------------

# Tests
facts("Does Smolyak Package Calculate Correct Values for:") do
	context("Basis functions") do
		for i in 1:sg.NumGrdPts, j in 1:sb.NumBF
			@fact isapprox(BF[i][j],sb.BF[i][j]) --> true
		end
	end
	context("1st derivatives") do
		for i in 1:sg.NumGrdPts, j in 1:sb.NumBF
			@fact isapprox(dBFdz1[i][j],sb.dBFdz[i][1][j]) --> true
			@fact isapprox(dBFdz2[i][j],sb.dBFdz[i][2][j]) --> true
		end
	end
	context("2nd derivatives") do 
		for i in 1:sg.NumGrdPts, j in 1:sb.NumBF
			@fact isapprox(d2BFdz11[i][j],sb.d2BFdz2[i][1][1][j]) --> true
			@fact isapprox(d2BFdz22[i][j],sb.d2BFdz2[i][2][1][j]) --> true
			@fact isapprox(d2BFdz12[i][j],sb.d2BFdz2[i][1][2][j]) --> true
		end
	end
end
