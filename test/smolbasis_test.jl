# TEST OF DERIVATIVES USING D = 2,  µ = 1 if using [-1,1]^2 

using Smolyak
D = 2
mu = ones(Int64,D)
sg = SmolyakGrid(mu)
NumDeriv = 2
sb = SmolyakBasis(sg,NumDeriv)
makeBasis!(sb)

# Polynomial
coef = rand(sb.NumBF)
sp = SmolyakPoly(sb,coef)
getValue!(sp,sb)
getGrad!(sp,sb)
getHess!(sp,sb)

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
	push!(BF,slice(InputBF,i,:))
end

InputdBFdz1  = hcat(zeros(sg.NumGrdPts), ones(sg.NumGrdPts), 4z1, zeros(sg.NumGrdPts), zeros(sg.NumGrdPts));
dBFdz1 = []
for i in 1:5
	push!( dBFdz1,slice(InputdBFdz1,i,:) )
end

InputdBFdz2  = hcat(zeros(sg.NumGrdPts), zeros(sg.NumGrdPts), zeros(sg.NumGrdPts), ones(sg.NumGrdPts), 4z2);
dBFdz2 = []
for i in 1:5
	push!( dBFdz2,slice(InputdBFdz2,i,:) )
end

Inputd2BFdz11 = hcat(zeros(sg.NumGrdPts), zeros(sg.NumGrdPts), 4ones(sg.NumGrdPts), zeros(sg.NumGrdPts), zeros(sg.NumGrdPts));
d2BFdz11 = []
for i in 1:5
	push!(d2BFdz11,slice(Inputd2BFdz11,i,:) )
end

Inputd2BFdz22 = hcat(zeros(sg.NumGrdPts), zeros(sg.NumGrdPts), zeros(sg.NumGrdPts), zeros(sg.NumGrdPts), 4ones(sg.NumGrdPts));
d2BFdz22 = []
for i in 1:5
	push!(d2BFdz22,slice(Inputd2BFdz22,i,:) )
end

d2BFdz12 = [[0.0 for i in 1:5] for j in 1:5]

# Tests
Passed = 1
for i in 1:5, j in 1:5
	Passed *= isapprox(BF[i][j],sb.BF[i][j])
	Passed *= isapprox(dBFdz1[i][j],sb.dBFdz[i][1][j])
	Passed *= isapprox(dBFdz2[i][j],sb.dBFdz[i][2][j])
	Passed *= isapprox(d2BFdz11[i][j],sb.d2BFdz2[i][1][1][j])
	Passed *= isapprox(d2BFdz22[i][j],sb.d2BFdz2[i][2][1][j])
	Passed *= isapprox(d2BFdz12[i][j],sb.d2BFdz2[i][1][2][j])
end

==(Passed,1) ? print("Passed the tests") : print("Problem: Didn't pass test")
