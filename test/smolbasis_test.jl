# TEST OF DERIVATIVES USING D = 2,  µ = 1 if using [-1,1]^2 
using Smolyak
D = 2
µ = 1
sg = SmolyakGrid(D,µ)
NumDeriv = 2
sb = SmolyakBasis(sg,NumDeriv)
makeBF!(sb)
coef = rand(sb.NumBF)
sp = SmolyakPoly(sb,coef)

#= Analytical answers =#
z1 = sg.zGrid[1,:];
z2 = sg.zGrid[2,:];

BF = [ones(1,sg.NumGrdPts), z1, 2z1.^2-1, z2, 2z2.^2-1];
dBFdz1  = [zeros(1,sg.NumGrdPts), ones(1,sg.NumGrdPts), 4z1, zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts)];
dBFdz2  = [zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), ones(1,sg.NumGrdPts), 4z2];
d2BFdz12 = [zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), 4ones(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts)];
d2BFdz22 = [zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), 4ones(1,sg.NumGrdPts)];
d2BFdz1dz2 = zeros(sg.NumGrdPts,sg.NumGrdPts);

# Tests
Passed = 1
Passed *= <=(maximum(abs2(BF - sb.BF)),1e-14)
Passed *= <=(maximum(abs2(dBFdz1 - sb.dBFdz[:,:,1])),1e-14)
Passed *= <=(maximum(abs2(dBFdz2 - sb.dBFdz[:,:,2])),1e-14)
Passed *= <=(maximum(abs2(d2BFdz12 - sb.d2BFdz2[:,:,1,1])),1e-14)
Passed *= <=(maximum(abs2(d2BFdz22 - sb.d2BFdz2[:,:,2,2])),1e-14)
Passed *= <=(maximum(abs2(d2BFdz1dz2 - sb.d2BFdz2[:,:,1,2])),1e-14)
Passed *= <=(maximum(abs2(d2BFdz1dz2 - sb.d2BFdz2[:,:,2,1])),1e-14)

==(Passed,1) ? print("Passed the tests") : print("Problem: Didn't pass test")
