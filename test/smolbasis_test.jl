# TEST OF DERIVATIVES USING D = 2,  µ = 1 if using [-1,1]^2 
using Smolyak	
D = 2
µ = 1
lb = -2ones(D)
ub = 2ones(D)
sg = SmolyakGrid(D,µ,lb,ub)
NumDeriv = 2
sb = SmolyakBasis(sg,NumDeriv)
==(sg.D,2) ? scatter(sg.zGrid[:,1],sg.zGrid[:,2]) : nothing

#= Analytical answers =#
x = sg.zGrid[:,1];
y = sg.zGrid[:,2];
psi = [ones(sg.NumGrdPts) x 2x.^2-1 y 2y.^2-1];
∂psi∂x  = [zeros(sg.NumGrdPts) ones(sg.NumGrdPts) 4x zeros(sg.NumGrdPts) zeros(sg.NumGrdPts)];
∂psi∂y  = [zeros(sg.NumGrdPts) zeros(sg.NumGrdPts) zeros(sg.NumGrdPts) ones(sg.NumGrdPts) 4y];
∂2psi∂x2 = [zeros(sg.NumGrdPts) zeros(sg.NumGrdPts) 4ones(sg.NumGrdPts) zeros(sg.NumGrdPts) zeros(sg.NumGrdPts)];
∂2psi∂y2 = [zeros(sg.NumGrdPts) zeros(sg.NumGrdPts) zeros(sg.NumGrdPts) zeros(sg.NumGrdPts) 4ones(sg.NumGrdPts)];
∂2psi∂x∂y = zeros(sg.NumGrdPts,sg.NumGrdPts);

# Tests
Passed = 1
Passed *= <=(maximum(abs2(psi - sb.Psi)),1e-14)
Passed *= <=(maximum(abs2(∂psi∂x - sb.∂Psi∂x[:,:,1])),1e-14)
Passed *= <=(maximum(abs2(∂psi∂y - sb.∂Psi∂x[:,:,2])),1e-14)
Passed *= <=(maximum(abs2(∂2psi∂x2 - sb.∂2Psi∂x2[:,:,1,1])),1e-14)
Passed *= <=(maximum(abs2(∂2psi∂y2 - sb.∂2Psi∂x2[:,:,2,2])),1e-14)
Passed *= <=(maximum(abs2(∂2psi∂x∂y - sb.∂2Psi∂x2[:,:,1,2])),1e-14)
Passed *= <=(maximum(abs2(∂2psi∂x∂y - sb.∂2Psi∂x2[:,:,2,1])),1e-14)

==(Passed,1) ? print("Passed the tests") : print("Problem: Didn't pass test")
