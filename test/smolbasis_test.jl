# TEST OF DERIVATIVES USING D = 2,  µ = 1 if using [-1,1]^2 
using Smolyak, PyPlot
D = 2
µ = 1
sg = SmolyakGrid(D,µ)
NumDeriv = 2
sb = SmolyakBasis(sg,false, NumDeriv)

#= Analytical answers =#
z1 = sg.zGrid[1,:];
z2 = sg.zGrid[2,:];
==(sg.D,2) ? scatter(z1,z2) : nothing

Ψ = [ones(1,sg.NumGrdPts), z1, 2z1.^2-1, z2, 2z2.^2-1];
∂Ψ∂z1  = [zeros(1,sg.NumGrdPts), ones(1,sg.NumGrdPts), 4z1, zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts)];
∂Ψ∂z2  = [zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), ones(1,sg.NumGrdPts), 4z2];
∂2Ψ∂z12 = [zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), 4ones(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts)];
∂2Ψ∂z22 = [zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), zeros(1,sg.NumGrdPts), 4ones(1,sg.NumGrdPts)];
∂2Ψ∂z1∂z2 = zeros(sg.NumGrdPts,sg.NumGrdPts);

# Tests
Passed = 1
Passed *= <=(maximum(abs2(Ψ - sb.Ψ)),1e-14)
Passed *= <=(maximum(abs2(∂Ψ∂z1 - sb.∂Ψ∂z[:,:,1])),1e-14)
Passed *= <=(maximum(abs2(∂Ψ∂z2 - sb.∂Ψ∂z[:,:,2])),1e-14)
Passed *= <=(maximum(abs2(∂2Ψ∂z12 - sb.∂2Ψ∂z2[:,:,1,1])),1e-14)
Passed *= <=(maximum(abs2(∂2Ψ∂z22 - sb.∂2Ψ∂z2[:,:,2,2])),1e-14)
Passed *= <=(maximum(abs2(∂2Ψ∂z1∂z2 - sb.∂2Ψ∂z2[:,:,1,2])),1e-14)
Passed *= <=(maximum(abs2(∂2Ψ∂z1∂z2 - sb.∂2Ψ∂z2[:,:,2,1])),1e-14)

==(Passed,1) ? print("Passed the tests") : print("Problem: Didn't pass test")
