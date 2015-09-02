
using Smolyak, PyPlot

D = 2
µ = 3
lb = [-2. , -3.]
ub = [2., 3.]

sg = SmolyakGrid(D,µ,lb,ub)
==(sg.D,2) ? scatter(sg.zGrid[:,1],sg.zGrid[:,2]) : nothing

NumDeriv = 2
SpOut = 1
sb = SmolyakBasis(sg,NumDeriv,SpOut)
Xpts = -2 + 4*rand(400,D);
DomCheck = 1;
sbx = SmolyakBasis(Xpts,sg,NumDeriv,SpOut)

coef = 2ones(sb.NumBasisFun);
sp = SmolyakPoly(sb,coef);
sp = SmolyakPoly(sb,coef,1);
sp = SmolyakPoly(sb,coef,0);

∂z∂x = 
