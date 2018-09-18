# Check the correlation of different polynomials

using Statistics, Plots
sp = SpreadPoly(collect(0.01:.01:1.;), 10, [[0.,1.0] for d in 1:1]);
makeBF!(sp);
CM = Statistics.cor(VVtoMatrix(sp.BF)[2:end,2:end]);
display("****************************************");
display("Correlation of Basis Fun.: Spread Poly");
display("****************************************");
round.(CM,2)
spy_spread = spy(CM, title="Spread", ticks = nothing);
basis_spread = plot(VVtoMatrix(sp.BF), title="Spread", legend=false, xticks = nothing);

cheb = ChebyshevPoly(collect(-1.0:.02:1.0;), 10, [[-1.0,1.0] for d in 1:1]);
makeBF!(cheb);
CM = Statistics.cor(VVtoMatrix(cheb.BF)[2:end,2:end]);
display("****************************************");
display("Correlation of Basis Fun.: Chebyshev Poly");
display("****************************************");
round.(CM,2)
spy_cheb = spy(CM, title="Chebyshev", ticks = nothing);
basis_cheb = plot(VVtoMatrix(cheb.BF),  title="Chebyshev", legend=false, xticks = nothing);

op = OrdinaryPoly(collect(0.01:.01:.99;), 10);
makeBF!(op);
CM = Statistics.cor(VVtoMatrix(op.BF)[2:end-1,2:end-1]);
display("****************************************");
display("Correlation of Basis Fun.: Ordinary Poly");
display("****************************************");
round.(CM,2)
spy_op = spy(CM, title="Ordinary", ticks = nothing);
basis_op = plot(VVtoMatrix(op.BF), title="Ordinary", legend=false, xticks = nothing);

plot(spy_spread, basis_spread, spy_cheb, basis_cheb, spy_op, basis_op, layout=(3,2))