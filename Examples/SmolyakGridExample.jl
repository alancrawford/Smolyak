using Smolyak, Plots
gr()


# Isotrophic Smolyak Grid in 2-dimensions

D = 2;
plotbin = [];
for level in 1:4
	# Create a Smolyak Kernel 
	sk = SmolyakKernel(D, level)

	# Create Smolyak grid corresponding to Smolyak Kernel above
	sg = SmolyakGrid(sk)

	# Convert to a matrix (optional - useful for plotting)
	grid = VVtoMatrix(sg.grid)

	push!(plotbin, scatter(view(grid,:,1), view(grid,:, 2), 
							title="μ=$(level)", 
							legend=false, 
							)
		)
end 
plot(plotbin... , layout = (2, 2))

# savefig("../Examples/IsotrophicSmolyakGridExample.png")

# Anisotrophic Smolyak Grid in 2-dimensions

plotbin = [];
for level in 1:4
	# Create a Smolyak Kernel 
	sk = SmolyakKernel([1,level])

	# Create Smolyak grid corresponding to Smolyak Kernel above
	sg = SmolyakGrid(sk)

	# Convert to a matrix (optional - useful for plotting)
	grid = VVtoMatrix(sg.grid)

	push!(plotbin, scatter(view(grid,:,1), view(grid,:, 2), 
							title="μ=[2, $(level)]", 
							legend=false, 
							)
		)
end 
plot(plotbin... , layout = (2, 2))

# savefig("../Examples/AnisotrophicSmolyakGridExample.png")