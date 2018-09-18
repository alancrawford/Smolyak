# High-Dim Smolyak Utilities

# This is the GridIndex to create Binds
function makeHDSmolIdx(mu::Vector{Int64})
	mu_max = maximum(mu)
	D = length(mu)
	parts = Vector{Int64}[]
	for d in D:D+mu_max
		push!(parts,collect(partitions(d,D))...)
	end

	Out = ones(Int64,D)'
	for i in eachindex(parts)
		N = findfirst(parts[i],1)-1
		X = unique(collect(permutations(parts[i][1:N])))
		if N==0
		 	continue
		elseif N==1
			A = ones(Int64,D,D)
			A[diagind(A,0)] = X[1][]*ones(Int64,D)
			append!(Out, A) 
		elseif N==2
			for n in eachindex(X), i in eachindex(X[n]), j in 1:D 
				A = ones(Int64,D,D)
				A[diagind(A,0)] = X[n][1]*ones(Int64,D)
				A[:,j] = X[n][2]*ones(Int64,D)
				append!(Out, A) 	
			end
		elseif N==3
			for n in eachindex(X), i in eachindex(X[n]), j_1 in 1:D, j_2 in j_1:D
				A = ones(Int64,D,D)
				A[diagind(A,0)] = X[n][1]*ones(Int64,D)
				A[:,j_1] = X[n][2]*ones(Int64,D)
				A[:,j_2] = X[n][3]*ones(Int64,D)
				append!(Out, A)  	
			end
		else 
			throw("WARNING: No code for mu>3 for high dimensions (ie. D>20)") 
		end
	end

	GI = unique(Out[find(sum(Out,2).<=D+mu_max),:],1)
	anstrphc_indx = ones(Int64,size(GI,1))
	for i in CartesianRange(size(GI))
		if anstrphc_indx[i.I[1]]==0
			continue
		else 
		>(GI[i]-1,mu[i.I[2]]) ? anstrphc_indx[i.I[1]] = 0 : nothing
		end
	end
	GridIdx = GI[findin(anstrphc_indx,1),:]
	Res = Vector{Int64}[]
	for i in 1:size(GridIdx,1)
		push!(Res, view(GridIdx,i,:))
	end	

	return Res
end

@generated function makeNumGridPts(GridIdx::VV{Int64},mu::NTuple{N,Int64}) where N
	quote
		max_mu=0
		for i in 1:$N
			max_mu = max(max_mu,mu[i])
		end
		A = A_pidx(max_mu+1)
		NGP = 0
		for i = 1:length(GridIdx)
			NGP += length(Iterators.product((@ntuple $N k->A[GridIdx[i][k]])...))
		end
		return NGP
	end
end



