
function VVtoMatrix(A::VV{T}) where T<:Real
	N = length(A)
	M = length(A[1])
	Mat = []
	for a in A
		append!(Mat, a)
	end
	return reshape(Mat, M, N)'
end



