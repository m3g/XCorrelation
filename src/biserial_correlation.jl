function biserial_correlation(
	discreteVariable,
	continuousArray :: Array{Float64})

	## Sanity check
	
	discreteArray = discreteVariable .* 1.

	length(discreteArray) == length(continuousArray) || throw(DimensionMismatch("Size of discrete Array not equal to size of continuous Array"))

	## Computation

	Ntrue = sum(discreteArray .== 1)
	Nfalse = sum(discreteArray .== 0)

	SmeanTrue = Statistics.mean(continuousArray[discreteArray .== 1])
	SmeanFalse = Statistics.mean(continuousArray[discreteArray .== 0])

	sdev = Statistics.std(continuousArray, corrected=false)

	bisCorr = ((SmeanTrue - SmeanFalse)/sdev)*sqrt(Ntrue*Nfalse/(length(discreteArray)^2))

	(bisCorr == NaN) && (bisCorr = -1.)

	return( bisCorr )

end
