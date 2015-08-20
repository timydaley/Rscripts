LinearInterpolation <- function(lower.x, upper.x, inter.x, lower.y, upper.y){
	return(lower.y + (upper.y - lower.y)*(inter.x - lower.x)/(upper.x - lower.x))
}

ComplexityCurveLinearInterpolation <- function(input.x, input.y, reference.x){
	stopifnot(length(input.x) == length(input.y))
	stopifnot(tail(input.x, 1) > tail(reference.x, 1))
	stopifnot(input.x[1] == 0)
	stopifnot(reference.x[1] == 0)
	

	output.y = c(0, lapply(reference.x[2:length(reference.x)],
	                       function(x) 
	                       LinearInterpolation(input.x[tail(which(input.x < x), 1)],
	                                           input.x[head(which(input.x >= x), 1)],
	                                           x, 
	                                           input.y[tail(which(input.x < x), 1)],
	                                           input.y[head(which(input.x >= x), 1)])))
	                                          
    return(as.vector(output.y, mode = "numeric"))	                                          
}