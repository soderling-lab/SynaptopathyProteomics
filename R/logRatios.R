logRatios <- function(x){

	# Calculate the log ratio of all comparisons of two 
	# elements in vector x.

	# Declare a simple function to calculate logRatio of x and y.
	logRatio <- function(x,y) { 
		return(log2(y) - log2(x)) 
	}

	# Generate all combinations of n elements.
	choose_from <- seq(x)
	combn_dm <- combn(choose_from,m=2)
	combn_list <- split(combn_dm,rep(c(1:dim(combn_dm)[2]),each=2))
	contrasts <- lapply(combn_list,function(x) {
				    names(x) <- c("x","y")
				    return(x) 
			    })

	# For every contrast, calculate log ratio of x and y.
	ratios <- sapply(contrasts,function(idx) {
				 do.call(logRatio,as.list(x[idx]))
			    })

	return(ratios)
}
