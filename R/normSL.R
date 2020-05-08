normSL <- function(tp,groupBy="Sample"){
	suppressPackageStartupMessages({
		library(dplyr)
		library(data.table)
	})
	# Create an expression that will be evaluated to dynamically 
	# group data based on user provided grouping factors. Default's
	# to Sample -- every replicate.
	tp <- ungroup(tp)
	cmd <- paste0("group_by(tp, ",paste(groupBy,collapse=", "),")")
	# Calculate column sums for grouped samples.
	data_SL <- eval(parse(text=cmd)) %>% 
		summarize(Total=sum(Intensity,na.rm=TRUE)) %>%
		group_split()
	data_SL <- tp %>% group_by(Sample) %>%
		summarize(Total=sum(Intensity,na.rm=TRUE)) %>%
		group_split()
	# Calculate normalization factors.
	data_SL <- lapply(data_SL,function(x) {
			  x$Mean <- mean(x$Total,na.rm=TRUE)
			  x$NormFactor <- x$Mean/x$Total
			  x$Norm <- x$Total*x$NormFactor
			  return(x) })
	data_SL <- do.call(rbind,data_SL)
	# Collect normalization factors as named vector.
	normFactors <- data_SL$NormFactor
	names(normFactors) <- as.character(data_SL$Sample)
	tp$Intensity <- tp$Intensity * normFactors[as.character(tp$Sample)]
	tp <- ungroup(tp)
	tp <- as.data.table(tp)
	return(tp)
}
