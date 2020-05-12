#' eBLM_regression
#'
#' A function for performing moderated empirical bayes regression for removal
#' of covariate effects.
#'
#' @param
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @import
#'
#' @export
#'
#' @examples
#' eBLM_regression(data_in, traits)
eBLM_regression <- function(tp, traits, ignore = NULL) {
	# Prepare the expression data.
	tp_in <- as.data.table(tp)
	tp <- tp %>% filter(Treatment != ignore) %>% as.data.table()
	dt <- tp %>% dcast(Accession ~ Sample, value.var="Intensity")
	dm <- as.matrix(dt,rownames="Accession") %>% log2() %>% t()
	# The data cannot contain any missing values.
	n_missing <- sum(is.na(dm))
	if (n_missing > 0) { stop("The data cannot contain missing values.") }
	# Prepare the design df.
	traits_sub <- traits %>% filter(Sample %in% rownames(dm))
	sex <- as.factor(traits_sub$Sex)
	age <- as.numeric(traits_sub$Age)
	batch <- as.factor(traits_sub$PrepDate)
	status <- traits_sub$Treatment
	design <- as.data.frame(cbind(status, batch, sex, age))
	# Define covariates to be removed. Which model to choose?
	covariates <- cbind(design$batch, design$sex) # Batch + Sex
	# Correct for the batch effect using empiricalBayesLM from WGCNA package.
	fit_eblm <- empiricalBayesLM(dm, removedCovariates = covariates,
				     fitToSamples = design$status == "WT")
	data_eblm <- fit.eblm$adjustedData
	# Return regressed data, un-logged.
	data_out <- 2^(t(data_eblm))
	tp_eblm <- as.data.table(data_out,keep.rownames="Accession") %>%
		reshape2::melt(id.var="Accession", variable.name="Sample",
			       value.name="Intensity")
	tp_eblm$Sample <- as.character(tp_eblm$Sample)
	# Merge regressed data with input data.
	tp_result <- full_join(tp_in,tp_eblm,by=c("Sample","Accession"))
	# Replace missing values.
	idx <- is.na(tp_result$Intensity.y)
	tp_result$Intensity.y[idx] <- tp_result$Intensity.x[idx]
	tp_result$Intensity.x <- NULL
	# Fix column names.
	idy <- which(colnames(tp_result) == "Intensity.y")
	colnames(tp_result)[idy] <- "Intensity"
	return(tp_result)
}
