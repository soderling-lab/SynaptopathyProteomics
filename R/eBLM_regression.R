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
eBLM_regression <- function(prot_dt, groups=c("Genotype","Treatment"),
			    combine.group = "Treatment",combine.var="WT",
			    batch="Genotype",fit.to="WT",
		            value.var="Intensity",
			    treatment.ignore = NULL){

	suppressPackageStartupMessages({
		library(dplyr)
		library(WGCNA)
		library(data.table)
	})

	# Prepare the expression data.
	tp_in <- as.data.table(prot_dt)
	tp_filt <- tp_in %>% filter(Treatment != treatment.ignore) %>% 
		as.data.table()
	dm <- tp_filt %>% dcast(Accession ~ Sample, value.var=value.var) %>%
		as.matrix(rownames="Accession") %>% log2()

	# The data cannot contain any missing values.
	n_missing <- sum(is.na(dm))
	if (n_missing > 0) { stop("The data cannot contain missing values.") }

	# Sample to batch mapping.
	all_batches <- tp_filt %>% dplyr::select(all_of(batch)) %>% 
		interaction() %>% as.character()
	names(all_batches) <- tp_filt$Sample

	# Sample to treatment group mapping.
	all_groups <- tp_filt %>% dplyr::select(all_of(groups)) %>% 
		interaction() %>% as.character()
	names(all_groups) <- tp_filt$Sample

	# Combine treatment groups?
	all_groups <- tp_filt %>% dplyr::select(all_of(groups)) %>% 
		interaction() %>% as.character()
	names(all_groups) <- tp_filt$Sample
	all_groups[grep(combine.var,all_groups)] <- combine.var

	# Design dt.
	design <- data.table(batch = as.factor(all_batches[colnames(dm)]),
			     treatment = as.factor(all_groups[colnames(dm)]))

	# Correct for the batch effect using empiricalBayesLM from WGCNA package.
	fit_eblm <- empiricalBayesLM(t(dm), removedCovariates = design$batch,
				     fitToSamples = design$treatment == fit.to)

	# Extract adjusted data.
	data_eblm <- 2^(t(fit_eblm$adjustedData))

	# Merge adjusted data with input data.
	tp_eblm <- as.data.table(data_eblm,keep.rownames="Accession") %>%
		reshape2::melt(id.var="Accession", variable.name="Sample",
			       value.name="Abundance")
	tp_eblm$Sample <- as.character(tp_eblm$Sample)
	tp_result <- left_join(tp_in,tp_eblm,by=c("Sample","Accession"))

	return(tp_result)
}
