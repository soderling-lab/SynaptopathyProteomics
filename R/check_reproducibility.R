#' check_reproducibility

# Define a function that checks reproducibility of a protein.
check_reproducibility <- function(prot_list,protein,value.var="Abundance",
				  treatment.subset="WT",
				  fun="bicor",threshold=0.8) {
	# Which proteins are highly reproducible between replicates.
	# For a given protein, cast the data into a matrix: Fraction ~ Replicate.
	dt <- prot_list[[protein]] %>% 
		filter(Treatment == treatment.subset) %>%
		as.data.table() 
        dt_temp <- dcast(dt, Genotype ~ Treatment + Channel, 
		      value.var=value.var) 
	value_cols <- grep("_",colnames(dt_temp))
	dm <- dt_temp %>% dplyr::select(all_of(value_cols)) %>% 
		as.matrix(rownames.value=dt_temp$Genotype)
	# missing values can arise in dm if outlier sample was removed.
	# Calculate the pairwise correlation between samples.
	opts <- "pairwise.complete.obs"
	if (fun == "bicor") { cormat <- WGCNA::bicor(log2(dm),
						     use=opts) }
	if (fun == "pearson") { cormat <- cor(log2(dm),method="pearson",
					      use=opts) }
	if (fun == "spearman") { cormat <- cor(log2(dm),method="spearman",
					       use=opts) }
	# Ignore self- and duplicate- comparisons.
	diag(cormat) <- NA
	cormat[cormat==0] <- NA
	cormat[lower.tri(cormat)] <- NA
	# Melt into a vector of correlation values.
	x <- reshape2::melt(cormat,na.rm=TRUE,value.name="cor")[["cor"]]
	# Check if all values are > threshold.
	check <- sum(x > threshold)
	return(check)
}
