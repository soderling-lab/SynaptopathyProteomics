filtProt <- function(tp,controls,remove.protein.outliers, ignore = NULL,
		     rowmax=0.5,nbins=5,nSD=4, ohw=FALSE, summary=TRUE) {

	# FIXME: Split into multiple functions.

	# Store a copy of the input data.
	tp_in <- tp

	# Remove proteins that were only identified by a single peptide.
	out1 <- unique(tp$Accession[tp$Peptides==1])
	n_ohw <- length(out1)

	# Remove proteins with too many missing values.
	missing_limit <- rowmax * sum(!grepl(controls,unique(tp$Sample)))
	data_filt <- tp %>% filter(Treatment != controls) %>% 
		group_by(Accession) %>% 
		summarize(N = length(Intensity),
			  nMissing = sum(is.na(Intensity)),
			  Remove = sum(is.na(Intensity)) > missing_limit)
	out2 <- unique(data_filt$Accession[data_filt$Remove])
	n_filt_missing <- length(out2)

	# Remove proteins with any missing QC values.
	data_QC <- tp %>% filter(Treatment==controls) %>% 
		group_by(Accession,Genotype,Treatment) %>%
		summarize(N = length(Intensity),
			  nMissing = sum(is.na(Intensity)),
			  Remove=any(is.na(Intensity))) 
	out3 <- unique(data_QC$Accession[data_QC$Remove])
	n_filt_qc <- length(out3)

	# Combine proteins out, and then remove from data.
	proteins_out <- unique(c(out1,out2,out3))
	tp_filt <- tp %>% filter(Accession %notin% proteins_out)
	tp_filt <- as.data.frame(tp_filt)

	# Remove outlier protein measurements.
	if (remove.protein.outliers) { 

		# Split tidyprot dt into list of proteins.
		# Ignore any specified proteins.
		tp_list <- tp_filt %>% 
			filter(Treatment != controls, Accession %notin% ignore) %>%
			group_by(Accession,Treatment) %>% 
			group_split()
		names(tp_list) <- sapply(tp_list,function(x) unique(x$Accession))

		# Loop to calculate MEAN ratio of log2 intensities for all
		# replicate comparisons. This can take a little time.
		tp_list <- lapply(tp_list,function(x) {
					  x$Ratio <- mean(logRatios(x$Intensity))
					  return(x)
			  })

		# Bind results together as df.
		ratio_data <- do.call(rbind,tp_list)

		# Bin data by intensity.
		breaks <- quantile(ratio_data$Intensity,
				   seq(0,1,length.out=nbins+1),
				   names=FALSE, na.rm=TRUE)
		ratio_data$Bin <- cut(ratio_data$Intensity,breaks,
				      labels=FALSE,include.lowest=TRUE)

		# Summarize bins.
		ratio_df <- ratio_data %>% group_by(Bin) %>% 
			summarize("Median"= median(Intensity),
				  "Mean"= mean(Ratio,na.rm=TRUE),
				  "Std" = sd(Ratio),
				  "N" = sum(!is.na(Ratio)),
				  "Min" = mean(Ratio,na.rm=TRUE)-(nSD*sd(Ratio)),
				  "Max" = mean(Ratio,na.rm=TRUE)+(nSD*sd(Ratio)))

		# Determine if measurement is outside percision limits.
		ratio_data$Min <- ratio_df$Min[ratio_data$Bin]
		ratio_data$Max <- ratio_df$Max[ratio_data$Bin]
		out_low <- ratio_data$Ratio < ratio_data$Min 
		out_high <- ratio_data$Ratio > ratio_data$Max
		out <- out_low | out_high
		ratio_data$isOutlier <- out

		# Summarize number of protein outliers per bin.
		nOutliers <- ratio_data %>% group_by(Bin) %>%
			summarize(n=sum(isOutlier))
		ratio_df$nOutliers <- nOutliers[["n"]]

		# Collect protein outliers.
		idx <- ratio_data$isOutlier
		protein_outliers <- unique(ratio_data$Accession[idx])
		n_outliers <- length(protein_outliers)

	} else {
		# Don't remove any protein measurements.
		protein_outliers <- NULL
	}

	# Combine proteins out, and then remove from data.
	all_out <- unique(c(proteins_out,protein_outliers))
	n_out <- length(all_out)
	tp <- tp_in %>% filter(Accession %notin% all_out) %>% as.data.table

	# Status report.
	n_prot <- length(unique(tp$Accession))
	if (summary) {
	message(paste("Removed proteins identified by one peptide.....:",
		      formatC(n_ohw,big.mark=",")))
	message(paste("Removed proteins with too many missing values..:",
		      formatC(n_filt_missing,big.mark=",")))
	message(paste("Removed proteins with any missing QC values....:", 
		      formatC(n_filt_qc,big.mark=",")))
	if (remove.protein.outliers) {
		message(paste("Removed proteins with outlier measurements.....:",
			      formatC(n_outliers,big.mark=",")))
	}
	message(paste("Final number of reproducibly quantified proteins:",
		      formatC(n_prot,big.mark=",")))
	}

	# Return tidy protein.
	return(tp)
}
