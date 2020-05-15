#!/usr/bin/env Rscript
intrabatch_combat <- function(data_in,samples,group,batch,
			      ignore="QC",r2_threshold=0.3) {
	# Subset data and cast into matrix.
	idx <- grepl(group,data_in$Sample) & !grepl(ignore,data_in$Sample)
	dm <- data_in %>% 
		filter(idx) %>% 
		as.data.table() %>%
		dcast(Accession ~ Sample,value.var="Intensity") %>%
		as.matrix(rownames="Accession")
	# Ignore rows with missing values.
	rows_out <- apply(dm, 1, function(x) any(is.na(x)))
	subdm <- dm[!rows_out, ]
	# Get Traits info.
	subtraits <- samples %>% filter(Sample %in% colnames(subdm))
	# Check: is there more than one batch?
	if (length(unique(subtraits$PrepDate)) == 1) { 
		warning(paste("For group:",group, "-",
			      "Cannot perform ComBat with one batch!"),
			call.=FALSE)
		return(data_in)
	}
	# Define batch covariate.
	batch_cov <- as.numeric(as.factor((subtraits[[batch]]))) - 1
	names(batch_cov) <- subtraits$Sample
	# Check the correlation between batch and PC1.
	pc1 <- prcomp(t(log2(subdm)))$x[, 1]
	r2_1 <- cor(batch_cov[names(pc1)], pc1)
	# Check: is there evidence of batch effect?
	if (abs(r2_1) < r2_threshold) {
		warning(paste("For group:",group, "-",
			      "No detectable batch effect!",
			      "Not applying ComBat."), call.=FALSE)
		return(data_in)
	}
	# Insure that traits data is in matching order.
	rownames(subtraits) <- subtraits$Sample
	subtraits <- subtraits[colnames(subdm),]
	# Create ComBat model.
	combat_model <- model.matrix(~ as.factor(subtraits$Treatment),
				     data = as.data.frame(log2(subdm)))
	# Apply ComBat.
	message(paste("Performing ComBat for group:",group))
	data_combat <- suppressMessages({
		sva::ComBat(dat = log2(subdm),
			    batch = as.vector(subtraits$PrepDate), 
			    mod = combat_model, mean.only = FALSE)
	})
	# Correlation between batch and PC1 post-ComBat.
	pc1 <- prcomp(t(data_combat))$x[, 1]
	r2_2 <- cor(batch_cov[names(pc1)], pc1)
	message(paste("Initial coorelation between batch and samples:",
		      round(r2_1,3)))
	message(paste("Final coorelation between batch and samples:",
		      round(r2_1,3),"\n"))
	# Put data back together again.
	data_out <- reshape2::melt(2^data_combat,value.name="Intensity")
	colnames(data_out)[c(1,2)] <- c("Accession","Sample")
	data_out$Accession <- as.character(data_out$Accession)
	data_out$Sample <- as.character(data_out$Sample)
	tidy_out <- left_join(data_in,data_out,by=c("Accession","Sample"))
	# Add regressed data.
	idx <- !is.na(tidy_out$Intensity.y)
	tidy_out$Intensity.x[idx] <- tidy_out$Intensity.y[idx]
	tidy_out$Intensity.y <- NULL
	colnames(tidy_out)[which(colnames(tidy_out)=="Intensity.x")] <- "Intensity"
	return(tidy_out)
}
