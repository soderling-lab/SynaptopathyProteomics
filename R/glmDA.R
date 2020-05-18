glmDA <- function(tp,comparisons=c("Genotype","Treatment"), samples, gene_map){
	# Asses Differential Abundance with EdgeR GLM.

	# Imports.
	suppressPackageStartupMessages({
		library(data.table)
		library(dplyr)
		library(edgeR)
		library(tibble)
	})

	# Cast tp into data matrix for EdgeR.
	tp_in <- tp <- as.data.table(tp)
	dm <- tp %>% 
		dcast(Accession ~ Sample, value.var="Intensity") %>% 
		as.matrix(rownames=TRUE)

	# Create dge object.
	dge <- DGEList(counts=dm)

	# Perform TMM normalization.
	dge <- calcNormFactors(dge,method="TMM")

	# Create sample groupings given contrasts of interest.
	groups <- as.character(interaction(samples %>% select(comparisons)))
	names(groups) <- samples$Sample

	# Combine WTs?
	groups[grep("WT",groups)] <- "WT"

	# Annotate dge object with sample groups.
	dge$samples$group <- as.factor(groups[rownames(dge$samples)])

	# Annotate dge object with genetic backgrounds.
	backgrounds <- samples$Genotype
	names(backgrounds) <- samples$Sample
	group_order <- c("Shank2","Shank3","Syngap1","Ube3a")
	dge$samples$background <- factor(backgrounds[rownames(dge$samples)],
					 levels=group_order)

	# Create a design matrix for GLM.
	design <- model.matrix(~background + group, data = dge$samples)

	# Estimate dispersion.
	dge <- estimateDisp(dge, design, robust = TRUE)

	# Fit a general linear model.
	fit <- glmQLFit(dge, design, robust = TRUE)

	# Extracted fitted values.
	dm_fit <- fit$fitted.values

	# Evaluate differences for contrasts specified by design by
	# calling glmQLFTest(fit,coef=contrast)
	contrasts <- colnames(design)
	qlf <- lapply(contrasts,function(x) glmQLFTest(fit, coef=x))

	# Extract the results.
	getTopTags <- function(qlf) { topTags(qlf, n = Inf)[["table"]] }
	glm_results <- lapply(qlf, getTopTags)
	names(glm_results) <- contrasts

	# Insure first column is Accession.
	glm_results <- lapply(glm_results,function(x) {
				      Accession <- rownames(x)
				      x <- add_column(x,Accession,.before=1)
				      rownames(x) <- NULL
				      return(x)
			  })

	# Add percent WT and sort by pvalue.
	glm_results <- lapply(glm_results,function(x) {
				      x$logCPM <- 2^x$logFC
				      idy <- grep("logCPM",colnames(x))
				      colnames(x)[idy] <- "PercentWT"
				      x <- x[order(x$PValue,decreasing=FALSE),]
				      return(x)
			  })

	# Define a function to annotate results with gene ids.
	add_ids <- function(x,gene_map) {
		Uniprot <- x$Accession
		idx <- match(Uniprot,gene_map$uniprot)
		Symbol <- gene_map[["symbol"]][idx]
		Entrez <- gene_map[["entrez"]][idx]
		x <- tibble::add_column(x,Symbol,.after=1)
		x <- tibble::add_column(x,Entrez,.after=2)
		rownames(x) <- NULL
		return(x)
	}

	# Annotate with gene ids.
	glm_results <- lapply(glm_results,function(x) add_ids(x,gene_map))

	# Add fitted values to results..
	# FIXME:
	for (i in 1:length(glm_results)){
		df <- glm_results[[i]]
		contrast <- contrasts[i]
		keep <- names(which(design[,contrast] != 0))
		dm <- log2(dge$counts[,keep])
		dt <- as.data.table(dm,keep.rownames="Accession")
		glm_results[[i]] <- left_join(df,dt,by="Accession")
	}

	# Return list of normalized data and results.
	return(list(data=dm_fit,results=glm_results))
}
