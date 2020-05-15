glmDA <- function(tp,comparisons, samples, gene_map,
		  samples_to_ignore, alpha=0.1){
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
	dge <- calcNormFactors(dge)

	# Extract normalized data from dge object.
	dm_dge <- dge$counts %>% as.data.table(keep.rownames="Accession")
	tp_dge <- melt(dm_dge,id.vars="Accession",
		       variable.name="Sample",value.name="Intensity")
	tp_dge$Sample <- as.character(tp_dge$Sample)

	# Merge with meta data from input.
	tp$Intensity <- NULL
	tp_dge <- left_join(tp,tp_dge,by=c("Sample","Accession")) %>% 
		as.data.table()

	# Drop samples to ignore from dge object.
	# We dont want to utilize QC when estimating dispersion or 
	# performing statistical testing.
	drop <- grepl(samples_to_ignore,rownames(dge$samples))
	dge$samples <- dge$samples[!drop,]
	drop <- grepl(samples_to_ignore,colnames(dge$counts))
	dge$counts <- dge$counts[,!drop]

	# Create sample groupings given contrasts of interest.
	subsamples <- samples %>% filter(Sample %in% colnames(dge$counts))
	contrasts <- unlist(strsplit(comparisons,"\\."))
	subsamples$Group <- apply(subsamples[,contrasts],1,paste,collapse=".")
	groups <- subsamples$Group
	names(groups) <- subsamples$Sample
	groups <- as.factor(groups[rownames(dge$samples)])
	dge$samples$group <- groups

	# Create a design matrix for GLM.
	design <- model.matrix(~ 0 + groups, data = dge$samples)
	colnames(design) <- levels(dge$samples$group)

	# Estimate dispersion.
	dge <- estimateDisp(dge, design, robust = TRUE)

	# Fit a general linear model.
	fit <- glmQLFit(dge, design, robust = TRUE)

	# Generate pairwise contrasts list.
	g <- as.character(groups)
	idx <- sapply(strsplit(g,"\\."),"[",1)
	contrasts_list <- split(g,idx)
	contrasts_list <- lapply(contrasts_list,unique)
	contrasts_list <- lapply(contrasts_list,function(x) x[order(x)])
	contrasts_list <- lapply(contrasts_list,function(x) { 
					 paste(x,collapse="-") })

	# Loop to generate contrasts for edgeR::glm.
	for (i in 1:length(contrasts_list)) {
		contrast <- contrasts_list[[i]]
		cmd <- paste0("makeContrasts(",contrast,", levels = design)")
		contrasts_list[[i]] <- eval(parse(text=cmd))
	}

	# Call glmQLFTest() to evaluate differences in contrasts.
	qlf <- lapply(contrasts_list,function(x) glmQLFTest(fit, contrast = x))

	# Call topTags to add FDR. Gather tabularized results.
	glm_results <- lapply(qlf, function(x) {
				      topTags(x, n = Inf, sort.by = "none")$table
			  })

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

	# Add expression data to results.
	for (i in 1:length(glm_results)){
		df <- glm_results[[i]]
		contrast <- contrasts_list[[i]]
		idx <- dge$samples$group %in% names(which(contrast[,1] != 0))
		keep <- rownames(dge$samples)[idx]
		dm <- log2(dge$counts[,keep])
		dt <- as.data.table(dm,keep.rownames="Accession")
		glm_results[[i]] <- left_join(df,dt,by="Accession")
	}

	# Return list of normalized data and results.
	return(list(data=tp_dge,results=glm_results))
}
