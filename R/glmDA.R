glmDA <- function(tp,value.var="Abundance",
		  comparisons=c("Genotype","Treatment"), 
		  model){

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
	all_groups <- tp_in %>% dplyr::select(all_of(comparisons)) %>% 
		interaction() %>% as.character()
	names(all_groups) <- tp_in$Sample

	# Combine WTs.
	all_groups[grep("WT",all_groups)] <- "WT"

	# Annotate dge object with sample groups.
	dge$samples$group <- as.factor(all_groups[rownames(dge$samples)])

	# Annotate dge object with genetic backgrounds.
	all_backgrounds <- tp_in[["Genotype"]]
	names(all_backgrounds) <- tp_in$Sample
	group_order <- c("Shank2","Shank3","Syngap1","Ube3a")
	dge$samples$background <- factor(all_backgrounds[rownames(dge$samples)],
					 levels=group_order)

	# Create a design matrix for GLM -- using blocking model with
	# genetic background as covariate.
	#design <- model.matrix(~background + group, data = dge$samples)
	cmd <- paste0("model.matrix(",model, ", data = dge$samples)")
	design <- eval(parse(text=cmd))

	# Estimate dispersion.
	dge <- estimateDisp(dge, design, robust = TRUE)

	# Fit a general linear model.
	fit <- glmQLFit(dge, design, robust = TRUE)

	# Evaluate differences for contrasts specified by design by
	# calling glmQLFTest(fit,coef=contrast)
	contrasts <- colnames(design)
	qlf <- lapply(contrasts,function(x) glmQLFTest(fit, coef=x))
	names(qlf) <- contrasts

	# Extract the results.
	getTopTags <- function(qlf) { topTags(qlf, n = Inf)[["table"]] }
	glm_results <- lapply(qlf, getTopTags)

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

	glm_results <- lapply(glm_results,as.data.table)
	return(glm_results)
}
