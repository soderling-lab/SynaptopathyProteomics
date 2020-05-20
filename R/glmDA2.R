glmDA <- function(tp,treatment.ignore = "QC", value.var="Intensity",
		  comparisons=c("Genotype","Treatment"), 
		  combine.WT=FALSE){

	# Assess Differential Abundance with EdgeR GLM.
	# Imports.
	suppressPackageStartupMessages({
		library(data.table)
		library(dplyr)
		library(edgeR)
		library(tibble)
	})

	# Perform the analysis WITHOUT combining WT groups.
	if (!combine.WT) {

		# Cast tp into data matrix for EdgeR.
		tp_in <- tp %>% as.data.table() %>% 
			filter(Treatment != treatment.ignore) %>% as.data.table()
		dm <- tp_in %>% 
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

		# Annotate dge object with sample groups.
		dge$samples$group <- as.factor(all_groups[rownames(dge$samples)])

		# Create a design matrix for GLM. Each group is seperate.
		design <- model.matrix(~ 0 + group, data = dge$samples)

		# Estimate dispersion.
		dge <- estimateDisp(dge, design, robust = TRUE)

		# Create all pairwise contrasts.
		# NOTE: makeContrasts doesn't work when passed env variables 
		# when not run interactively.
contrasts <- list(
makeContrasts('groupCortex.Shank2.KO - groupCortex.Shank2.WT',levels=design),
makeContrasts('groupCortex.Shank3.KO - groupCortex.Shank3.WT',levels=design),
makeContrasts('groupCortex.Syngap1.HET - groupCortex.Syngap1.WT',levels=design),
makeContrasts('groupCortex.Ube3a.KO - groupCortex.Ube3a.WT',levels=design),
makeContrasts('groupStriatum.Shank2.KO - groupStriatum.Shank2.WT',levels=design),
makeContrasts('groupStriatum.Shank3.KO - groupStriatum.Shank3.WT',levels=design),
makeContrasts('groupStriatum.Syngap1.HET - groupStriatum.Syngap1.WT',levels=design),
makeContrasts('groupStriatum.Ube3a.KO - groupStriatum.Ube3a.WT',levels=design))

		# Fit a general linear model.
		fit <- glmQLFit(dge, design, robust = TRUE)

		# Get observed values.
		new.names <- c("Accession","Sample","Obs.Intensity")
		suppressWarnings({
			obs_dt <- cpm(dge, log=FALSE, 
				      prior.count=fit$prior.count) %>%
			reshape2::melt() %>% 
			setNames(nm=new.names) %>%
			left_join(tp_in,by=c("Accession","Sample")) %>%
			as.data.table()
		})
		new.names <- c("Accession","Sample","Fit.Intensity")
		suppressWarnings({
			fit_dt <- cpm(fit, log=FALSE) %>%
				reshape2::melt() %>% 
				setNames(nm=new.names) %>%
				left_join(tp_in,by=c("Accession","Sample")) %>%
				as.data.table()
		})

		# Evaluate differences for all specified contrasts.
		qlf <- lapply(contrasts,function(x) glmQLFTest(fit,contrast=x))

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
					      x <- x[order(x$PValue,
							   decreasing=FALSE),]
					      return(x)
				  })
		glm_results <- lapply(glm_results,as.data.table)

		# Annotate with gene IDs.
		annotate_results <- function(df){
			idx <- match(df$Accession,tp_in$Accession)
			Symbol <- tp_in$Symbol[idx]
			Entrez <- tp_in$Entrez[idx]
			df <- tibble::add_column(df,Symbol,.after="Accession")
			df <- tibble::add_column(df,Entrez,.after="Symbol")
			return(df)
		}
		glm_results <- lapply(glm_results,annotate_results)

	} else if (combine.WT) {

	# Do analysis WITH combined WT samples.
	combine.var = "WT"

	# Cast tp into data matrix for EdgeR.
	tp_in <- tp %>% as.data.table() %>% 
		filter(Treatment != treatment.ignore) %>% as.data.table()
	dm <- tp_in %>% 
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
	message(paste("Combining",combine.var,"group!"))
	all_groups[grep(combine.var,all_groups)] <- combine.var

	# Annotate dge object with sample groups.
	dge$samples$group <- as.factor(all_groups[rownames(dge$samples)])

	# Annotate dge object with tissue groups.
	all_tissues <- tp_in[["Tissue"]]
	names(all_tissues) <- tp_in$Sample
	dge$samples$tissue <- factor(all_tissues[rownames(dge$samples)],
					 levels=c("Cortex","Striatum"))

	# Annotate dge object with genetic backgrounds.
	all_backgrounds <- tp_in[["Genotype"]]
	names(all_backgrounds) <- tp_in$Sample
	group_order <- c("Shank2","Shank3","Syngap1","Ube3a")
	dge$samples$background <- factor(all_backgrounds[rownames(dge$samples)],
					 levels=group_order)

	# Create a design matrix for GLM -- using blocking model with
	# genetic background as covariate.
	design <- model.matrix(~background + tissue + group, data = dge$samples)
	#colnames(design)[1] <- "Intercept"

	# Estimate dispersion.
	dge <- estimateDisp(dge, design, robust = TRUE)

	# Fit a general linear model.
	fit <- glmQLFit(dge, design, robust = TRUE)

	# Get observed data.
	new.names <- c("Accession","Sample","Obs.Intensity")
	suppressWarnings({
		obs_dt <- cpm(dge, log=FALSE, prior.count=fit$prior.count) %>%
			reshape2::melt() %>% 
			setNames(nm=new.names) %>%
			left_join(tp_in,by=c("Accession","Sample")) %>%
			as.data.table()
	})
	# Get fitted data.
	new.names <- c("Accession","Sample","Fit.Intensity")
	suppressWarnings({
	fit_dt <- cpm(fit, log=FALSE) %>%
		reshape2::melt() %>% 
		setNames(nm=new.names) %>%
		left_join(tp_in,by=c("Accession","Sample")) %>%
		as.data.table()
	})

	# Evaluate differences for contrasts specified by design by
	# calling glmQLFTest(fit,coef=contrast)
	contrasts <- colnames(design)
	qlf <- lapply(contrasts,function(x) glmQLFTest(fit,coef=x))

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

	# Annotate with gene IDs.
	annotate_results <- function(df){
		idx <- match(df$Accession,tp_in$Accession)
		Symbol <- tp_in$Symbol[idx]
		Entrez <- tp_in$Entrez[idx]
		df <- tibble::add_column(df,Symbol,.after="Accession")
		df <- tibble::add_column(df,Entrez,.after="Symbol")
		return(df)
	}
	glm_results <- lapply(glm_results,annotate_results)

	# Final clean-up.
	glm_results <- lapply(glm_results,as.data.table)
	names(glm_results) <- contrasts
	#glm_results <- glm_results[c(5,6,7,8)]
	#names(glm_results) <- c("Shank3","Syngap1","Ube3a","Shank2")
	#glm_results <- glm_results[c("Shank2","Shank3","Syngap1","Ube3a")]

	}

	return(list(results=glm_results,obs.vals=obs_dt,fit.vals=fit_dt))
}
