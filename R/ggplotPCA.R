#' ggplotPCA

ggplotPCA <- function(tp,value.var="Abundance",groups=c("Genotype","Treatment"),
		      log=FALSE, combine.group = NULL, combine.var = NULL, 
		      treatment.ignore = "QC", group.order=NULL) {

	suppressPackageStartupMessages({
		library(ggplot2)
	})

	# Cast data into a matrix.
	dm <- tp %>% filter(Treatment != treatment.ignore) %>%
		as.data.table() %>%
		dcast(Accession ~ Sample,value.var=value.var) %>%
		as.matrix(rownames=TRUE)

	# Create sample to group mapping.
	all_groups <- tp %>% dplyr::select(all_of(groups)) %>% 
		interaction() %>% as.character()
	names(all_groups) <- tp$Sample

	# Combine groups if specified.
	if (!is.null(combine.var) & !is.null(combine.group)) {
		idx <- grepl(combine.var,tp[[combine.group]])
		all_groups[idx] <- combine.var
	}

	# Check for missing values.
	missing_vals <- any(is.na(dm))
	if (missing_vals) { stop("Data contains missing values!") }

	# Perform PCA.
	if (log) {
		pca <- prcomp(t(log2(dm)))[["x"]][,1:2]
	} else {
		pca <- prcomp(t(dm))[["x"]][,1:2]
	}

	# Data for plotting.
	pca_dt <- as.data.table(pca,keep.rownames="Sample")
	pca_dt$group <- as.factor(all_groups[pca_dt$Sample])

	# Generate plot.
	plot <- ggplot(pca_dt, aes(x=PC1,y=PC2)) + 
		geom_text(aes(label=pca_dt$group))

	return(plot)
}
