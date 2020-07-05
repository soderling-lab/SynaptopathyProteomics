#!/usr/bin/env Rscript

#' ---
#' title:
#' description: Preprocessing of TMT data.
#' authors: Tyler W Bradshaw
#' ---

## Optional Parameters:

#---------------------------------------------------------------------
## Misc function - getrd().
#---------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
	in_root <- function(h=here, dir=dpat) { 
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE))) 
		return(check)
	}
	# Loop to find root.
	while (!in_root(here)) { 
		here <- dirname(here) 
	}
	root <- here
	return(root)
}

# Parse the command line arguments.
parse_args <- function(default="Cortex", args=commandArgs(trailingOnly=TRUE)){
	# Input must be Cortex or Striatum.
	msg <- c("Please specify a tissue type to be analyzed:\n",
	 "Choose either 'Cortex' or 'Striatum'.")
	# If interactive, return default tissue.
	if (interactive()) { 
		return("Cortex") 
	} else {
		# Check arguments.
		check <- !is.na(match(args[1], c("Cortex", "Striatum")))
		if (length(args == 1) & check) { 
			tissue  <- args[1]
			start <- Sys.time()
			message(paste("Starting analysis at:", start))
			message(paste0("Analyzing ", tissue,"..."))
		} else {
			stop(msg) 
		}
		return(tissue)
	}
}

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. Load custom functions and prepare
# the project directory for saving output files.

# Parse input.
tissue <- parse_args()

# Load the R env.
root <- getrd()
renv::load(root,quiet=TRUE)

# Load required packages.
suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(tibble)
  library(data.table)
})

# Load project specific data and functions.
devtools::load_all(quiet=TRUE)

# Project directories.
datadir <- file.path(root,"data") # Input/Output data.
rdatdir <- file.path(root,"rdata") # Temporary data files.

# Load the data.
data(cortex) # tidy_prot
#data(striatum)

data(gene_map)

#--------------------------------------------------------------------
## EdgeR glm
#--------------------------------------------------------------------

# Cast tp into data matrix for EdgeR. Don't log!
groups <- unique(tidy_prot$Genotype)

# Loop:
results_list <- list()
for (geno in groups){
	dm <- tidy_prot %>% filter(Genotype==geno) %>%
		dcast(Accession ~ Sample, value.var="Intensity") %>% 
		as.matrix(rownames=TRUE)
	# Create dge object.
	dge <- DGEList(counts=dm)
	# Perform TMM normalization.
	dge <- calcNormFactors(dge)
	# Sample mapping.
	samples <- rownames(dge$samples)
	idx <- match(samples,tidy_prot$Sample)
	genotype <- tidy_prot$Genotype[idx]
	treatment <- tidy_prot$Treatment[idx]
	dge$samples$group <- interaction(genotype,treatment)
	# Basic design matrix for GLM -- all groups treated seperately.
	design <- model.matrix(~ 0 + group, data = dge$samples)
	# Estimate dispersion:
	dge <- estimateDisp(dge, design, robust = TRUE)
	# Fit a general linear model.
	fit <- glmQLFit(dge, design, robust = TRUE)
	# Create contrast.
	wt <- colnames(design)[grepl("WT",colnames(design))]
	mut <- colnames(design)[grepl("HET|KO",colnames(design))]
	contr <- limma::makeContrasts(paste(mut,wt,sep="-"),levels=design)
	# Assess differences.
	qlf <- glmQLFTest(fit,contrast=contr)
	# Call topTags to add FDR. Gather tabularized results.
	glm_results <- topTags(qlf, n = Inf, sort.by = "PValue")$table
	# Use gene map to annotate glm_results with entrez Ids and gene symbols.
	idx <- match(rownames(glm_results), gene_map$uniprot)
	glm_results <- add_column(glm_results, "Accession"=rownames(glm_results),
				  .before = 1)
	glm_results <- add_column(glm_results, "Entrez" = gene_map$entrez[idx], 
				  .after = 1)
	glm_results <- add_column(glm_results, "Symbol" = gene_map$symbol[idx], 
				  .after = 2)
	# Add expression data.
	wt_cols <- grep(strsplit(wt,"\\.")[[1]][2],colnames(dm))
	mut_cols <- grep(strsplit(mut,"\\.")[[1]][2],colnames(dm))
	dt <- as.data.table(dm[,c(wt_cols,mut_cols)],keep.rownames="Accession")
	glm_results <- left_join(glm_results,dt,by="Accession")
	results_list[[geno]] <- glm_results
}

