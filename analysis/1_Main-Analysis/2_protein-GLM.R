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
  library(data.table)
  library(tibble)
})

# Load project specific data and functions.
suppressWarnings({ devtools::load_all() })

# Project directories.
datadir <- file.path(root,"data") # Input/Output data.
rdatdir <- file.path(root,"rdata") # Temporary data files.

# Load the data.
data(cortex) # tidy_prot

# Cast tp into data matrix for EdgeR. Don't log!
dm <- tidy_prot %>% 
	dcast(Accession ~ Sample, value.var="Intensity") %>% 
	as.matrix(rownames=TRUE)

# Create dge object.
dge <- DGEList(counts=dm)

# Perform TMM normalization.
dge <- calcNormFactors(dge)

# Create sample groupings given contrasts of interest.
colnames(tidy_prot)
dt <- tidy_prot %>% dplyr::select(Sample,Genotype,Treatment) %>% unique() 
idx <- match(rownames(dge$samples),dt$Sample)
genotype <- dt$Genotype[idx]
treatment <- dt$Treatment[idx]
dge$samples$group <- as.character(interaction(genotype,treatment))[idx]

# Create a design matrix for GLM.
design <- model.matrix(~ 0 + group, data = dge$samples)

# Estimate dispersion.
dge <- estimateDisp(dge, design, robust = TRUE)

# Fit a general linear model.
fit <- glmQLFit(dge, design, robust = TRUE)

# Evaluate differences within genotypes.
genotypes = unique(genotype)
qlf_list <- list()
for (geno in genotypes){
	qlf <-  glmQLFTest(fit,coef=grep(geno,colnames(design)))
	qlf_list[[geno]] <- qlf
}

# Call topTags to add FDR. Gather tabularized results.
glm_results <- lapply(qlf_list, function(x) {
			      topTags(x, n = Inf, sort.by = "PValue")$table })

# Insure first column is Accession.
glm_results <- lapply(glm_results,function(x) {
			      as.data.table(x,keep.rownames="Accession") })

# Add percent WT and sort by pvalue.
glm_results <- lapply(glm_results,function(x) {
			      x[["Percent Control"]] <- 2^x$logFC
			      return(x) })

# Create gene map.
uniprot <- unique(glm_results[[1]]$Accession)
entrez <- mgi_batch_query(ids=uniprot)
symbols <- getPPIs::getIDs(entrez,from="entrez",to="symbol",species="mouse")
gene_map <- as.data.table(keep.rownames="uniprot",entrez)
gene_map$symbol <- symbols[as.character(gene_map$entrez)]

# Add gene annotations.
results_list <- list()
for (i in c(1:length(glm_results))) {
	df <- glm_results[[i]]
	idx <- match(df$Accession,gene_map$uniprot)
	Entrez <- gene_map$entrez[idx]
	Symbol <- gene_map$symbol[idx]
	df <- tibble::add_column(df,Entrez,.after=1)
	df <- tibble::add_column(df,Symbol,.after=2)
	results_list[[i]] <- df
}

names(results_list) <- gsub("\\(|\\)","",colnames(design))

# Save.
myfile <- file.path(root,"tables",paste0(tissue,"_GLM_Results.xlsx"))
write_excel(results_list,file=myfile)
