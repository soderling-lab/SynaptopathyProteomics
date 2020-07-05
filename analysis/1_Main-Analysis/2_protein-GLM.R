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
#for (geno in groups){
geno = "Shank2"
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

# Assess differences.
qlf <- glmQLFTest(fit)

## Determine number of significant results with decideTests().
summary_table <- summary(decideTests(qlf))

# Call topTags to add FDR. Gather tabularized results.
qlf <- list(qlf)
glm_results <- lapply(qlf, function(x){ 
			      topTags(x, n = Inf, sort.by = "none")$table })

# Convert logCPM column to percent WT and annotate with candidate column.
glm_results <- lapply(glm_results, function(x) annotateTopTags(x))

# Use gene map to annotate glm_results with entrez Ids and gene symbols.
for (i in 1:length(glm_results)) {
  x <- glm_results[[i]]
  idx <- match(rownames(x), gene_map$ids)
  x <- tibble::add_column(x, "Gene|Uniprot" = gene_map$ids[idx], .before = 1)
  x <- tibble::add_column(x, "Uniprot" = gene_map$uniprot[idx], .after = 1)
  x <- tibble::add_column(x, "Entrez" = gene_map$entrez[idx], .after = 2)
  x <- tibble::add_column(x, "Symbol" = gene_map$gene[idx], .after = 3)
  glm_results[[i]] <- x
}

# Add expression data.
for (i in 1:length(glm_results)) {
  namen <- names(glm_results)[i]
  df <- glm_results[[i]]
  comparison <- contrasts[[namen]]
  groups <- rownames(comparison)[!comparison == 0]
  samples <- traits$SampleID[traits$group %in% groups]
  dat <- cleanDat[, samples]
  colnames(dat) <- traits$ColumnName[match(colnames(dat), traits$SampleID)]
  dat <- dat[, c(grep("HET|KO", colnames(dat)), grep("WT", colnames(dat)))]
  out <- merge(df, log2(dat), by = "row.names")
  glm_results[[i]] <- out
}

# Sort by pvalue.
glm_results <- lapply(glm_results, function(x) x[order(x$PValue), ])

# Reorder by genotype.
idx <- c(
  grep("Shank2", names(glm_results)),
  grep("Shank3", names(glm_results)),
  grep("Syngap1", names(glm_results)),
  grep("Ube3a", names(glm_results))
)
glm_results <- glm_results[idx]

# Final renaming.
namen <- unlist({
  lapply(
    lapply(strsplit(gsub("HET.|KO.", "", names(glm_results)), "\\."), rev),
    function(x) paste(x, collapse = " ")
  )
})
names(glm_results) <- namen

# Remove Row.names column.
f <- function(x) {
  x$Row.names <- NULL
  return(x)
}
glm_results <- lapply(glm_results, f)


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

#---------------------------------------------------------------------
## EdgeR statistical comparisons post-TAMPOR.
#---------------------------------------------------------------------

message("\nPerforming statistical testing with EdgeR glm...")

# Statistical comparisons are KO/HET versus all WT of a tissue type.

# Prepare data for EdgeR.
# Data should NOT be log2 transformed.
# Remove QC samples prior to passing data to EdgeR.
out <- alltraits$SampleType[match(colnames(cleanDat), 
				  rownames(alltraits))] == "QC"
data_in <- cleanDat[, !out]

# Number of proteins...
nprots <- formatC(dim(data_in)[1],big.mark=",")
nsamples <- dim(data_in)[2]
message(paste("\nQuantified", nprots, "proteins from", 
	      nsamples, "samples."))

# Create DGEList object...
y_DGE <- DGEList(counts = data_in)

# TMM Normalization.
y_DGE <- calcNormFactors(y_DGE)

# Create sample mapping to Tissue.Genotype.
# Group WT Cortex samples and WT Striatum samples together.
traits <- subset(alltraits, rownames(traits) %in% colnames(data_in))
traits <- traits[match(colnames(data_in), rownames(traits)), ]
if (!all(traits$SampleID == colnames(data_in))) { stop() }
# Save results to file.
myfile <- file.path(rdatdir,paste0(output_name,"_GLM_Results.RData"))
saveRDS(glm_results, myfile)

# Save results to file as spreadsheet.
myfile <- file.path(tabsdir,paste0(output_name,"_TMT_Results.xlsx"))
write_excel(glm_results, myfile)
