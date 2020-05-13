#!/usr/bin/env Rscript

#' ---
#' title: 1_TAMPOR.R
#' description: TAMPOR Normalization of preprocessed TMT data.
#' authors: Tyler W Bradshaw, Eric B Dammer (TAMPOR).
#' ---

## Inputs:
# Input data should be in root/rdata:

# Sample meta data.
input_samples = list("Cortex" = "4227_TMT_Cortex_Combined_traits.csv",
		     "Striatum" = "4227_TMT_Striatum_Combined_traits.csv") 

# Preprocessed expression data.
input_data = list("Cortex" = "Cortex_preprocessed.RData",
		  "Striatum" = "Striatum_preprocessed.RData")

# Gene mapping data.
input_maps = list("Cortex" = "Cortex_gene_map.RData",
		  "Striatum" = "Striatum_gene_map.RData")

## Other parameters:
output_name = "Combined"
alpha_threshold = 0.1 # FDR threshold for signficance.
n_threads = parallel::detectCores() - 1 # Number of cores for parallel processing.

## Main Outputs:
# Stored in root/tables/
# 0. [output_name]_GLM_Results.xlsx - Statistical results.

## Output for downstream analysis:
# Stored in root/rdata/
# 0. [output_name]_gene_map.RData   - gene identifier map.
# 1. [output_name]_tidy_protein.csv - tidy, final, normalized protein data.
# 2. [output_name]_
# 3. [output_name]_

## Order of data processing operations:

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. 
# Load custom functions and prepare the project directory for saving 
# output files.

# Load renv.
renv::load(getrd())

# Load required packages.
suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(tibble)
  library(data.table)
})

# Load additional functions:
TBmiscr::load_all()

# Set any other directories.
root <- getrd()
funcdir <- file.path(root, "R")
figsdir <- file.path(root,"figs")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#---------------------------------------------------------------------
## Combine gene mapping data.
#---------------------------------------------------------------------

message("\nCombining gene identifier maps...")

# Load gene mapping data.
myfiles <- sapply(input_maps, function(f) file.path(rdatdir, f))
gene_maps <- lapply(myfiles,readRDS)

# Combine.
gene_map <- dplyr::bind_rows(gene_maps) %>% unique()

#---------------------------------------------------------------------
## Combine cortex and striatum traits.
#---------------------------------------------------------------------
# We will utilize TAMPOR to combine the Cortex and Striatum datasets.
# Merge the preprocessed data and traits files.

message("\nPreparing sample meta data for TAMPOR normalization...")

# Load the sample info into a list, traits.
myfiles <- sapply(input_samples,function(f) file.path(datadir, f))
traits <- lapply(myfiles, fread)

# Bind the elements of the list into a df.
traits <- do.call(rbind, traits)

# SampleIDs are batch.channel.
# This doesn't work???
batch <- paste0("b",as.numeric(interaction(traits$Tissue,traits$Genotype)))
#batch <- paste0("b",as.numeric(as.factor(paste(traits$Tissue,traits$Genotype))))
channel <- traits$Channel
traits$SampleID <- paste(batch,channel,sep=".")

# Traits must include Batch column:
traits$Batch <- batch

# Insure that rownames are new sampleIDS
rownames(traits) <- traits$SampleID

#---------------------------------------------------------------------
## Combine cortex and striatum expression data.
#---------------------------------------------------------------------
# We will utilize TAMPOR to combine the Cortex and Striatum datasets.
# Merge the preprocessed data and traits files.

message("\nPreparing protein expression data TAMPOR normalization...")

# Load the Cortex and Striatum cleanDat.
myfiles <- sapply(input_data,function(f) file.path(rdatdir, f))
data_list <- lapply(myfiles, readRDS)

# Bind data frames together.
data_merge <- dplyr::bind_rows(data_list) %>% as.data.table()

# Add SampleID annotation.
data_merge$SampleID <- traits$SampleID[match(data_merge$Sample,traits$Sample)]

# Cast into a data.matrix.
dm <- data_merge %>% 
	dcast(Accession ~ SampleID, value.var="Intensity") %>% 
	as.matrix(rownames="Accession")

# Drop QC samples.
qc_samples <- traits %>% filter(Treatment == "QC")
out <- colnames(dm) %in% qc_samples$SampleID
dm <- dm[,!out]

# Remove rows with any missing values.
missing_vals <- apply(dm,1,function(x) any(is.na(x)))
message(paste("\nRemoving", sum(missing_vals), 
	      "rows that contain missing values."))
combDat <- dm[!missing_vals,]

#---------------------------------------------------------------------
## Perform TAMPOR normalization.
#---------------------------------------------------------------------

message(paste("\nPerforming TAMPOR normalization to merge Cortex and",
	      "Striatum datasets..."))

# Insure than any samples that were removed from cleanDat are removed from
# traits (outliers and QC samples).
sub_traits <- traits[rownames(traits) %in% colnames(combDat), ]

# Rownames of traits must match column_names(data).
rownames(sub_traits) <- sub_traits$SampleID

# GIS index for TAMPOR normalization is all WT samples.
# WT Cortex and WT Striatum samples will be scaled by TAMPOR to be 
# ~ equal.
controls <- sub_traits$SampleID[grepl("WT", sub_traits$Treatment)]

# For some reason, traits can only contain these columns:
sub_traits <- as.data.frame(sub_traits %>% select(SampleID,Batch))
rownames(sub_traits) <- sub_traits$SampleID

# Perform TAMPOR normalization.
# Ignore warnings about closing connections.
results <- TAMPOR(
  dat = combDat,
  traits = sub_traits,
  batchPrefixInSampleNames = TRUE,
  samplesToIgnore = "None",
  GISchannels = controls,
  parallelThreads = n_threads
)

# Collect normalized, relative abundance data.
cleanDat <- results$cleanRelAbun

#---------------------------------------------------------------------
## Identify and remove any sample outliers.
#---------------------------------------------------------------------

# There are none at threshold = -3.0

#---------------------------------------------------------------------
## Protein differential abundance.
#---------------------------------------------------------------------

# Status.
message(paste("\nAnalyzing protein differential abundance..."))

# We will analyze cleanDat.
dm <- cleanDat

# Create dge object, perform final TMM normalization.
dge <- DGEList(counts=dm)
dge <- calcNormFactors(dge)

# Extract final normalized data from dge object.

# SampleID to group mapping.
groups <- paste(traits$Tissue,traits$Treatment,traits$Genotype,sep=".")
groups[grep("Cortex.WT",groups)] <- "Cortex.WT"
groups[grep("Striatum.WT",groups)] <- "Striatum.WT"
names(groups) <- traits$SampleID

# Annotate DGE object with sample groups.
groups <- groups[rownames(dge$samples)]
dge$samples$group <- as.factor(groups)

# Create design matrix.
design <- model.matrix(~ 0 + groups, data = dge$samples)
colnames(design) <- levels(dge$samples$group)

# Estimate dispersion.
dge <- estimateDisp(dge,design,robust=TRUE)

# Fit a glm.
fit <- glmQLFit(dge, design, robust = TRUE)

# Make contrasts for statistical comparisons.
g1 <- colnames(design)[grepl("Cortex", colnames(design))][-5]
g2 <- colnames(design)[grepl("Striatum", colnames(design))][-5]

cont1 <- makePairwiseContrasts(list(g1), list("Cortex.WT"))
cont2 <- makePairwiseContrasts(list(g2), list("Striatum.WT"))

# NOTE: For some reason, loops or lapply doesn't work with makeContrasts()
contrasts <- list(
  makeContrasts(cont1[1], levels = design),
  makeContrasts(cont1[2], levels = design),
  makeContrasts(cont1[3], levels = design),
  makeContrasts(cont1[4], levels = design),
  makeContrasts(cont2[1], levels = design),
  makeContrasts(cont2[2], levels = design),
  makeContrasts(cont2[3], levels = design),
  makeContrasts(cont2[4], levels = design)
)
names(contrasts) <- unlist({
  lapply(contrasts, function(x) sapply(strsplit(colnames(x), " "), "[", 1))
})

# Call glmQLFTest() to evaluate differences in contrasts.
qlf <- lapply(contrasts, function(x) glmQLFTest(fit, contrast = x))
names(qlf) <- names(contrasts)

## Determine number of significant results with decideTests().
summary_list <- lapply(qlf, function(x) summary(decideTests(x)))

# Table summarizing DA proteins.
overall <- t(matrix(unlist(summary_list), nrow = 3, ncol = 8))
rownames(overall) <- unlist(lapply(contrasts, function(x) colnames(x)))
colnames(overall) <- c("Down", "NS", "Up")
overall <- as.data.frame(overall)
row_names <- sapply(strsplit(rownames(overall), " - "), "[", 1)
row_names <- gsub(".KO.|.HET.", " ", row_names)
overall <- add_column(overall, Experiment = row_names, .before = 1)
overall <- overall[, c(1, 3, 2, 4)]
overall$"Total Sig" <- rowSums(overall[, c(3, 4)])
overall <- overall[c(2, 6, 3, 7, 1, 5, 4, 8), ] # Reorder.
rownames(overall) <- NULL

# Table:
message(paste("\nSummary of differentially abundant proteins (FDR<",
	      alpha_threshold,"):"))
knitr::kable(overall)

# Call topTags to add FDR. Gather tabularized results.
f <- function(x) { return(topTags(x, n = Inf, sort.by = "none")$table) }
glm_results <- lapply(qlf, f)

# Convert logCPM column to percent WT.
f <- function(x) { x$logCPM <- 2^x$logFC; return(x) }
glm_results <- lapply(glm_results, f)

# Rename logCPM column.
f <- function(x) { colnames(x)[2] <- "PercentWT"; return(x) }
glm_results <- lapply(glm_results, f)

# Function to annotate results with gene ids.
add_ids <- function(x,gene_map) {
	Uniprot <- rownames(x)
	x <- tibble::add_column(x,Uniprot,.after=0)
	idx <- match(rownames(x),gene_map$uniprot)
	Symbol <- gene_map[["symbol"]][idx]
	Entrez <- gene_map[["entrez"]][idx]
	x <- tibble::add_column(x,Symbol,.after=1)
	x <- tibble::add_column(x,Entrez,.after=2)
	rownames(x) <- NULL
	return(x)
}

# Annotate with gene ids.
glm_results <- lapply(glm_results,function(x) add_ids(x,gene_map))

# Sort by PValue.
f <- function(x) { x[order(x$PValue),] }
glm_results <- lapply(glm_results, f)

# Add candidate column.
f <- function(x) { x$Candidate <- x$FDR < alpha_threshold ; return(x) }
glm_results <- lapply(glm_results, f)

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:
# Stored in root/rdata/
# 0. [output_name]_gene_map.RData   -- gene identifier map.
# 1. [output_name]_tidy_peptide.csv -- raw peptide data.
# 2. [output_name]_cleanDat.RData   -- data for TAMPOR.

## Save key results.
message("\nSaving data for downstream analysis...")

# 0. gene_map.RData - gene identifier map.
myfile <- file.path(rdatdir,paste(output_name,"gene_map.RData",sep="_"))
saveRDS(gene_map,myfile)

# 1. EdgeR statistical results.
myfile <- file.path(tabsdir,paste(output_name,"GLM_Results.xlsx",sep="_"))
write_excel(glm_results,myfile)

# 2. [output_name]_tidy_protein.csv -- tidy, final normalized data.
#myfile <- file.path(rdatdir,paste(output_name,"tidy_peptide.csv",sep="_"))
#fwrite(tidy_peptide, myfile)

# 3. [output_name]_norm_protein.csv -- final normalized data matrix. 
#myfile <- file.path(rdatdir,paste(output_name,"cleanDat.RData",sep="_"))
#fwrite(tidy_protein, myfile)

message("\nDone!")
