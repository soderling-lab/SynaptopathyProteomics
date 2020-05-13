#!/usr/bin/env Rscript

#' ---
#' title: 1_TAMPOR.R
#' description: TAMPOR Normalization of preprocessed TMT data.
#' authors: Tyler W Bradshaw, Eric B Dammer (TAMPOR).
#' ---

## Inputs:
# Input data should be in root/rdata:
input_samples = list("Cortex" = "4227_TMT_Cortex_Combined_traits.csv",
		     "Striatum" = "4227_TMT_Striatum_Combined_traits.csv") 
input_data = list("Cortex" = "Cortex_cleanDat.RData",
		  "Striatum" = "Striatum_cleanDat.RData")

## Other parameters:
output_name = "Combined"

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
  library(tibble)
  library(data.table)
  library(ggplot2)
  library(edgeR)
  library(openxlsx)
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
## Combine cortex and striatum traits.
#---------------------------------------------------------------------
# We will utilize TAMPOR to combine the Cortex and Striatum datasets.
# Merge the preprocessed data and traits files.

# Load the sample info into a list, traits.
myfiles <- sapply(input_samples,function(f) file.path(datadir, f))
traits <- lapply(myfiles, fread)

# Bind the elements of the list into a df.
traits <- do.call(rbind, traits)

# SampleIDs are batch.channel.
#batch <- paste0("b",as.numeric(interaction(traits$Tissue,traits$Genotype)))
batch <- paste0("b",as.numeric(as.factor(paste(traits$Tissue,traits$Genotype))))
channel <- traits$Channel
traits$SampleID <- paste(batch,channel,sep=".")

# Traits must include:
traits$Batch <- batch

# Insure that rownames are new sampleIDS
rownames(traits) <- traits$SampleID

#---------------------------------------------------------------------
## Combine cortex and striatum expression data.
#---------------------------------------------------------------------
# We will utilize TAMPOR to combine the Cortex and Striatum datasets.
# Merge the preprocessed data and traits files.

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

# Remove rows with missing values.
missing_vals <- apply(dm,1,function(x) any(is.na(x)))
message(paste("Removing", sum(missing_vals), 
	      "rows that contain missing values."))
cleanDat <- dm[!missing_vals,]

# Save merged traits file.
myfile <- file.path(rdatdir, paste(output_name,"traits.RData",sep="_"))
saveRDS(traits, file = myfile)

#---------------------------------------------------------------------
## Perform TAMPOR normalization.
#---------------------------------------------------------------------

# Insure than any samples that were removed from cleanDat are removed from
# traits (outliers and QC samples).
sub_traits <- traits[rownames(traits) %in% colnames(cleanDat), ]

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
  dat = cleanDat,
  traits = sub_traits,
  batchPrefixInSampleNames = TRUE,
  samplesToIgnore = "None",
  GISchannels = controls,
  parallelThreads = 8
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

# We will analyze cleanDat.
dm <- cleanDat
class(dm)
dim(dm)

# Create dge object.
dge <- DGEList(counts=dm)
dge <- calcNormFactors(dge)

# SampleID to group mapping.
groups <- paste(traits$Tissue,traits$Treatment,traits$Genotype,sep=".")
groups[grep("Cortex.WT",groups)] <- "Cortex.WT"
groups[grep("Striatum.WT",groups)] <- "Striatum.WT"
#groups <- groups[!grepl("QC",groups)]
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

g1 <- colnames(design)[grepl("Cortex", colnames(design))][-5]
g2 <- colnames(design)[grepl("Striatum", colnames(design))][-5]

cont1 <- makePairwiseContrasts(list(g1), list("Cortex.WT"))
cont2 <- makePairwiseContrasts(list(g2), list("Striatum.WT"))

# Make contrasts for EdgeR.
# For some reason loops or lapply dont work with the makeContrasts function.
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

# Table of DA candidates.
# Modify tables theme to change font size.
# Cex is a scaling factor relative to the defaults.
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params = list(cex = 0.75)),
  colhead = list(fg_params = list(cex = 0.75)),
  rowhead = list(fg_params = list(cex = 0.75))
)

# Create table and add borders.
mytable <- tableGrob(overall, rows = NULL, theme = mytheme)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 2, b = nrow(mytable), l = 1, r = ncol(mytable)
)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 1, l = 1, r = ncol(mytable)
)

# Call topTags to add FDR. Gather tabularized results.
f <- function(x) { return(topTags(x, n = Inf, sort.by = "none")$table) }
glm_results <- lapply(qlf, f)

# Convert logCPM column to percent WT and annotate with candidate column.
f <- function(x) { x$logCPM <- 
glm_results <- lapply(glm_results, function(x) annotateTopTags(x))

# Use protmap to annotate glm_results with entrez Ids and gene symbols.
for (i in 1:length(glm_results)) {
  x <- glm_results[[i]]
  idx <- match(rownames(x), protmap$ids)
  x <- add_column(x, "Gene|Uniprot" = protmap$ids[idx], .before = 1)
  x <- add_column(x, "Uniprot" = protmap$uniprot[idx], .after = 1)
  x <- add_column(x, "Entrez" = protmap$entrez[idx], .after = 2)
  x <- add_column(x, "Symbol" = protmap$gene[idx], .after = 3)
  glm_results[[i]] <- x
}




# Collect the results.
tidy_protein <- results$data
glm_results <- results$results

# Summary of DA proteins:
message(paste0("Summary of differentially abundant proteins ",
	      "in each subceulluar fraction (FDR < ",alpha,"):"))
results$summary

# Proteins that are commonly dysregulated:
combined_results <- rbindlist(glm_results,idcol="Fraction")
idx <- match(combined_results$Accession,gene_map$uniprot)
combined_results$Gene <- gene_map$symbol[idx]
df <- combined_results %>% group_by(Accession) %>% 
	summarize(Gene = unique(Gene),
		  nSig = sum(FDR<0.1),
		  nFractions = length(FDR)) %>% 
	filter(nSig==nFractions)
message(paste("Proteins that are differentially abundant in all fractions:\n",
	      paste(df$Gene,collapse=", ")))

