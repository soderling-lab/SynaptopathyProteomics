#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Parse command line arguments:
if (interactive()) {
	## If interactive:
	# User defined parameters (you only need to change these two):
	analysis_type = "Striatum" # Tissue type for analysis.
} else if (!interactive()) {
	## If not interactive, check that only 1 arg is passed.
	args <- commandArgs(trailingOnly=TRUE)
	if (length(args) == 1) { 
		analysis_type = commandArgs(trailingOnly=TRUE)[1]
		message(paste("\nAnalyzing",analysis_type,"..."))
	} else { 
		stop("Specify either 'Cortex' or 'Striatum'.",call.=FALSE) 
	}
}

## Optional parameters:
KW_alpha = 0.1
DT_alpha = 0.1

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root)

# Global imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(WGCNA)
	library(data.table)
})

# Functions.
suppressWarnings({ devtools::load_all() })

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#--------------------------------------------------------------------
## Load the data.
#--------------------------------------------------------------------

# Load the cortex and striatum data from root/data.
data(cortex_data)
data(striatum_data)

# Load the graph partitions:
data(cortex_partition)
data(striatum_partition)

# Grab the data for the tissue type we are analyzing.
data_list <- list("Cortex"=cortex_data,
		  "Striatum"=striatum_data)[[analysis_type]]
partition <- list("Cortex"=cortex_partition,
		  "Striatum"=striatum_partition)[[analysis_type]]

# Reset partition index for self-preserved modules.
partition <- reset_index(partition)

# Data matrix:
dm <- data_list$Data

# Load the sample meta data.
data(samples)

# Load adjacency matrix--coerce to a data.matrix.
adjm <- data_list$Adjm

# Load network--coerce to a matrix.
netw <- data_list$Netw

# Load glm statistical results.
myfile <- file.path(rdatdir, "GLM_Stats.RData")
glm_stats <- readRDS(myfile)

#---------------------------------------------------------------------
## Collect all modules in a list.
#---------------------------------------------------------------------

# Create list of modules.
modules <- split(partition,partition)
names(modules) <- paste0("M",names(modules))

# Remove M0.
modules <- modules[-which(names(modules) == "M0")]

# Module sizes.
module_sizes <- sapply(modules,length)

# Total Number of modules.
# P.values will be corrected for n comparisions.
nModules <- length(modules)
message(paste("Total number of modules:", nModules))

# Drop proteins assigned to M0 from partition.
part <- partition[!partition == 0]

# Percent clustered.
percent_clust <- length(part)/length(partition)
message(paste0("Percent proteins clustered: ",
	      round(100*percent_clust,2),"%."))

#---------------------------------------------------------------------
## Explore changes in module summary expression.
#---------------------------------------------------------------------

## Calculate Module Eigengenes.
# NOTE: Soft power does not influence MEs.
# NOTE: Do not need to sort partition to be in the same order as dm!
MEdata_list <- moduleEigengenes(t(dm), colors = partition, 
			    excludeGrey = TRUE, softPower = 1, impute = FALSE)
ME_dm <- as.matrix(MEdata_list$eigengenes)

# Create list of MEs.
ME_list <- setNames(object = lapply(seq(ncol(ME_dm)), function(x) ME_dm[, x]),
		    nm = names(modules))

# Module membership (KME).
data_KME <- signedKME(t(dm), MEdata_list$eigengenes, 
		      corFnc = "bicor", outputColumnName = "M")
KME_list <- lapply(seq(ncol(data_KME)), function(x) {
  v <- vector("numeric", length = nrow(data_KME))
  names(v) <- rownames(data_KME)
  v[] <- data_KME[[x]]
  v <- v[order(v, decreasing = TRUE)]
  return(v)
})
names(KME_list) <- colnames(data_KME)

# Get Percent Variance explained (PVE).
PVE <- as.numeric(MEdata_list$varExplained)
names(PVE) <- names(modules)
medianPVE <- median(PVE)
message(paste0(
  "Median module coherence (PVE): ",
  round(100 * medianPVE, 2), "%."
))

# Sample to group mapping for statistical testing.
all_groups <- as.character(interaction(samples$Model,
				       samples$SampleType))
all_groups <- trimws(all_groups)
names(all_groups) <- samples$SampleID
groups <- all_groups[rownames(ME_dm)]

# Combine WT.
groups[grep("WT",groups)] <- "WT"

# Perform Kruskal Wallis tests to identify modules whose summary
# expression profile is changing.
KW_list <- lapply(ME_list, function(x) {
			  # x <- ME_list[[1]]
			  kruskal.test(x ~ as.factor(groups[names(x)]))
})

# Combine into a data.table.
KW_dt <- as.data.table(do.call(rbind, KW_list),keep.rownames="Module") %>% 
	select(c(1,2,3,4)) # Drop 'method' and 'data.name'
colnames(KW_dt)[4] <- "pval"

# Correct p-values for n module comparisons.
KW_dt$padj <- as.numeric(KW_dt$pval) * nModules
KW_dt$padj[KW_dt$padj > 1.0] <- 1

# Clean-up column names.
colnames(KW_dt)[c(-1)] <- paste0("KW.",colnames(KW_dt)[c(-1)])

# Significant modules.
KWsig_dt <- KW_dt %>% filter(KW.padj < KW_alpha)
sigModules <- KWsig_dt[["Module"]]
nSigModules <- length(sigModules)

# Status.
message(paste0(
  "Number of modules with significant (p.adj < ", KW_alpha, ")",
  " Kruskal-Wallis test: ", nSigModules, "."
))

# Perform Dunnetts test for post-hoc comparisons.
# Note: P-values returned by DunnettTest have already been adjusted for
# multiple comparisons!

# Define control group and levels (order) for DunnettTest.
control_group <- unique(groups[grep("WT",groups)])

# Loop to perform DTest. 
# NOTE: This takes several seconds.
DT_list <- lapply(ME_list, function(x) {
  # x <- ME_list[[1]]
  group_order <- c("WT","Shank2.KO","Shank3.KO","Syngap1.HET","Ube3a.KO")
  g <- as.factor(groups[names(x)])
  levels(g) <- group_order
  DT_test <- DescTools::DunnettTest(x ~ g, control = control_group)
  DT_test$data.name <- NULL
  DT_result <- as.data.table(do.call(rbind,DT_test),keep.rownames="Contrast")
  return(DT_result)
})

# Collect DT results as data.table.
DT_dt <- bind_rows(DT_list,.id="Module")

# For every column, add "DT." to its name, except for the columns 
# named 'Module' and 'Contrast'.
colnames(DT_dt) <- gsub("WT.","",colnames(DT_dt))
idx <- which(colnames(DT_dt) %notin% c("Module","Contrast"))
colnames(DT_dt)[idx] <- paste0("DT.",colnames(DT_dt)[idx])

# Summarize number of sig tests.
nSigDT <- sapply(DT_list, function(x) sum(x$pval < DT_alpha))
message(paste0(
  "Number of significant (p.adj < ", DT_alpha, ")",
  " Dunnett's test post-hoc test(s): ", nSigModules, "."
))
idx <- nSigDT[sigModules] > 0
knitr::kable(t(nSigDT[sigModules][idx]))

# Combine modules stats.
results <- left_join(KW_dt,DT_dt,by="Module") %>% as.data.table()

# Module Size.
results$"Size" <- module_sizes[results$Module]

# PVE.
results$"PVE" <- PVE[results$Module]

# NSigDT
results$"nSig DT" <- nSigDT[results$Module]

# Sig?
is_sig <- results$KW.padj < KW_alpha & results$DT.pval < DT_alpha
results$"SigKW & SigDT" <- is_sig

#--------------------------------------------------------------------
## Save results.
#--------------------------------------------------------------------

# Annotate results with proteins.
results$Proteins <- sapply(modules[results$Module],function(x) {
				   paste(names(x),collapse=";") })

# Save to file.
myfile <- file.path(tabsdir,paste0(analysis_type,"_Module_stats.csv"))
results %>% as.data.table() %>% fwrite(myfile)

message("\nDone!")
