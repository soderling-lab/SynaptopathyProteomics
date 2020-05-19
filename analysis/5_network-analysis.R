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
	analysis_type = "Cortex" # Tissue type for analysis.
} else if (!interactive()) {
	## If not interactive, check that only 1 arg is passed.
	args <- commandArgs(trailingOnly=TRUE)
	if (length(args) == 1) { 
		analysis_type = commandArgs(trailingOnly=TRUE)[1]
	} else { 
		stop("Specify either 'Cortex' or 'Striatum'.",call.=FALSE) 
	}
}

## Optional parameters:
KW_alpha = 0.1
DT_alpha = 0.1

## Input data in root/rdata:
input_data <- list("Cortex" = list(
				   adjm = "Cortex_Adjm.csv",
				   netw = "Cortex_NE_Adjm.csv",
				   ppis = "Cortex_PPI_Adjm.csv",
				   stat = "Cortex_glm_stats.csv",
				   gmap = "Cortex_gene_map.RData",
			   	   data = "Cortex_final_protein.csv",
				   part = "Cortex_NE_SurpriseVertexPartition.csv",
				   pres = "Cortex_partition_self_preservation_enforced.csv"),
		   "Striatum" = list(
				     adjm = "Striatum_Adjm.csv",
				     netw = "Striatum_NE_Adjm.csv",
				     ppis = "Striatum_PPI_Adjm.csv",
				     stat = "Striatum_glm_stats.csv",
				     gmap = "Striatum_gene_map.RData",
				     data = "Striatum_final_protein.csv",
				     part = "Striatum_NE_SurpriseVertexPartition.csv",
				     pres = "Striatum_partition_self_preservation_enforced.csv")
		   )[[analysis_type]]

## Sample meta data in root/data:
input_meta <- list("Cortex" = "Cortex_Samples.csv",
		   "Striatum" = "Striatum_Samples.csv")[[analysis_type]]

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
TBmiscr::load_all()

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#--------------------------------------------------------------------
## Load the data.
#--------------------------------------------------------------------

message(paste0("\nAnalyzing ",analysis_type,"..."))

# Load protein expression data:
# Load the data, subset to remove QC data, coerce to matrix, 
# Log2 transform, and finally transpose such that rows = samples 
# and columns = proteins.
myfile <- file.path(rdatdir, input_data[['data']])
dm <- fread(myfile) %>% filter(Treatment != "QC") %>% 
	as.data.table() %>%
	dcast(Accession ~ Sample,value.var="Intensity") %>%
	as.matrix(rownames="Accession") %>% log2() %>% t()

# Load adjacency matrix--coerce to a data.matrix.
myfile <- file.path(rdatdir, input_data[['adjm']])
adjm <- fread(myfile) %>% as.matrix(rownames="Accession")

# Load network--coerce to a matrix.
myfile <- file.path(rdatdir, input_data[['netw']])
netw <- fread(myfile) %>% as.matrix(rownames="Accession")

# Load Leidenalg graph partition.
# This is the intial partition of the graph.
myfile <- file.path(rdatdir, input_data[['part']])
part_dt <- fread(myfile, drop=1)
la_partition <- as.numeric(part_dt) + 1
names(la_partition) <- colnames(part_dt)
n_resolutions <- nrow(part_dt)

# Load graph partition after enforcing module self-preservation.
myfile <- file.path(rdatdir, input_data[['pres']])
part_dt <- fread(myfile)
partition <- as.numeric(part_dt)
names(partition) <- colnames(part_dt)

# Reset partition index for self-preserved modules.
partition <- reset_index(partition)

# Load gene identifier map.
#myfile <- file.path(rdatdir, input_data[['gmap']])
#gene_map <- readRDS(myfile)

# Load glm statistical results.
myfile <- file.path(rdatdir, input_data[["stat"]])
glm_stats <- fread(myfile)

# Load sample info.
myfile <- file.path(datadir, input_meta)
samples <- fread(myfile)

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
MEdata_list <- moduleEigengenes(dm, colors = partition, 
			    excludeGrey = TRUE, softPower = 1, impute = FALSE)
ME_dm <- as.matrix(MEdata_list$eigengenes)


# Create list of MEs.
ME_list <- setNames(object = lapply(seq(ncol(ME_dm)), function(x) ME_dm[, x]),
		    nm = names(modules))

# Module membership (KME).
data_KME <- signedKME(dm, MEdata_list$eigengenes, 
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
all_groups <- as.character(interaction(samples$Genotype,
				       samples$Treatment))
names(all_groups) <- samples$Sample
groups <- all_groups[rownames(ME_dm)]

# Combine WT.
groups[grep("WT",groups)] <- "WT"

# Coerce to a factor.
group_levels <- c("WT","Shank2.WT","Shank3.WT","Syngap1.HET","Ube3a.KO")

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
control_group <- "WT"

# Loop to perform DTest. 
# NOTE: This takes several seconds.
DT_list <- lapply(ME_list, function(x) {
  # x <- ME_list[[1]]
  g <- as.factor(groups[names(x)])
  levels(g) <- group_levels
  return({ 
	  DescTools::DunnettTest(x ~ g, control = control_group)[[control_group]] %>% 
		  as.data.table(keep.rownames="Contrast") })
})

# Collect DT results as data.table.
DT_dt <- bind_rows(DT_list,.id="Module")

# For every column, add "DT." to its name, except for the column named `Module`.
idx <- which(colnames(DT_dt) %notin% c("Module","Contrast"))
colnames(DT_dt)[idx] <- paste0("DT.",colnames(DT_dt)[idx])

# Summarize number of sig tests.
DT_summary <- DT_dt %>% group_by(Module) %>% 
	summarize(nSig = sum(DT.pval < KW_alpha)) %>%  as.data.table()
#knitr::kable(DT_summary)

# Number of modules with significant KW + DT changes.
nSigDT <- sapply(DT_list, function(x) sum(x$pval < DT_alpha))
#message("Summary of Dunnett's test changes for significant modules:")
message(paste0(
  "Number of significant (p.adj < ", DT_alpha, ")",
  " Dunnett's test post-hoc test(s): ", nSigModules, "."
))
knitr::kable(t(nSigDT[sigModules]))

# Combine modules stats.
results <- left_join(KW_dt,DT_dt,by="Module")

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
##
#--------------------------------------------------------------------
save.image()
quit()

#--------------------------------------------------------------------
## Save results.
#--------------------------------------------------------------------

# Annotate results with proteins.
ids <- paste(glm_stats$Symbol,glm_stats$Accession,sep="|")
names(ids) <- glm_stats$Accession
results$Proteins <- sapply(modules[results$Module], function(x) {
				   paste(ids[names(x)],collapse=";") })

# Save to file.
myfile <- file.path(tabsdir,paste0(analysis_type,"_Module_stats.csv"))
results %>% as.data.table() %>% fwrite(myfile)

message("\nDone!")
