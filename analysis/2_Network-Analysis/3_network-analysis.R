#!/usr/bin/env Rscript

#' ---
#' title: Network Analysis
#' description: Protein co-expression network analysis
#' authors: Tyler W Bradshaw
#' ---

## User parameters to change:
analysis_type = "Cortex"
root = "/mnt/d/projects/SynaptopathyProteomics" 

## Optional parameters:
alpha_KW = 0.1
alpha_DT = 0.05

## Input data in root/rdata:
input_data <- list("Cortex" = list(
				   adjm = "Cortex_Adjm.csv",
				   netw = "Cortex_NE_Adjm.csv",
				   stat = "Cortex_glm_stats.csv",
				   gmap = "Cortex_gene_map.RData",
			   	   data = "Cortex_norm_protein.csv",
				   part = "Cortex_NE_SurpriseVertexPartition.csv",
				   pres = "Cortex_partition_self_preservation_enforced.csv"),
		   "Striatum" = list(
				     adjm = "Striatum_Adjm.csv",
				     netw = "Striatum_NE_Adjm.csv",
				     stat = "Striatum_glm_stats.csv",
				     gmap = "Striatum_gene_map.RData",
				     data = "Striatum_norm_protein.csv",
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

# Calculate Module Eigengenes.
# NOTE: Soft power does not influence MEs.
# NOTE: Do not need to sort partition to be in the same order as dm!
data_ME <- moduleEigengenes(dm, colors = partition, 
			    excludeGrey = TRUE, softPower = 1, impute = FALSE)
MEs <- as.matrix(data_ME$eigengenes)

# Create list of MEs.
ME_list <- lapply(seq(ncol(MEs)), function(x) MEs[, x])
names(ME_list) <- names(modules)

# Module membership (KME).
data_KME <- signedKME(dm, data_ME$eigengenes, 
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
PVE <- as.numeric(data_ME$varExplained)
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
groups <- all_groups[rownames(MEs)]

# Perform Kruskal Wallis tests to identify modules whose summary
# expression profile is changing.
KW_list <- lapply(ME_list, function(x) {
  kruskal.test(x ~ groups[names(x)])
})
KW_dt <- as.data.table(do.call(rbind, KW_list),keep.rownames="Module") %>% 
	select(c(1,2,3,4))
colnames(KW_dt)[4] <- "pval"

# Correct p-values for n module comparisons.
KW_dt$padj <- as.numeric(KW_dt$pval) * nModules
KW_dt$padj[KW_dt$padj > 1.0] <- 1

# Clean-up column names.
colnames(KW_dt)[c(-1)] <- paste0("KW.",colnames(KW_dt)[c(-1)])

# Significant modules.
sigKW_dt <- KW_dt %>% filter(KW.padj < alpha_KW)
sigModules <- sigKW_dt[["Module"]]
nSigModules <- length(sigModules)

# Status.
message(paste0(
  "Number of modules with significant (p.adj < ", alpha_KW, ")",
  " Kruskal-Wallis test: ", nSigModules, "."
))

# Perform Dunnetts test for post-hoc comparisons.
# Note: P-values returned by DunnettTest have already been adjusted for
# multiple comparisons!

# Define control group and levels (order) for DunnettTest.
group_order <- c("Shank2.WT", "Shank2.KO", "Shank3.WT", "Shank3.KO",
		 "Syngap1.WT","Syngap1.HET","Ube3a.WT","Ube3a.HET")

control_groups <- paste(c("Shank2","Shank3","Syngap1","Ube3a"),"WT", sep = ".")

# Loop to perform DTest. 
# NOTE: This takes several seconds.
DT_list <- lapply(ME_list, function(x) {

  g <- factor(groups[names(x)], levels = group_order)

  dt_res <- DescTools::DunnettTest(x ~ g, control = control_groups)

  result <- as.data.table(dt_res[[control_group]],keep.rownames="Contrast")
  return(result)
})

# Collect DT results as data.table.
DT_dt <- bind_rows(DT_list,.id="Module")
colnames(DT_dt)[c(-1)] <- paste0("DT.",colnames(DT_dt)[c(-1)])

# Number of modules with significant KW + DT changes.
nSig_tests <- sapply(DT_list, function(x) sum(x$pval < alpha_DT))
message("Summary of Dunnett's test changes for significant modules:")
nSig_tests[sigModules]

# Combine modules stats.
results <- left_join(KW_dt,DT_dt,by="Module")
colnames(results)

# Module Size.
results$"Size" <- module_sizes[results$Module]

# PVE.
results$"PVE" <- PVE[results$Module]

# Sig?
results$"Sig" <- results$KW.padj < alpha_KW & results$DT.pval < alpha_DT

#--------------------------------------------------------------------
## Save results.
#--------------------------------------------------------------------

# Annotate results with proteins.
ids <- paste(glm_stats$Symbol,glm_stats$Accession,sep="|")
names(ids) <- glm_stats$Accession
results$Proteins <- sapply(modules[results$Module], function(x) {
				   paste(ids[names(x)],collapse=";") })

# Save to file.
myfile <- file.path(tabsdir)
results %>% as.data.table() %>% fwrite(myfile)

message("\nDone!")
