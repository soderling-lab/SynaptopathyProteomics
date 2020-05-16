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
alpha_KW <- 0.1
alpha_DT <- 0.05

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
input_meta <- list("Cortex" = "4227_TMT_Cortex_Combined_traits.csv",
		   "Striatum" = "4227_TMT_Striatum_Combined_traits.csv")[[analysis_type]]


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
root <- TBmiscr::getrd()
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

#--------------------------------------------------------------------
## Load the data.
#--------------------------------------------------------------------

# Load expression data:
# Load the data, subset to remove QC data, coerce to matrix, 
# Log2 transform, and finally transpose such that rows = samples 
# and columns = proteins.
myfile <- file.path(rdatdir, input_data[['data']])
dm <- fread(myfile) %>% filter(Treatment != "QC") %>% 
	as.data.table() %>%
	dcast(Accession ~ Sample,value.var="Intensity") %>%
	as.matrix(rownames="Accession") %>% log2() %>% t()

# Load adjmatrix--coerce to a matrix.
myfile <- file.path(rdatdir, input_data[['adjm']])
adjm <- fread(myfile) %>% as.matrix(rownames="Accession")

# Load network--coerce to a matrix.
myfile <- file.path(rdatdir, input_data[['netw']])
netw <- fread(myfile) %>% as.matrix(rownames="Accession")

# Load Leidenalg graph partition.
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
myfile <- file.path(rdatdir, input_data[['gmap']])
gene_map <- readRDS(myfile)

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
sapply(modules,length)

# Total Number of modules.
# P.values will be corrected for n comparisions.
nModules <- length(modules)
message(paste("Total number of modules:", nModules))

# Drop proteins assigned to M0 from partition.
part <- partition[!partition == 0]

# Percent clustered.
percent_clust <- length(part)/length(partition)
message(paste("Percent proteins clustered:",round(percent_clust,3)))

#---------------------------------------------------------------------
## Explore changes in module summary expression.
#---------------------------------------------------------------------

# Calculate Module Eigengenes.
# Note: Soft power does not influence MEs.
# Note: Do not need to sort partition to be in the same order!
data_ME <- moduleEigengenes(dm, colors = partition, 
			    excludeGrey = TRUE, softPower = 1, impute = FALSE)
MEs <- as.matrix(data_ME$eigengenes)

# Create list of MEs.
ME_list <- lapply(seq(ncol(MEs)), function(x) MEs[, x])
names(ME_list) <- names(modules)

# Module membership (KME).
data_KME <- signedKME(dm, data_ME$eigengenes, corFnc = "bicor", outputColumnName = "M")
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
message(paste(
  "Median module coherence (PVE):",
  round(100 * medianPVE, 2), "(%)."
))

# Sample to group mapping for statistical testing.
all_groups <- as.character(interaction(samples$Treatment,
				       samples$Genotype,
				       samples$Tissue))
names(all_groups) <- samples$Sample
groups <- all_groups[rownames(MEs)]

# Group all WT samples from a tissue type together.
groups[grepl("WT.*.Cortex", groups)] <- "WT.Cortex"
groups[grepl("WT.*.Striatum", groups)] <- "WT.Striatum"

# Perform Kruskal Wallis tests to identify modules whose summary
# expression profile is changing.
KW_list <- lapply(ME_list, function(x) {
  kruskal.test(x ~ groups[names(x)])
})
data_KW <- as.data.table(do.call(rbind, KW_list)) %>% select(c(1,2,3))

# Correct p-values for n comparisons.
data_KW$p.adj <- as.numeric(data_KW$p.value) * nModules
data_KW$p.adj[data_KW$p.adj > 1.0] <- 1

# Significant modules.
alpha_KW = 0.05
sigModules <- rownames(data_KW)[data_KW$p.adj < alpha_KW]
nSigModules <- length(sigModules)
message(paste0(
  "Number of modules with significant (p.adj < ", alpha_KW, ")",
  " Kruskal-Wallis test: ", nSigModules, "."
))
sigModules <- paste0("M",sigModules)

# Perform Dunnetts test for post-hoc comparisons.
# Note: P-values returned by DunnettTest have already been adjusted for
# multiple comparisons!

# Define control group and levels (order) for DunnettTest.
group_order <- c("WT", "KO.Shank2", "KO.Shank3", "HET.Syngap1", "KO.Ube3a")
group_levels <- paste(group_order, analysis_type, sep = ".")
control_group <- paste("WT", analysis_type, sep = ".")

# Loop to perform DTest. NOTE: This takes several seconds.
DT_list <- lapply(ME_list, function(x) {
  g <- factor(groups[names(x)], levels = group_levels)
  result <- DescTools::DunnettTest(x ~ g, control = control_group)[[control_group]]
  return(as.data.table(result))
})

# Number of modules with significant KW + DT changes.
nSig_tests <- sapply(DT_list, function(x) sum(x$pval < alpha_DT))
message("Summary of Dunnett's test changes for significant modules:")
nSig_tests[sigModules]

# What are the proteins???!?!
foo = modules[sigModules]
x = foo[[1]]
f <- function(x) { gene_map$symbol[match(names(x),gene_map$uniprot)] }
lapply(foo,f)
