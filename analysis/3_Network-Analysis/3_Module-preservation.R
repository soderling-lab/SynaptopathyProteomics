#!/usr/bin/env Rscript

#' ---
#' title: Module Self-Preservation
#' description: Evaluate Module self-preservation by permutation testing.
#' authors: Tyler W Bradshaw
#' ---

# Parse command line input:
# Analysis (tissue) type: cortex (1) or striatum(2).
args <- commandArgs(trailingOnly = TRUE)
msg <- c("Please specify a tissue type to be analyzed:\n",
	 "       Choose either 'Cortex' or 'Striatum'.")
if (!length(args == 1)) { 
	stop(msg) 
} else { 
	type <- match(args[1],c("Cortex","Striatum"))
	tissue <- c("Cortex", "Striatum")[type]
	message(paste0("\nAnalyzing ",tissue,"..."))
}

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# User parameters to change:
stats = c(1,2,6,7) # Module statistics to use for permutation testing.
self = test = tissue
strength = "strong" # Criterion for preservation: strong or weak, see NOTE:.
data_file = tissue # Which data to use?
net_file = tissue # Which networks to test self preservation in?
adjm_file = tissue # Which correlation (adjm) network to use?
partition_file = tissue # Which partition file to use?
replace_negative = "zero" # How should negative weights be handled?
min_size = 5 # minimum allowable size for a module.
verbose = FALSE
nThreads = parallel::detectCores() - 1
# NOTE: strong = all preservation statistics must be significant. Weak = any.

## Permutation Statistics:
# 1. avg.weight
# 2. coherence
# 3. cor.cor - Don't use for correlation matrices! 
#              Sensitive to small changes in edge weight.
# 4. cor.degree - Don't use! Sensitive to small changes.
# 5. cor.contrib - Don't use! Sensitive to small changes.
# 6. avg.cor
# 7. avg.contrib

## Permutation Statistics:
# 1. avg.weight (average edge weight) - Calculated from network. Assumes edge
#    weights are positive.
# 2. coherence - Calculated from the data. Quantifies the percent variance
#    explained by a modules summary vector.
# 3. cor.cor (concordance of correlation structure) - Calculated from
#    correlation matrix. Sensitive to small changes in correlation coefficients.
# 4. cor.degree (concordance of weighted degree) - Calculated from network. 
#    Assumes edge weights are positive. Sensitive to small changes in 
#    weighted degree.
# 5. cor.contrib (concordance of node contribution) - Calculated from the data.
#    Sensitve to small changes in node contribution. 
# 6. avg.cor (density of correlation structure) - Calculated from correlation
#    matrix.
# 7. avg.contrib (average node contribution) - Quantifies how similar nodes are
#    to summary profile.

# Load renv.
root <- getrd()
renv::load(root)

# Global options and imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(NetRep)
  library(data.table)
})

# Additional functions.
devtools::load_all()

# Directories.
rdatdir <- file.path(root, "rdata")

# Load expression data. Transpose -> rows = samples; columns = genes.
myfile <- file.path(rdatdir, paste0(data_file, "_cleanDat.RData"))
data <- readRDS(myfile)
colNames <- rownames(data)
data <- t(data)
colnames(data) <- colNames

# Load adjmatrix.
myfile <- file.path(rdatdir, paste0(adjm_file, "_Adjm.RData"))
adjm <- as.matrix(readRDS(myfile))
rownames(adjm) <- colnames(adjm)

# Load network.
myfile <- file.path(rdatdir, paste0(net_file, "_NE_Adjm.RData"))
net <- as.matrix(readRDS(myfile))
rownames(net) <- colnames(net)

# Load Leidenalg graph partitions from 2_la-clustering.
myfiles <- c("Cortex" = "Cortex_SurpriseVertexPartition.csv",
	     "Striatum" = "Striatum_SurpriseVertexPartition.csv")
part_dt <- fread(file.path(rdatdir,myfiles[partition_file]), 
		    header=TRUE,drop = 1)
n_res <- nrow(part_dt)

# Check that all columns in the data are in adjm and network.
out1 <- colnames(data) %notin% colnames(adjm)
out2 <- colnames(data) %notin% colnames(net)
if (sum(out1)>0 | sum(out2)>0) { 
	message("Warning: removing columns from data that are not in network.") 
}

# Enforce consistent dimensions between data and adjm.
idx <- idy <- match(colnames(data),colnames(adjm))
adjm <- adjm[idx,idy]
check <- all(colnames(data) == colnames(adjm))
if (!check) { message("Problem: data doesn't match correlation matrix!") }

# Enforce consistent dimensions between data and partitions.
idy <- match(colnames(data),colnames(part_dt))
part_dt <- as.data.frame(part_dt)[,idy]
check <- all(colnames(data) == colnames(part_dt))
if (!check) { message("Problem: data doesn't match part_dt!") }

# Enforce consistent dimensions between data and network.
idz <- match(colnames(data),colnames(net))
net <- net[idz,idz]
check <- all(colnames(data) == colnames(net))
if (!check) { message("Problem: data doesn't match network!") }

# Final check.
check <- all(colnames(part_dt) == colnames(data) & colnames(data) == colnames(net))

#-------------------------------------------------------------------------------
## Permutation testing.
#-------------------------------------------------------------------------------

# Input for NetRep:
# 1. Expression data.
# 2. Correlation matrix or interaction network (co-expr, GO, or PPI).
# 3. Interaction network - Edges must be positive! Should network be weighted?
# 4. Network partitions.
data_list <- list(self = data)
correlation_list <- list(self = adjm)

# Networks (edges) should be positive...
if (replace_negative == "absolute value") {
	# Replace negative edges as absolute value.
	network_list <- list(self = abs(net))
} else if (replace_negative == "zero") {
	net[net<0] <- 0
	network_list <- list(self = net)
}

# Module preservation stats.
module_stats <- paste(c(
  "avg.weight", "coherence", "cor.cor", "cor.degree",
  "cor.contrib", "avg.cor", "avg.contrib"
)[stats], collapse = ", ")

# Status report:
message(paste("Evaluating self-preservation of",
		self, "modules in the",test,"network."))
message(paste0(
  "Module statistic(s) used to evaluate module preservation: ",
  module_stats), ".")
message(paste0(
  "Criterion for module preservation/divergence: ",
  strength, ".", "\n"
))

# Loop through partitions, evaluating self-preservation.
results <- list()
for (resolution in seq(n_res)) {
  message(paste("Working on partition", resolution, 
		"of", n_res, "..."))
  # Get partition.
  partition <- as.integer(part_dt[resolution, ]) 
  # Add 1 so that all module assignments >0.
  if (min(partition)==0) { partition <- partition + 1 }
  partition <- reset_index(partition)
  names(partition) <- colnames(part_dt)
  # Remove modules less than min size.
  too_small <- which(table(partition) < min_size)
  partition[partition %in% too_small] <- 0 
  module_list <- list(self = partition)
  # Perform permutation test for module self-preservation.
  suppressWarnings({
    selfPreservation <- NetRep::modulePreservation(
      network = network_list,
      data = data_list,
      correlation = correlation_list,
      moduleAssignments = module_list,
      modules = NULL,
      backgroundLabel = 0,
      discovery = "self",
      test = "self",
      selfPreservation = TRUE,
      nThreads = nThreads,
      nPerm = NULL, # Determined automatically by the function.
      null = "overlap",
      alternative = "greater", # Greater for self-preservation.
      simplify = TRUE,
      verbose = verbose
    )
  })
  # Remove NS modules--set NS modules to 0.
  preservedParts <- check_modules(selfPreservation, strength, stats)
  nModules <- length(unique(partition))
  out <- names(preservedParts)[preservedParts == "ns"]
  partition[partition %in% out] <- 0
  nPreserved <- nModules - length(out)
  message(paste("...", nPreserved, "of", nModules, "modules are preserved."))
  # Return results.
  results[[resolution]] <- partition
  # Save to Rdata.
  if (resolution == n_res) {
    output_name <- paste0(tissue, "_Self_Preservation.RData")
    saveRDS(results, file.path(rdatdir, output_name))
    message("Done!")
  }
} # Ends loop.
