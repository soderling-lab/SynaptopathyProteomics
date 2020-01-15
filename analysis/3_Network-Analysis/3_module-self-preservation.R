#!/usr/bin/env Rscript

#' ---
#' title: Module Self-Preservation
#' description: Evaluate Module self-preservation by permutation testing.
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change:
stats <- c(1,2,6,7) # Module statistics to use for permutation testing.
strength <- "strong" # Criterion for preservation: strong = ALL, weak = ANY sig stats.
self <- "Striatum" # Which networks to test self preservation in? #self = c("wt","ko","Cortex","Striatum","Sombined", "PPI", "GO")
nres <- 100 # Total number of resolutions to be anlyzed.
verbose <- FALSE

# NetRep input:
# 1. Expression data.
# 2. Correlation matrix or interaction network (GO or PPI).
# 3. Interaction network - Edges must be positive! Should network be weighted?
# 4. Network partitions.

## Permutation Statistics:
# 1. avg.weight
# 2. coherence
# 3. cor.cor - Don't use! Sensitive to small changes.
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

# Is this a slurm job?
slurm <- any(grepl("SLURM", names(Sys.getenv())))
if (slurm) {
  # SLURM job notes - sent to job_*.info
  nThreads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
  jobID <- as.integer(Sys.getenv("SLURM_JOBID"))
  info <- as.matrix(Sys.getenv())
  idx <- grepl("SLURM", rownames(info))
  myfile <- file.path("./out", paste0("job_", jobID, ".info"))
  write.table(info[idx, ], myfile, col.names = FALSE, quote = FALSE, sep = "\t")
} else {
  nThreads <- 8
  jobID <- ""
}

# Global options and imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(WGCNA)
  library(NetRep)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
datadir <- file.path(root, "rdata")
funcdir <- file.path(root, "R")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load expression data. Transpose -> rows = samples; columns = genes.
myfile <- file.path(datadir, paste0("3_", self, "_cleanDat.RData"))
data <- readRDS(myfile)
colNames <- rownames(data)
data <- t(data)
colnames(data) <- colNames

# Load adjmatrix or interaction network.
myfile <- file.path(datadir, paste0("3_", self, "_Adjm.RData"))
adjm <- as.matrix(readRDS(myfile))
rownames(adjm) <- colnames(adjm)

# Load Leidenalg graph partitions from 2_la-clustering.
myfiles <- c("Cortex" = file.path(datadir,"147731383_Cortex_CPMVertexPartition_partitions.csv"),
	    "Striatum" = file.path(datadir,"148436673_Striatum_CPMVertexPartition_partitions.csv"))
partitions <- data.table::fread(myfiles[self], header=TRUE,drop = 1)

# Enforce consistent dimensions between data and adjm.
# Remove duplicate column from data.
out <- colnames(data) %notin% colnames(adjm)
data <- data[,!out]

# Enforce consistent order between data, adjm, and partitions.
idx <- idy <- match(colnames(data),colnames(adjm))
adjm <- adjm[idx,idy]
check <- all(colnames(data) == colnames(adjm))
if (!check) { message("Problem!") }
idy <- match(colnames(data),colnames(partitions))
partitions <- as.data.frame(partitions)[,idy]
if (!check) { message("Problem!") }
check <- all(colnames(data) == colnames(partitions))

#-------------------------------------------------------------------------------
## Permutation testing.
#-------------------------------------------------------------------------------

# Input for NetRep:
# Networks (edges) should be positive -> AbsoluteValue()
data_list <- list(self = data)
correlation_list <- list(self = adjm)
network_list <- list(self = abs(adjm))

# Module preservation stats.
module_stats <- paste(c(
  "avg.weight", "coherence", "cor.cor", "cor.degree",
  "cor.contrib", "avg.cor", "avg.contrib"
)[stats], collapse = ", ")

# Status report:
message(paste("Evaluating self-preservation of",
		self, "modules."))
message(paste0(
  "Module statistic(s) used to evaluate module preservation: ",
  module_stats), ".")
message(paste0(
  "Criterion for module preservation/divergence: ",
  strength, ".", "\n"
))

# Loop through partitions, evaluating self-preservation.
results <- list()
for (i in 1:nres) {
  message(paste("Working on partition", i, "of", nres, "..."))
  # Get partition--adding 1 so that all module assignments >0.
  partition <- as.integer(partitions[i, ]) + 2
  names(partition) <- colnames(adjm)
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
      verbose
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
  results[[i]] <- partition
  # Save to Rdata.
  if (i == nres) {
    output_name <- paste0(jobID, "_", self, "_Module_Self_Preservation.RData")
    saveRDS(results, file.path(datadir, output_name))
    message("Done!")
  }
} # Ends loop.
