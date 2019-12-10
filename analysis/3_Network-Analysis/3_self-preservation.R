#!/usr/bin/env Rscript

#' ---
#' title: Module Self-Preservation
#' description: Evaluate Module self-preservation by permutation testing.
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# NetRep input:
# 1. Expression data.
# 2. Correlation matrix.
# 3. Interaction network - Edges must be positive! Should network be weighted?
# 4. Network partitions.

## Permutation Statistics:
# 1. avg.weight (average edge weight) - Calculated from network. Assumes edge
#    weights are positive.
# 2. coherence - Calculated from the data. Quantifies the percent variance
#    explained by a modules summary vector.
# 3. cor.cor (concordance of correlation structure) - Calculated from
#    correlation matrix.
# 4. cor.degree (concordance of weighted degree) - Calculated from network. Assumes edge weights are
#    positive.
# 5. cor.contrib (average node contribution) - Calculated from the data.
# 6. avg.cor (density of correlation structure) - Calculated from correlation
#    matrix.
# 7. avg.contrib (average node contribution) - Quantifies how similar nodes are
#    to summary profile.

# User parameters to change:
stats <- c(1:7) # Module statistics to use for permutation testing.
strength <- "strong" # Criterion for preservation: strong = ALL, weak = ANY sig stats.
weighted <- FALSE # Weighted or unweighted. If TRUE, then appropriate soft-power will be calculated.
self <- "Striatum" # Which networks to test self preservation in? #self = c("wt","ko","cortex","striatum","combined")
nres <- 100 # Total number of resolutions to be anlyzed.

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
#invisible(sapply(myfun, source))
sapply(myfun,function(x) source(x,encoding = "UTF-8"))

quit()

# Load expression data. Transpose -> rows = samples; columns = genes.
myfile <- file.path(datadir, paste0("3_", self, "_cleanDat.RData"))
data <- t(readRDS(myfile))

# Load adjmatrix.
myfile <- file.path(datadir, paste0("3_", self, "_Adjm.RData"))
adjm <- as.matrix(readRDS(myfile))
rownames(adjm) <- colnames(adjm)

# Load Leidenalg graph partitions from 2_la-clustering.
myfile <- file.path(datadir, paste0("3_", self, "_partitions.csv"))
partitions <- data.table::fread(myfile, drop = 1, skip = 1)
colnames(partitions) <- colnames(adjm)

# Weighted or unweighted?
# If weighted, then calculate power for approximate scale free fit.
if (weighted) {
  message("Calculating soft-power for weighting co-expression graph!")
  sft <- silently({
    pickSoftThreshold(data,
      corFnc = "bicor",
      networkType = "signed",
      RsquaredCut = 0.8
    )$powerEstimate
  })
  names(sft) <- self
} else {
  # Power = 1 == Unweighted
  sft <- 1
  names(sft) <- self
}

#-------------------------------------------------------------------------------
## Permutation testing.
#-------------------------------------------------------------------------------

# Input for NetRep:
# Networks (edges) should be positive -> AbsoluteValue()
data_list <- list(self = data)
correlation_list <- list(self = adjm)
network_list <- list(self = abs(adjm^sft[self]))

# Loop through partitions, evaluating self-preservation.
results <- list()
for (i in 1:nres) {
  # Status report.
  message(paste("Working on partition", i, "of", nres, "..."))
  # Get partition--adding 1 so that all module assignments >0.
  partition <- as.integer(partitions[i, ]) + 1
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
      # nPerm = 100000, # Determined automatically by the function.
      null = "overlap",
      alternative = "greater", # Greater for self-preservation.
      simplify = TRUE,
      verbose = FALSE
    )
  })
  # Remove NS modules--set NS modules to 0.
  preservedParts <- check_modules(selfPreservation, strength, stats)
  nModules <- length(unique(partition))
  out <- names(preservedParts)[preservedParts == "ns"]
  partition[partition %in% out] <- 0
  nPreserved <- nModules - length(out)
  message(paste("...", nPreserved, "modules of", nModules, "modules are preserved."))
  # Return results.
  results[[i]] <- partition
  # Save to Rdata.
  if (i == nres) {
    output_name <- paste0(jobID, "_", self, "_Module_Self_Preservation.RDS")
    saveRDS(results, file.path(datadir, output_name))
    message("Done!")
  }
} # Ends loop.
