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
self <- "Cortex" # Which networks to test self preservation in?
partition_file <- "Cortex_MCL" # Which partition file to use?
resolutions <- seq(1:100) # Resolutions to be anlyzed.
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
  jobID <- Sys.Date()
}

# Global options and imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(WGCNA)
  library(NetRep)
  library(parallel)
  library(doParallel)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatdir <- file.path(root, "rdata")
funcdir <- file.path(root, "R")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load expression data. Transpose -> rows = samples; columns = genes.
myfile <- file.path(rdatdir, paste0("3_", self, "_cleanDat.RData"))
data <- readRDS(myfile)
colNames <- rownames(data)
data <- t(data)
colnames(data) <- colNames

# Load adjmatrix or interaction network.
myfile <- file.path(rdatdir, paste0("3_", self, "_Adjm.RData"))
adjm <- as.matrix(readRDS(myfile))
rownames(adjm) <- colnames(adjm)

# Load Leidenalg graph partitions from 2_la-clustering.
myfiles <- c("Cortex" = file.path(rdatdir,"147731383_Cortex_CPMVertexPartition_partitions.csv"),
	     "Cortex_MCL" = file.path(rdatdir,"3_Cortex_MCL_partitions.csv"),
	    "Striatum" = file.path(rdatdir,"148436673_Striatum_CPMVertexPartition_partitions.csv"),
	    "Striatum_MCL" = file.path(rdatdir,"3_Striatum_MCL_partitioncs.csv"))
partitions <- data.table::fread(myfiles[partition_file], header=TRUE,drop = 1)

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
# Networks (edges) should be positive -> abs()
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

# Serial execution.
results <- list()
t0 <- Sys.time()
for (i in seq(1,2)){
	results[[i]] <- module_self_preservation(partitions,resolution=1)
}
t1 <- Sys.time()
dt <- t1-t0
message(paste("Time to execute in series:",dt))

# Register some nodes to do work.
workers <- makeCluster(c(rep("localhost", nThreads)),type = "SOCK")
registerDoParallel(workers)

# Parallel execution with %dopar%:
t0 <- Sys.time()
results <- foreach(resolution=seq(1,2)) %dopar% {
	module_self_preservation(partitions,resolution)
}
t1 <- Sys.time()
dt <- t1-t0
print(dt)
message(paste("Time to execute in parallel:",dt))

# Close parallel connections.
suppressWarnings(stopCluster(workers))

quit()

# Save to Rdata.
output_name <- paste0(jobID, "_", partition_file, 
		  "_Module_Self_Preservation.RData")
saveRDS(results, file.path(rdatdir, output_name))
