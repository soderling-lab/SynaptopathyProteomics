#!/usr/bin/env Rscript

#' ---
#' title: Module Self-Preservation
#' description: Evaluate Module self-preservation by permutation testing.
#' authors: Tyler W Bradshaw
#' ---

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# User parameters to change:
stats <- c(1, 2, 6, 7) # Module statistics to use for permutation testing.
self <- "Striatum"
test <- "PPI"
strength <- "strong" # Criterion for preservation: strong = ALL, weak = ANY sig stats.
data_file <- "Striatum" # Which data to use?
net_file <- "PPI" # Which networks to test self preservation in?
adjm_file <- "Striatum" # Which correlation (adjm) network to use?
partition_file <- "Striatum_Surprise" # Which partition file to use?
replace_negative <- "zero" # How should negative weights be handled?
min_size <- 5 # minimum allowable size for a module.
verbose <- FALSE

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
  nThreads <- parallel::detectCores() - 1
  jobID <- Sys.Date()
}

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(NetRep)
})

# Additional functions.
suppressWarnings({
  devtools::load_all()
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatdir <- file.path(root, "rdata")
funcdir <- file.path(root, "R")

# Load expression data. Transpose -> rows = samples; columns = genes.
myfile <- file.path(rdatdir, paste0("3_", data_file, "_cleanDat.RData"))
data <- readRDS(myfile)
colNames <- rownames(data)
data <- t(data)
colnames(data) <- colNames

# Load adjmatrix.
myfile <- file.path(rdatdir, paste0("3_", adjm_file, "_Adjm.RData"))
adjm <- as.matrix(readRDS(myfile))
rownames(adjm) <- colnames(adjm)

# Load network.
myfile <- file.path(rdatdir, paste0("3_", net_file, "_Adjm.RData"))
net <- as.matrix(readRDS(myfile))
rownames(net) <- colnames(net)

# Load Leidenalg graph partitions from 2_la-clustering.
myfiles <- c(
  "Cortex" = "147731383_Cortex_CPMVertexPartition_partitions.csv",
  "Cortex_MCL" = "3_Cortex_MCL_partitions.csv",
  "Striatum" = "148436673_Striatum_CPMVertexPartition_partitions.csv",
  "Striatum_MCL" = "3_Striatum_MCL_partitions.csv",
  "Cortex_Surprise" = "3_Cortex_SurpriseVertexPartition_partitions.csv",
  "Striatum_Surprise" = "3_Striatum_SurpriseVertexPartition_partitions.csv"
)
partitions <- fread(file.path(rdatdir, myfiles[partition_file]),
  header = TRUE, drop = 1
)
resolutions <- nrow(partitions)

# Check that all columns in the data are in adjm and network.
out1 <- colnames(data) %notin% colnames(adjm)
out2 <- colnames(data) %notin% colnames(net)
if (sum(out1) > 0 | sum(out2) > 0) {
  message("Warning: removing columns from data that are not in network.")
}

# Enforce consistent dimensions between data and adjm.
idx <- idy <- match(colnames(data), colnames(adjm))
adjm <- adjm[idx, idy]
check <- all(colnames(data) == colnames(adjm))
if (!check) {
  message("Problem: data doesn't match correlation matrix!")
}

# Enforce consistent dimensions between data and partitions.
idy <- match(colnames(data), colnames(partitions))
partitions <- as.data.frame(partitions)[, idy]
check <- all(colnames(data) == colnames(partitions))
if (!check) {
  message("Problem: data doesn't match partitions!")
}

# Enforce consistent dimensions between data and network.
idz <- match(colnames(data), colnames(net))
net <- net[idz, idz]
check <- all(colnames(data) == colnames(net))
if (!check) {
  message("Problem: data doesn't match network!")
}

# Final check.
check <- all(colnames(partitions) == colnames(data) & colnames(data) == colnames(net))

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
  net[net < 0] <- 0
  network_list <- list(self = net)
}

# Module preservation stats.
module_stats <- paste(c(
  "avg.weight", "coherence", "cor.cor", "cor.degree",
  "cor.contrib", "avg.cor", "avg.contrib"
)[stats], collapse = ", ")

# Status report:
message(paste(
  "Evaluating self-preservation of",
  self, "modules in the", test, "network."
))
message(paste0(
  "Module statistic(s) used to evaluate module preservation: ",
  module_stats
), ".")
message(paste0(
  "Criterion for module preservation/divergence: ",
  strength, ".", "\n"
))

# Loop through partitions, evaluating self-preservation.
results <- list()
for (resolution in resolutions) {
  message(paste(
    "Working on partition", resolution,
    "of", length(resolutions), "..."
  ))
  # Get partition.
  partition <- as.integer(partitions[resolution, ])
  # Add 1 so that all module assignments >0.
  if (min(partition) == 0) {
    partition <- partition + 1
  }
  partition <- reset_index(partition)
  names(partition) <- colnames(partitions)
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
  if (resolution == length(resolutions)) {
    output_name <- paste0(
      jobID, "_", partition_file, "_", net_file,
      "_Module_Self_Preservation.RData"
    )
    saveRDS(results, file.path(rdatdir, output_name))
    message("Done!")
  }
} # Ends loop.
