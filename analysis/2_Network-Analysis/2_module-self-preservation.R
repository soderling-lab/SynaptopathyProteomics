#!/usr/bin/env Rscript

#' ---
#' title: Module Self-Preservation
#' description: Evaluate Module self-preservation by permutation testing.
#' authors: Tyler W Bradshaw
#' ---

## User parameters:
analysis_type = "Cortex"
root = "/mnt/d/projects/SynaptopathyProteomics"

## Permutation test options:
min_size = 5
verbose = FALSE
strength = "strong" 
stats = c(1, 2, 6, 7)
n_permutations = NULL
replace_negative = "zero"
self = test = analysis_type
n_threads = parallel::detectCores() - 1

## Input data in root/rdata:
input_data <- list("Cortex" = list(
				   adjm = "Cortex_Adjm.csv",
				   netw = "Cortex_NE_Adjm.csv",
			   	   data = "Cortex_norm_protein.csv",
				   part = "Cortex_NE_SurpriseVertexPartition.csv"),
		   "Striatum" = list(
				     adjm = "Striatum_Adjm.csv",
				     netw = "Striatum_NE_Adjm.csv",
				     data = "Striatum_norm_protein.csv",
				     part = "Striatum_NE_SurpriseVertexPartition.csv")
		   )[[analysis_type]]

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

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# Load renv.
renv::load(root)

# Global options and imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(NetRep)
	library(data.table)
})

# Additional functions in root/R.
TBmiscr::load_all()

# Directories.
root <- getrd()
rdatdir <- file.path(root, "rdata")

# Load expression data:
# Load the data, subset, coerce to matrix, Log2 transform, and 
# finally transpose such that rows = samples and columns = proteins.
myfile <- file.path(rdatdir, input_data[['data']])
dm <- fread(myfile) %>%
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
resolutions <- nrow(part_dt)

# Enforce consistent dimensions between data and adjm.
idx <- idy <- match(colnames(dm), colnames(adjm))
adjm <- adjm[idx, idy]
check <- all(colnames(dm) == colnames(adjm))
if (!check) {
  message("Problem: data doesn't match correlation matrix!")
}

# Enforce consistent dimensions between data and partition.
idy <- match(colnames(dm), colnames(part_dt))
partitions <- as.data.frame(part_dt)[, idy]
check <- all(colnames(dm) == colnames(partitions))
if (!check) {
  message("Problem: data doesn't match partitions!")
}

# Enforce consistent dimensions between data and network.
idx <- idy <- match(colnames(dm), colnames(netw))
netw <- netw[idx, idy]
check <- all(colnames(dm) == colnames(netw))
if (!check) {
  message("Problem: data doesn't match network!")
}

#-------------------------------------------------------------------------------
## Permutation testing.
#-------------------------------------------------------------------------------

# Input for NetRep:
# 1. Expression data.
# 2. Correlation matrix or interaction network (co-expr, GO, or PPI).
# 3. Interaction network - Edges must be positive! Should network be weighted?
# 4. Network partitions.
data_list <- list(self = dm)
correlation_list <- list(self = adjm)

# Networks (edges) should be positive...
if (replace_negative == "absolute value") {
  # Replace negative edges as absolute value.
  network_list <- list(self = abs(netw))
} else if (replace_negative == "zero") {
  netw[netw < 0] <- 0
  network_list <- list(self = netw)
}

# Module preservation stats to be used for enforcing preservation.
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
      nThreads = n_threads,
      nPerm = n_permutations, # Determined automatically by the function.
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
	  myfile <- file.path(rdatdir,paste0(self,"_partition_self_preservation_enforced.csv"))
	  do.call(rbind,results) %>% as.data.frame() %>% fwrite(myfile)
	  message("Done!")
  }
} # Ends loop.
