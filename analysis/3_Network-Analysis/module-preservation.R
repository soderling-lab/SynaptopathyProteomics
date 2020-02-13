#!/usr/bin/env Rscript

#' ---
#' title: Module Self-Preservation
#' description: Evaluate module preservation by permutation testing.
#' authors: Tyler W Bradshaw
#' ---

## User parameters to change:
stats <- c(1,2,6,7) 
strength <- "strong" # Criterion for preservation.
replace_negative <- "zero" # How to handle negative edge weights.

## Input data:
netw_files <- c(self="3_Cortex_Adjm.RData",
	        test="3_Striatum_Adjm.RData") 
data_files <- c(self="3_Cortex_cleanDat.RData",
		test="3_Striatum_cleanDat.RData") 
corr_files <- c(self="3_Cortex_Adjm.RData",
	        test="3_Striatum_Adjm.RData")
part_files <- c(self="2020-02-10_Cortex_Surprise_Module_Self_Preservation.RData",
		test="2020-02-10_Striatum_Surprise_Module_Self_Preservation.RData")

## Other NetRep parameters:
self <- "Cortex"
test <- "Striatum"
backgroundLabel <- 0
self_preservation <- FALSE
null <- "overlap"
alternative <- "greater"
verbose <- FALSE

## NetRep Permutation Statistics:
# 1. avg.weight (average edge weight) - Calculated from network. Assumes edge
#    weights are positive.
# 2. coherence - Calculated from the data. Quantifies the percent variance
#    explained by a modules summary vector.
# 3. cor.cor (concordance of correlation structure) - Calculated from
#    correlation matrix. Sensitive to small changes in correlation coefficients.
#    DONT USE -- sensitve to small changes in edge weight --false positives.
# 4. cor.degree (concordance of weighted degree) - Calculated from network. 
#    Assumes edge weights are positive. Sensitive to small changes in 
#    weighted degree. DONT USE. Sensitive to small changes -- false positives.
# 5. cor.contrib (concordance of node contribution) - Calculated from the data.
#    Sensitve to small changes in node contribution. DONT USE! Same as 3,4.
# 6. avg.cor (density of correlation structure) - Calculated from correlation
#    matrix.
# 7. avg.contrib (average node contribution) - Quantifies how similar nodes are
#    to summary profile.

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

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

#---------------------------------------------------------------------
## Collect input for permutation testing.
#---------------------------------------------------------------------

# Input for NetRep:
# 1. Expression data.
data <- lapply(data_files,function(x) readRDS(file.path(rdatdir,x)))
data <- lapply(data,t)
data_list <- data
names(data_list) <- c(self,test)

# 2. Correlation matrix or interaction network (co-expr, GO, or PPI).
adjm <- lapply(corr_files,function(x) readRDS(file.path(rdatdir,x)))
adjm <- lapply(adjm,as.matrix)
adjm <- lapply(adjm,function(x) { rownames(x) <- colnames(x) ; return(x) })
correlation_list <- adjm
names(correlation_list) <- c(self,test)

# 3. Interaction network - Edges must be positive! Should network be weighted?
net <- lapply(netw_files,function(x) readRDS(file.path(rdatdir,x)))
net <- lapply(net,as.matrix)
net <- lapply(net,function(x) { rownames(x) <- colnames(x) ; return(x) })
# Network (edges) should be positive.
# Replace negative edges.
if (replace_negative == "absolute value") {
	# Replace negative edges as absolute value.
	net <- lapply(net,abs)
} else if (replace_negative == "zero") {
	network_list <- lapply(net,function(x) { x[x<0] <- 0; return(x) })
}
names(network_list) <- c(self,test)

# 4. Network partitions.
parts <- lapply(part_files,function(x) readRDS(file.path(rdatdir,x)))
parts <- lapply(parts,"[[",1)
parts <- lapply(parts,reset_index)
module_list <- lapply(parts,reset_index)
names(module_list) <- c(self,test)

#---------------------------------------------------------------------
## Permutation testing.
#---------------------------------------------------------------------

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

# Perform permutation test for module self-preservation.
# Suppress warnings about missing values.
suppressWarnings({
selfPreservation <- NetRep::modulePreservation(
      network = network_list,
      data = data_list,
      correlation = correlation_list,
      moduleAssignments = module_list,
      modules = NULL,
      backgroundLabel = backgroundLabel,
      discovery = c(self,test),
      test = c(test,self),
      selfPreservation = self_preservation,
      nThreads = nThreads,
      nPerm = NULL, # Determined automatically by the function.
      null = null,
      alternative = alternative, # Greater for self-preservation.
      simplify = TRUE,
      verbose = verbose
      )
})

# Remove NS modules--set NS modules to 0.
check_preservation <- function(selfPreservation,self,test){
	partition <- module_list[[self]]
	results <- selfPreservation[[self]]
	preservedParts <- check_modules(results, strength, stats)
	nModules <- length(unique(partition[which(partition != 0)]))
	out <- names(preservedParts)[preservedParts == "ns"]
	partition[partition %in% out] <- 0
	nPreserved <- nModules - length(out)
	message(paste("...", nPreserved, "of", nModules, self,
		      "modules are preserved in the",test,"network."))
	return(partition)
}

# Check how many modules are preserved.
results <- list(self = check_preservation(selfPreservation,self="Cortex",test="Striatum"),
		test = check_preservation(selfPreservation,self="Striatum",test="Cortex"))

# Save to Rdata.
output_name <- paste0(jobID, "_", self, "_", test,
		      "_Module_Self_Preservation.RData")
saveRDS(results, file.path(rdatdir, output_name))
message("Done!")
