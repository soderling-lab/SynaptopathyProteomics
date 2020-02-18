#!/usr/bin/env Rscript

#' ---
#' title: Module Preservation
#' description: Evaluate module preservation by permutation testing.
#'              This is a more general form of the module-self-preservation
#'              script. Should be able to handle different permutation
#'              test setups.
#' authors: Tyler W Bradshaw
#' ---

## User parameters to change:
stats <- c(1,2,6,7) 
#stats = 1 # For PPI preservation, just consider edge weight.
strength = "strong" # Criterion for preservation.
replace_negative = "zero" # How to handle negative edge weights. abs or zero.

## Other NetRep parameters:
nPerm = NULL
verbose = FALSE
null = "overlap"
backgroundLabel = 0
alternative = "greater"

# Organization of the permutation test.
# The data in data_list will be used, see below.
# Define the appropriate self and test for all inputs.
perm_test <- list(discovery = "Striatum", 
		  test = "Cortex",
		  self_preservation = FALSE,
		  network =     c(self="Striatum",test="Cortex"),
		  data =        c(self="Striatum",test="Cortex"),
		  correlation = c(self="Striatum",test="Cortex"),
		  module =      c(self="Striatum",test="Cortex"))

## Input data.
# Data should be in rdata/.
data_files <- list(
		   network = list(Cortex = "3_Cortex_Adjm.RData",
			       Striatum = "3_Striatum_Adjm.RData",
		               PPI = "3_PPI_Adjm.RData"),
		   data = list(Cortex = "3_Cortex_cleanDat.RData",
			       Striatum = "3_Striatum_cleanDat.RData"),
		   correlation = list(Cortex = "3_Cortex_Adjm.RData",
		               Striatum = "3_Striatum_Adjm.RData"),
		   module = list(Cortex = "2020-02-10_Cortex_Surprise_Module_Self_Preservation.RData",
			       Striatum ="2020-02-10_Striatum_Surprise_Module_Self_Preservation.RData")
		   )

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

# Check if performing self-preservation test.
if (perm_test$self_preservation == TRUE) {
	self_preservation <- TRUE
	self <- perm_test$discovery
	test <- perm_test$discovery
} else {
	self_preservation <- FALSE
	self <- perm_test$discovery
	test <- perm_test$test
}

# Input for NetRep:
# 1. Expression data.
myfiles <- data_files$data[perm_test$data]
data <- lapply(myfiles,function(x) readRDS(file.path(rdatdir,x)))
data_list <- lapply(data,t)
if (self_preservation) { names(data_list) <- self }

# 2. Correlation matrix or interaction network (co-expr, GO, or PPI).
myfiles <- data_files$correlation[perm_test$correlation]
adjm <- lapply(myfiles,function(x) readRDS(file.path(rdatdir,x)))
adjm <- lapply(adjm,as.matrix)
adjm <- lapply(adjm,function(x) { rownames(x) <- colnames(x) ; return(x) })
correlation_list <- adjm
if (self_preservation) { names(correlation_list) <- self }

# 3. Interaction network - Edges must be positive!
myfiles <- data_files$network[perm_test$network]
net <- lapply(myfiles,function(x) readRDS(file.path(rdatdir,x)))
net <- lapply(net,as.matrix)
net <- lapply(net,function(x) { rownames(x) <- colnames(x) ; return(x) })
# Network (edges) should be positive.
# Replace negative edges.
if (replace_negative == "abs") {
	# Replace negative edges as absolute value.
	network_list <- lapply(net,abs)
} else if (replace_negative == "zero") {
	network_list <- lapply(net,function(x) { x[x<0] <- 0; return(x) })
}
if (self_preservation) { names(network_list) <- self }

# 4. Network partitions.
myfiles <- data_files$module[perm_test$module]
parts <- lapply(myfiles,function(x) readRDS(file.path(rdatdir,x)))
parts <- lapply(parts,"[[",1)
parts <- lapply(parts,reset_index)
module_list <- lapply(parts,reset_index)
if (self_preservation) { names(module_list) <- self }

# Insure matching order of input data.
colNames <- colnames(data_list[[1]])
data_list <- lapply(data_list,function(x) sortInput(x,colNames))
network_list <- lapply(network_list,function(x) sortInput(x,colNames))
correlation_list <- lapply(correlation_list,function(x) sortInput(x,colNames))
module_list <- lapply(module_list,function(x) sortInput(x,colNames))

#---------------------------------------------------------------------
## Permutation testing.
#---------------------------------------------------------------------

# Module preservation stats.
module_stats <- paste(c(
  "avg.weight", "coherence", "cor.cor", "cor.degree",
  "cor.contrib", "avg.cor", "avg.contrib"
)[stats], collapse = ", ")

# Status report:
message(paste0(
  "Module statistic(s) used to evaluate module preservation: ",
  module_stats), ".")
message(paste0(
  "Criterion for module preservation/divergence: ",
  strength, ".", "\n"
))
message(paste("Evaluating preservation of",perm_test$discovery, 
	      "modules in the", perm_test$test, "network..."))

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
      discovery = self,
      test = test,
      selfPreservation = self_preservation,
      nThreads = nThreads,
      nPerm = nPerm, # If null then determined automatically by the function.
      null = null,
      alternative = alternative, # Greater for self-preservation.
      simplify = TRUE,
      verbose = verbose
      )
})

# Remove NS modules--set NS modules to 0.
check_preservation <- function(selfPreservation,module_list,self,test,stats){
	# Wrapper around check_modules()
	partition <- module_list[[self]]
	preservedParts <- check_modules(selfPreservation, strength, stats)
	nModules <- length(unique(partition[which(partition != 0)]))
	out <- names(preservedParts)[preservedParts == "ns"]
	partition[partition %in% out] <- 0
	nPreserved <- nModules - length(out)
	message(paste("...", nPreserved, "of", nModules, self,
		      "modules are preserved in the",test,"network.\n"))
	# Return partition with NS modules set to 0.
	return(partition)
}

# Check how many modules are preserved.
results <- check_preservation(selfPreservation,module_list,self,test,stats)

# Save to Rdata.
output_name <- paste0(jobID, "_", perm_test$discovery, "_", 
		      perm_test$test, "_Module_Self_Preservation.RData")
saveRDS(results, file.path(rdatdir, output_name))
