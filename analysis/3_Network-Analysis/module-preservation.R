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
stats <- c(1,2,6,7) # Module statistics to use for permutation testing.
self <- "Cortex"
test <- "Striatum"
strength <- "strong" # Criterion for preservation: strong = ALL, weak = ANY sig stats.
data_files <- c(self="Cortex",test="Striatum") # Which data to use?
net_files <- c(self="Cortex",test="Striatum") # Which networks to test self preservation in?
adjm_files <- c(self="Cortex",test="Striatum") # Which correlation (adjm) network to use?
partition_files <- c(self="Cortex_Surprise",test="Striatum_Surprise") # Which partition file to use?
replace_negative <- "zero" # How should negative weights be handled?
min_size <- 5 # minimum allowable size for a module.
self_preservation <- FALSE
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
devtools::load_all()

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatdir <- file.path(root, "rdata")
funcdir <- file.path(root, "R")

# Load expression data. Transpose -> rows = samples; columns = genes.
myfiles <- lapply(data_files,function(x) {
			  file.path(rdatdir,paste0("3_",x,"_cleanDat.RData"))
})
data <- lapply(myfiles,readRDS)
data <- lapply(data,t)

# Load adjmatrix.
myfiles <- lapply(adjm_files,function(x) {
			  file.path(rdatdir,paste0("3_",x,"_Adjm.RData"))
})
adjm <- lapply(myfiles,readRDS)
adjm <- lapply(adjm,as.matrix)
adjm <- lapply(adjm,function(x){ rownames(x) <- colnames(x); return(x) })

# Load network.
myfiles <- lapply(net_files,function(x) {
			  file.path(rdatdir,paste0("3_",x,"_Adjm.RData"))
})
net <- lapply(myfiles,readRDS)
net <- lapply(net,as.matrix)
net <- lapply(net,function(x){ rownames(x) <- colnames(x); return(x) })

# Load Leidenalg graph partitions from 2_la-clustering.
all_partitions <- c("Cortex" = "147731383_Cortex_CPMVertexPartition_partitions.csv",
		    "Cortex_MCL" = "3_Cortex_MCL_partitions.csv",
		     "Striatum" = "148436673_Striatum_CPMVertexPartition_partitions.csv",
		     "Striatum_MCL" = "3_Striatum_MCL_partitions.csv",
		     "Cortex_Surprise" = "3_Cortex_SurpriseVertexPartition_partitions.csv",
		     "Striatum_Surprise" = "3_Striatum_SurpriseVertexPartition_partitions.csv")
myfiles <- all_partitions[sapply(partition_files,function(x) grep(x,all_partitions))]
partitions <- lapply(myfiles,function(x) fread(file.path(rdatdir,x),drop=1))

# Coerce partitions dataframes into vectors.
partitiondf2vector <- function(parts) {
	parts <- split(parts,seq(nrow(parts)))
	parts <- lapply(parts,function(x) {
				# Coerce to numeric vector.
				y <- as.numeric(x)
				# Add 1 so non-zero index.
				if (min(y)==0) { y <- y + 1 }
				names(y) <- names(x)
				# Reset module index.
				y <- reset_index(y)
				# Remove small modules.
				too_small <- which(table(y) < min_size)
				y[y %in% too_small] <- 0 
				return(y)
			     })
	return(parts)
}
partitions <- lapply(partitions,partitiondf2vector)

# Number of resolutions.
resolutions <- min(unlist(lapply(partitions,length)))

#---------------------------------------------------------------------
## Permutation testing.
#---------------------------------------------------------------------

# Input for NetRep:
# 1. Expression data.
# 2. Correlation matrix or interaction network (co-expr, GO, or PPI).
# 3. Interaction network - Edges must be positive! Should network be weighted?
# 4. Network partitions.
data_list <- data
correlation_list <- net

# Networks (edges) should be positive...
if (replace_negative == "absolute value") {
	# Replace negative edges as absolute value.
	network_list <- lapply(net,abs)
} else if (replace_negative == "zero") {
	network_list <- lapply(net,function(x) { x[x<0] <- 0; return(x) })
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
for (resolution in resolutions) {
  message(paste("Working on partition", resolution, 
		"of", length(resolutions), "..."))
# Get list of modules.
parts <- lapply(partitions,"[[",resolution)
module_list <- list("self" = parts[[grep(self,names(parts))]],
		    "test" = parts[[grep(test,names(parts))]])
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
      test = "test",
      selfPreservation = self_preservation,
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
    output_name <- paste0(jobID, "_", self, "_", test,
			  "_Module_Self_Preservation.RData")
    saveRDS(results, file.path(rdatdir, output_name))
    message("Done!")
  }
} # Ends loop.
