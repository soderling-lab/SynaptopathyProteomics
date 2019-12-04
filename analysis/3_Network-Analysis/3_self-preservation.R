#!/usr/bin/env Rscript

#' ---
#' title:
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
stats = c(1:7)          # Module statistics to use for permutation testing.
strength = "strong"     # Criterion for preservation: strong = ALL, weak = ANY sig stats.
weighted = FALSE        # Weighted or unweighted. If TRUE, then appropriate soft-power will be calculated.
self = "combined"       # Which networks to test self preservation in? #self = c("wt","ko")     
nres <- 100             # Total number of resolutions to be anlyzed.

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
myfun <- list.files(funcdir, pattern = "silently.R", full.names = TRUE)
invisible(sapply(myfun, source))

# Load expression data. Transpose -> rows = samples; columns = genes.
wtDat <- t(readRDS(file.path(datadir, "3_WT_cleanDat.RData")))
koDat <- t(readRDS(file.path(datadir, "3_KO_cleanDat.RData")))
combDat <- t(readRDS(file.path(datadir, "3_Combined_cleanDat.RData")))

# Compute adjmatrix:
wtAdjm <- silently(WGCNA::bicor(wtDat))
koAdjm <- silently(WGCNA::bicor(koDat))
combAdjm <- silently(WGCNA::bicor(combDat))

# Weighted or unweighted?
if (weighted) {
	# Calculate power for approximate scale free fit.
	sft <- silently({
  sapply(list(wtDat, koDat,combDat), function(x) {
    pickSoftThreshold(x,
      corFnc = "bicor",
      networkType = "signed",
      RsquaredCut = 0.8
    )$powerEstimate
  })
})
	names(sft) <- c("wt", "ko","combined")
} else {
	# Power = 1 == Unweighted
	sft <- rep(1,3)
	names(sft) <- c("wt", "ko","combined")
}

# Load Leidenalg graph partitions.
myfiles <- list.files(datadir, pattern = "*partitions.csv", full.names = TRUE)
koParts <- data.table::fread(myfiles[grep("KO",myfiles)], drop = 1, skip = 1)
wtParts <- data.table::fread(myfiles[grep("WT",myfiles)], drop = 1, skip = 1)
#combParts <- data.table::fread(myfiles[grep("Combined",myfiles)], drop = 1, skip = 1)
colnames(koParts) <- colnames(wtParts) <- colnames(wtAdjm)
#colnames(combParts) <- colnames(combAdjm)

#-------------------------------------------------------------------------------
## Permutation testing.
#-------------------------------------------------------------------------------

# Input for NetRep:
# Networks (edges) should be positive -> AbsoluteValue()
data_list <- list(wt = wtDat, ko = koDat, combined = combDat)
correlation_list <- list(wt = wtAdjm, ko = koAdjm, combined = combDat)
network_list <- list(wt = abs(wtAdjm^sft["wt"]), 
		     ko = abs(koAdjm^sft["ko"]),
		     combined = abs(combAdjm^sft["combined"]))


# Loop through partitions, evaluating self-preservation.
results <- list()
for (i in 1:nres) {
  # Status report.
  message(paste("working on partition", i, "..."))
  # Get partition--adding 1 so that all module assignments >0.
  wtPartition <- as.integer(wtParts[i, ]) + 1
  koPartition <- as.integer(koParts[i, ]) + 1
  combPartition <- as.integer(combParts[i, ]) + 1
  names(wtPartition) <- names(koPartition) <- colnames(wtAdjm)
  names(combPartition) <- colnames(combAdjm)
  module_list <- list(wt = wtPartition, ko = koPartition, combined = combPartition)
  # Perform permutation test for module self-preservation.
  H0 <- as.list(self)
  suppressWarnings({
    selfPreservation <- lapply(H0, function(x) {
      NetRep::modulePreservation(
        network = network_list,
        data = data_list,
        correlation = correlation_list,
        moduleAssignments = module_list,
        modules = NULL,
        backgroundLabel = 0,
        discovery = x,
        test = x,
        selfPreservation = TRUE,
        nThreads = nThreads,
        # nPerm = 100000, # Determined automatically by the function.
        null = "overlap",
        alternative = "greater", # Greater for self-preservation.
        simplify = TRUE,
        verbose = FALSE
      )
    })
  })
  # Function to check module preservation/divergence.
  check_modules <- function(x) {
    # Collect observed values, nulls, and p.values -> p.adj.
    obs <- x$observed[, stats]
    nulls <- apply(x$nulls, 2, function(x) apply(x, 1, mean))[, stats]
    q <- apply(x$p.values, 2, function(x) p.adjust(x, "bonferroni"))[, stats]
    q[is.na(q)] <- 1
    # If testing more than one statistic, consider strong or weak preservation.
    fx <- c("strong" = "all", "weak" = "any")[strength]
    if (length(stats) > 1) {
      sig <- apply(q < 0.05, 1, eval(fx))
      greater <- apply(obs > nulls, 1, eval(fx))
      less <- apply(obs < nulls, 1, eval(fx))
    } else {
      # If testing a single statistic...
      sig <- q < 0.05
      greater <- obs > nulls
      less <- obs < nulls
    }
    # Define preserved, divergent, and ns modules.
    nModules <- length(x$nVarsPresent)
    v <- rep("ns", nModules)
    v[greater & sig] <- "preserved"
    v[less & sig] <- "divergent"
    names(v) <- names(x$nVarsPresent)
    return(v)
  }
  # Remove NS modules--set NS modules to 0.
  preservedPartitions <- lapply(selfPreservation, check_modules)
  out <- lapply(preservedPartitions, function(x) names(x)[x == "ns"])
  wtPartition[wtPartition %in% out[[1]]] <- 0
  koPartition[koPartition %in% out[[2]]] <- 0
  # Return results.
  results[[i]] <- list(wt = wtPartition, ko = koPartition)
} # END LOOP.

# Save to Rdata.
output_name <- paste0(jobID, "_Module_Self_Preservation.RDS")
saveRDS(results, file.path(datadir, output_name))
