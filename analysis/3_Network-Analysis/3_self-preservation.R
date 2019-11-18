#!/usr/bin/env Rscript

#' ---
#' title:
#' description: Evaluate WT and KO module self-preservation with perm test.
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# Module statistics:
# 1. avg.weight - (average edge weight) assumes positive edges - calculated 
#    from network.
# 2. coherence (module coherence) - calculated from modules summary profile.
# 3. cor.cor (concordance of correlation structure) - calculated from
#    correlation matrix.
# 4. cor.degree - assumes positive edges - calculated from network.
# 5. cor.contrib (concordance of node contribution) - calculated from modules 
#    summary profile. 
# 6. avg.cor (density of correlation structure) - calculate from correlation 
#    matrix.
# 7. avg.contrib (average node contribution) - calculated from modules summary
#    profile.

# User parameters to change:
stats <- c(1:7)  
strength <- "strong"   # Preservation criterion strong = all, or weak = any sig stats.

# Is this a slurm job?
slurm <- any(grepl("SLURM", names(Sys.getenv())))
if (slurm) {
  # SLURM job notes - sent to job_*.info
  nThreads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # Number of threads.
  job <- as.integer(Sys.getenv("SLURM_JOBID"))
  info <- as.matrix(Sys.getenv())
  idx <- grepl("SLURM", rownames(info))
  myfile <- file.path("./out", paste0("job_", job, ".info"))
  write.table(info[idx, ], myfile, col.names = FALSE, quote = FALSE, sep = "\t")
} else {
  nThreads <- 8
  job <- ""
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

# Compute adjmatrix:
wtAdjm <- silently(WGCNA::bicor(wtDat))
koAdjm <- silently(WGCNA::bicor(koDat))

# Calculate power for approximate scale free fit.
sft <- silently({
	sapply(list(wtDat,koDat),function(x) 
	       pickSoftThreshold(x,
				 corFnc="bicor",
				 networkType="signed",
				 RsquaredCut=0.8)$powerEstimate)
})
names(sft) <- c("wt","ko")

# Load partitions.
myfiles <- list.files(datadir, pattern = "*partitions.csv", full.names = TRUE)
koParts <- data.table::fread(myfiles[1], drop = 1, skip = 1)
wtParts <- data.table::fread(myfiles[2], drop = 1, skip = 1)
colnames(koParts) <- colnames(wtParts) <- colnames(wtAdjm)

#-------------------------------------------------------------------------------
## Permutation testing.
#-------------------------------------------------------------------------------

# Input for NetRep:
# Networks should be transformed with power for scale free fit, and positive.
data_list <- list(wt = wtDat, ko = koDat)
correlation_list <- list(wt = wtAdjm, ko = koAdjm)
network_list <- list(wt = abs(wtAdjm^sft["wt"]), ko = abs(koAdjm^sft["ko"])) 

# Loop through partitions, evaluating self-preservation.
nres <- dim(koParts)[1]
results <- list()

for (i in 1:nres) {
  # Status report.
  message(paste("working on partition", i, "..."))
  # Get partition--adding 1 so that all module assignments >0.
  wtPartition <- as.integer(wtParts[i, ]) + 1
  koPartition <- as.integer(koParts[i, ]) + 1
  names(wtPartition) <- names(koPartition) <- colnames(wtAdjm)
  module_list <- list(wt = wtPartition, ko = koPartition)
  # Perform permutation test for module self-preservation.
  # Use lapply to do for both wt and ko networks...
  self <- as.list(c("wt", "ko"))
  suppressWarnings({
    selfPreservation <- lapply(self, function(x) {
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
        alternative = "greater",
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
    fx <- c("strong"="all","weak"="any")[strength]
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
  preservedPartitions <- lapply(selfPreservation,check_modules)
  out <- lapply(preservedPartitions, function(x) names(x)[x == "ns"])
  wtPartition[wtPartition %in% out[[1]]] <- 0
  koPartition[koPartition %in% out[[2]]] <- 0
  # Return results.
  results[[i]] <- list(wt = wtPartition, ko = koPartition)
} # END LOOP.

# Save to Rdata.
output_name <- paste0(job, "_Module_Self_Preservation.RDS")
saveRDS(results, file.path(datadir, output_name))
