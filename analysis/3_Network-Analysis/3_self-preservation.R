#!/usr/bin/env Rscript
# Examine module self-preservation.

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
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

# Load partitions.
myfiles <- list.files(datadir, pattern = "*partitions.csv", full.names = TRUE)
koParts <- data.table::fread(myfiles[1], drop = 1, skip = 1)
wtParts <- data.table::fread(myfiles[2], drop = 1, skip = 1)
colnames(koParts) <- colnames(wtParts) <- colnames(wtAdjm)

#-------------------------------------------------------------------------------
## Permutation testing.
#-------------------------------------------------------------------------------

# Input for NetRep:
data_list <- list(wt = wtDat, ko = koDat)
correlation_list <- list(wt = wtAdjm, ko = koAdjm)
network_list <- list(wt = wtAdjm, ko = koAdjm)

# Loop through partitions, evaluating self-preservation.
n <- dim(koParts)[1]
results <- list()

for (i in 1:n) {
  # status
  message(paste("working on partition", i, "..."))
  # Get partition.
  wtPartition <- as.integer(wtParts[i, ]) + 1
  koPartition <- as.integer(koParts[i, ]) + 1
  names(wtPartition) <- names(koPartition) <- colnames(wtAdjm)
  module_list <- list(wt = wtPartition, ko = koPartition)
  # Perform permutation test for module self-preservation.
  # Done for both wt and ko networks...
  self <- as.list(c("wt", "ko"))
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
      nThreads = 24,
      # nPerm = 100000,
      null = "overlap",
      alternative = "greater",
      simplify = TRUE,
      verbose = FALSE
    )
  })

  # Function to get max pvalue.
  maxp <- function(preservation) {
    p <- apply(preservation$p.values, 1, function(x) max(x, na.rm = TRUE))
    q <- p.adjust(p, "bonferroni")
    return(q)
  }
  q <- lapply(selfPreservation, maxp)
  # Modules with NS preservation stats.
  out <- lapply(q, function(x) names(x)[x > 0.05])
  # For NS modules, set module membership to 0.
  wtPartition[wtPartition %in% out[[1]]] <- 0
  koPartition[koPartition %in% out[[2]]] <- 0
  # Return results.
  results[[i]] <- list(wt = wtPartition, ko = koPartition)
} # END LOOP.

# Save to Rdata.
saveRDS(results, file.path(datadir, "3_self_preservation_results.RDS"))
