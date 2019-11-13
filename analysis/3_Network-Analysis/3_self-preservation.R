#!/usr/bin/env Rscript
# Examine module self-preservation.

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

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
} else  {
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

# Compute TOM adjcacency matrices--this insures that all edges are positve.
wtTOM <- TOMsimilarity(wtAdjm,TOMType="signed",verbose=0)
koTOM <- TOMsimilarity(koAdjm,TOMType="signed",verbose=0)
rownames(wtTOM) <- colnames(wtTOM) <- colnames(wtAdjm)
rownames(koTOM) <- colnames(koTOM) <- colnames(koAdjm)

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
network_list <- list(wt = wtTOM, ko = koTOM) # Use adjmatrices!

# Loop through partitions, evaluating self-preservation.
n <- dim(koParts)[1]
results <- list()

for (i in 1:100) {
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
      # nPerm = 100000,
      null = "overlap",
      alternative = "greater",
      simplify = TRUE,
      verbose = FALSE
    )
  })
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
output_name <- paste0(job,"_self_preservation_results.RDS")
saveRDS(results, file.path(datadir,output_name))
