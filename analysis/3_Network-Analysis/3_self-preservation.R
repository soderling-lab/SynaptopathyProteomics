#!/usr/bin/env Rscript
# Self-preservation of modules identified in WT and KO graphs are examined by
# permutation testing. The bicor co-expression matrices are transformed into 
# signed toplographical overlap matrices so that all edges are positive.

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change:
stats <- c(2,3,5,6,7)  # dont use average wedge weight and 
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

# Compute TOM adjcacency matrices--this insures that all edges are positve.
#wtTOM <- TOMsimilarity(wtAdjm, TOMType = "signed", verbose = 0)
#koTOM <- TOMsimilarity(koAdjm, TOMType = "signed", verbose = 0)
#rownames(wtTOM) <- colnames(wtTOM) <- colnames(wtAdjm)
#rownames(koTOM) <- colnames(koTOM) <- colnames(koAdjm)

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
network_list <- list(wt = wtAdjm, ko = koAdjm) # Use adjmatrices!

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

  check_modules <- function(x) {
    # Collect observed values, nulls, and p.values -> p.adj.
    obs <- x$observed[, stats]
    nulls <- apply(x$nulls, 2, function(x) apply(x, 1, mean))[, stats]
    q <- apply(x$p.values, 2, function(x) p.adjust(x, "bonferroni"))[, stats]
    q[is.na(q)] <- 1
    # If testing more than one statistic.
    fx <- c("strong"="all","weak"="any")[strength]
    if (length(stats) > 1) {
      sig <- apply(q < 0.05, 1, eval(fx))
      greater <- apply(obs > nulls, 1, eval(fx))
      less <- apply(obs < nulls, 1, eval(fx))
    } else {
      # If testing a single statistic.
      sig <- q < 0.05
      greater <- obs > nulls
      less <- obs < nulls
    }
    # Preserved, divergent, and ns modules.
    n <- length(x$nVarsPresent)
    v <- rep("ns", n)
    v[greater & sig] <- "preserved"
    v[less & sig] <- "divergent"
    names(v) <- names(x$nVarsPresent)
    return(v)
  } # ENDS function
  # Set ns modules to 0.
  v <- lapply(selfPreservation,check_modules) 
  out <- lapply(v, function(x) names(x)[x == "ns"])
  wtPartition[wtPartition %in% out[[1]]] <- 0
  koPartition[koPartition %in% out[[2]]] <- 0
  results[[i]] <- list(wt = wtPartition, ko = koPartition)
} # END LOOP.

# Save to Rdata.
output_name <- paste0(job, "_Module_Self_Preservation.RDS")
saveRDS(results, file.path(datadir, output_name))
