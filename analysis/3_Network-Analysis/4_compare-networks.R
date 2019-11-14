#!/usr/bin/env Rscript
# Loop through all resolutions of WT and KO graphs, comparing modules with
# permutation test in order to identify perserved and divergent modules.

# For example:
# Given a network partition at a given resolution:
# Get observed statistics for a module in the opposite (test) graph (e.g. KO module avg.edge weight in WT graph).
# Get nulls from 10,000 randomizations of the test graph.
# Compare observed versus NULL distributions.
# Correct p.values for n comparisons.
# If observed statistic is significantly greater than null -> preserved.
# If observed statistic is significantly less than null --> divergent.
# WT Modules that are divergent in the KO graph coorespond to LOF.
# KO Modules that are divergent in the WT graph coorespond to GOF.

# Criterion for preservation/divergence can be weak or strong.
# If weak, then requirement is ANY significant stats.
# If strong, then requirement is ALL significant stats.
# User can select which stats to use.

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change:
stats <- c(1:7) # Which of the seven module statistics to use.
strength <- "any" # Preservation criterion strong = all, or weak = any sig stats.
res <- 1 #c(1:100)      # Resolutions to analyze.

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
  library(data.table)
  library(dplyr)
  library(WGCNA)
  library(NetRep)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Load expression data.
wtDat <- t(readRDS(list.files(rdatdir,
  pattern = "WT_cleanDat",
  full.names = TRUE
)))
koDat <- t(readRDS(list.files(rdatdir,
  pattern = "KO_cleanDat",
  full.names = TRUE
)))

# Load adjmatrices.
wtAdjm <- t(readRDS(list.files(rdatdir,
  pattern = "WT_Adjm.RData",
  full.names = TRUE
)))
koAdjm <- t(readRDS(list.files(rdatdir,
  pattern = "KO_Adjm.RData",
  full.names = TRUE
)))

# Compute TOM adjcacency matrices--this insures that all edges are positve.
wtTOM <- TOMsimilarity(wtAdjm, TOMType = "signed", verbose = 0)
koTOM <- TOMsimilarity(koAdjm, TOMType = "signed", verbose = 0)
rownames(wtTOM) <- colnames(wtTOM) <- colnames(wtAdjm)
rownames(koTOM) <- colnames(koTOM) <- colnames(koAdjm)

# Load network partitions. Self-preservation enforced.
myfile <- list.files(rdatdir, pattern = "5716254", full.names = TRUE)
partitions <- readRDS(myfile)

#------------------------------------------------------------------------------
# Loop through all resolutions and perform permutation test.
#------------------------------------------------------------------------------

# All module statistics.
module_stats <- paste(c(
  "avg.weight", "coherence", "cor.cor", "cor.degree",
  "cor.contrib", "avg.cor", "avg.contrib"
)[stats], collapse = ", ")
# Status report:
nres <- length(res)
message(paste("Comparing WT and KO graphs at", nres, "resolution(s)."))
message(paste0(
  "Module statistic(s) used to evaluate module preservation: ",
  module_stats
), ".")
message(paste0(
  "Criterion for module preservation/divergence: ",
  strength, ".", "\n"
))

# LOOP TO ANALYZE ALL RESOLUTIONS:
output <- list()
for (r in res) {
  # Status report.
  message(paste("Working on resolution:", r, "..."))
  # Extract from list.
  wtPartition <- partitions[[r]][["wt"]]
  koPartition <- partitions[[r]][["ko"]]
  # Split into modules.
  wtModules <- split(wtPartition, wtPartition)
  koModules <- split(koPartition, koPartition)
  # Total number of modules.
  nModules <- c(
    "wt" = sum(names(table(wtPartition)) != 0),
    "ko" = sum(names(table(koPartition)) != 0)
  )
  # Checks:
  if (!all(colnames(wtDat) == colnames(koDat))) {
    stop("Input data don't match!")
  }
  if (!all(colnames(wtAdjm) == colnames(koAdjm))) {
    stop("Input data don't match!")
  }
  if (!all(names(wtPartition) %in% colnames(wtDat))) {
    stop("Input data don't match!")
  }
  if (!all(names(koPartition) %in% colnames(koDat))) {
    stop("Input data don't match!")
  }
  # Input for NetRep:
  # Note the networks are what are used to calc the avg edge weight statistic.
  # Note that NetRep assumes all edges are positive in calculating
  # avg.edge.weight and cor.degree.
  data_list <- list(wt = wtDat, ko = koDat)
  correlation_list <- list(wt = wtAdjm, ko = koAdjm)
  network_list <- list(wt = wtTOM, ko = koTOM)
  module_list <- list(wt = koPartition, ko = wtPartition)
  # ^This is correct: given the WT data/graph and the KO modules,
  # are modules preserved (the same) or divergent (different) in the KO graph?
  # Hypothesis for self-preservation.
  h0 <- list(
    wt = c(discovery = "wt", test = "wt"),
    ko = c(discovery = "ko", test = "ko")
  )
  # Perform permutation testing.
  # Suppress warnings which arise from small modules;
  # NA p.vals are replaced with 1.
  suppressWarnings({
    preservation <- lapply(h0, function(x) {
      NetRep::modulePreservation(
        network = network_list,
        data = data_list,
        correlation = correlation_list,
        moduleAssignments = module_list,
        modules = NULL,
        backgroundLabel = 0,
        discovery = x["discovery"],
        test = x["test"],
        selfPreservation = TRUE,
        nThreads = nThreads,
        # nPerm = 100000,  # determined by the function.
        null = "overlap",
        alternative = "two.sided", # c(greater,less,two.sided)
        simplify = TRUE,
        verbose = FALSE
      )
    })
  })
  # Identify preserved and divergent modules.
  check_modules <- function(x) {
    # Collect observed values, nulls, and p.values -> p.adj.
    obs <- x$observed[, stats]
    nulls <- apply(x$nulls, 2, function(x) apply(x, 1, mean))[, stats]
    q <- apply(x$p.values, 2, function(x) p.adjust(x, "bonferroni"))[, stats]
    q[is.na(q)] <- 1
    # If testing more than one statistic.
    if (length(stats) > 1) {
      fx <- c("strong" = "all", "weak" = "any")[strength]
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
    return(v)
  } # ENDS function
  # Collect strong or weak changes...
  module_changes <- lapply(preservation, check_modules)
  names(module_changes) <- c("ko", "wt") # fix names bc they were confusing.
  # Calculate percent NS, divergent, preserved.
  wtProts <- wtPartition
  wtProts[wtProts == "0"] <- "not-clustered"
  wtProts[wtProts %in% names(wtModules)[module_changes$wt == "preserved"]] <- "preserved"
  wtProts[wtProts %in% names(wtModules)[module_changes$wt == "divergent"]] <- "divergent"
  wtProts[wtProts %in% names(wtModules)[module_changes$wt == "ns"]] <- "ns"
  koProts <- koPartition
  koProts[koProts == "0"] <- "not-clustered"
  koProts[koProts %in% names(koModules)[module_changes$ko == "preserved"]] <- "preserved"
  koProts[koProts %in% names(koModules)[module_changes$ko == "divergent"]] <- "divergent"
  koProts[koProts %in% names(koModules)[module_changes$ko == "ns"]] <- "ns"
  pdWT <- sum(wtProts == "divergent") / length(wtProts)
  ppWT <- sum(wtProts == "preserved") / length(wtProts)
  pnWT <- sum(wtProts == "ns") / length(wtProts)
  pncWT <- sum(wtProts == "not-clustered") / length(wtProts)
  pdKO <- sum(koProts == "divergent") / length(koProts)
  ppKO <- sum(koProts == "preserved") / length(koProts)
  pnKO <- sum(koProts == "ns") / length(koProts)
  pncKO <- sum(koProts == "not-clustered") / length(koProts)
  # WT status report.
  message(paste("... Total number of WT modules:", nModules["wt"]))
  message(paste0(
    "... ... Number of WT modules preserved in KO graph: ",
    sum(module_changes$wt == "preserved"), " (", round(100 * ppWT, 3), "% proteins)."
  ))
  message(paste0(
    "... ... Number of WT modules divergent in KO graph: ",
    sum(module_changes$wt == "divergent"), " (", round(100 * pdWT, 3), "% proteins)."
  ))
  # KO status report.
  message(paste("... Total number of KO modules:", nModules["ko"]))
  message(paste0(
    "... ... Number of KO modules preserved in WT graph: ",
    sum(module_changes$ko == "preserved"), " (", round(100 * ppKO, 3), "% proteins)."
    ))
  message(paste0(
    "... ... Number of KO modules divergent in WT graph: ",
    sum(module_changes$ko == "divergent"), " (", round(100 * pdKO), "% proteins)."
  ))
  # Return.
  output[[r]] <- NULL
} # ENDS LOOP.

# Save output to file.
output_name <- paste0(job, "_Network_Comparisons.RData")
myfile <- file.path(rdatdir, output_name)
saveRDS(output, myfile)
