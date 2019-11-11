#!/usr/bin/env Rscript
# Loop through all resolutions of WT and KO graphs, comparing modules with
# permutation test in order to identify perserved and divergent modules.

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(NetRep)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Load expression data.
wtDat <- t(readRDS(list.files(rdatdir, pattern = "WT_cleanDat", full.names = TRUE)))
koDat <- t(readRDS(list.files(rdatdir, pattern = "KO_cleanDat", full.names = TRUE)))

# Load adjmatrices.
wtAdjm <- t(readRDS(list.files(rdatdir, pattern = "WT_Adjm.RData", full.names = TRUE)))
koAdjm <- t(readRDS(list.files(rdatdir, pattern = "KO_Adjm.RData", full.names = TRUE)))

# Load network partitions.
myfiles <- list.files(rdatdir, pattern = "*.partitions.csv", full.names = TRUE)
koParts <- data.table::fread(myfiles[1], drop = 1, skip = 1)
wtParts <- data.table::fread(myfiles[2], drop = 1, skip = 1)
colnames(koParts) <- colnames(wtParts) <- colnames(wtAdjm)

# Use a loop to make a list of partitions.
# Enforce the same format as other partitions.RData object.
# Add one such that all module indices are non-zero.
partitions <- list()
for (i in 1:nrow(wtParts)) {
  partitions[[i]] <- list(
    "wt" = unlist(wtParts[i, ]) + 1,
    "ko" = unlist(koParts[i, ]) + 1
  )
}

# LOOP TO ANALYZE ALL RESOLUTIONS:
output <- list()
nres <- length(partitions)
strength <- 1 # c(strong, weak)

for (r in 1:nres) {
  # Status report.
  message(paste("Working on resolution:", r, "..."))
  r_wt <- r_ko <- r
  # Extract from list.
  wtPartition <- partitions[[r_wt]][["wt"]]
  koPartition <- partitions[[r_ko]][["ko"]]
  # Split into modules.
  wtModules <- split(wtPartition, wtPartition)
  koModules <- split(koPartition, koPartition)
  # Total number of modules.
  nModules <- c("wt" = length(wtModules), "ko" = length(koModules))
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
  data_list <- list(wt = wtDat, ko = koDat)
  correlation_list <- list(wt = wtAdjm, ko = koAdjm)
  network_list <- list(wt = wtAdjm, ko = koAdjm)
  module_list <- list(wt = koPartition, ko = wtPartition)
  # ^This is correct: given the WT data/graph and the KO modules,
  # are modules preserved (the same) or divergent (different) in the KO graph?
  # Hypothesis for self-preservation.
  h0 <- list(
    wt = c(discovery = "wt", test = "wt"),
    ko = c(discovery = "ko", test = "ko")
  )
  # Perform permutation testing.
  # Suppress warnings which arise from small modules; NA p.vals are replaced with 1. 
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
        nThreads = 48,
        # nPerm = 100000,  # determined by the function.
        null = "overlap",
        alternative = "two.sided", # c(greater,less,two.sided)
        simplify = TRUE,
        verbose = FALSE
      )
    })
  })
  ## Identify preserved and divergent modules.
  check_modules <- function(x) {
    ## Strong preservation. All stats.
    all_obs <- x$observed
    all_nulls <- apply(x$nulls, 2, function(x) apply(x,1,mean))
    pmax <- apply(x$p.values,1,function(x) max(x,na.rm=TRUE))
    qmax <- p.adjust(pmax, "bonferroni")
    qmax[is.na(qmax)] <- 1
    n <- length(x$nVarsPresent)
    vmax <- rep("ns", n)
    # PRESERVED MODULES = obs > NULL & q < 0.05
   vmax[apply(all_obs > all_nulls,1,all) & qmax < 0.05] <- "preserved"
    # DIVERGENT MODULES = obs < NULL & q < 0.05
   vmax[apply(all_obs < all_nulls,1,all) & qmax < 0.05] <- "divergent"
   ## Just average edge weight = 1.
    obs <- x$observed[, 1]
    nullx <- apply(x$nulls[, 1, ], 1, mean) 
    p <- x$p.values[, 1]
    q <- p.adjust(p, "bonferroni")
    q[is.na(q)] <- 1
    n <- length(x$nVarsPresent)
    v <- rep("ns", n)
    # PRESERVED MODULES = obs > NULL & q < 0.05
    v[obs > nullx & q < 0.05] <- "preserved"
    # DIVERGENT MODULES = obs < NULL & q < 0.05
    v[obs < nullx & q < 0.05] <- "divergent"
    return(list("weak"=v,"strong"=vmax))
  } # ENDS function
  # Identify strong or weak changes...
  module_changes <- lapply(preservation, check_modules)
  module_changes <- sapply(module_changes,"[",strenght)
  names(module_changes) <- c("wt","ko")
  # Status report.
  message(paste("... Total number of WT modules:", nModules["wt"]))
  message(paste(
    "... ... Number of WT modules preserved in KO graph:",
    sum(module_changes$wt == "preserved")
  ))
  message(paste(
    "... ... Number of WT modules divergent in KO graph:",
    sum(module_changes$wt == "divergent")
  ))
  message(paste("... Total number of KO modules:", nModules["ko"]))
  message(paste(
    "... ... Number of KO modules preserved in WT graph:",
    sum(module_changes$ko == "preserved")
  ))
  message(paste(
    "... ... Number of KO modules divergent in WT graph:",
    sum(module_changes$ko == "divergent")
  ))
  # Return resolution, total number of modules, and module changes.
  output[[r]] <- list("resolution" = r, "nModules" = nModules, "Changes" = module_changes)
} # ENDS LOOP.

# Save output to file.
myfile <- file.path(rdatdir, "Network_Comparisons.RData")
saveRDS(output, myfile)
