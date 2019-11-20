#!/usr/bin/env Rscript

#' ---
#' title:
#' description: 
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change:
stats <- c(1:7)        # Which module statistics (7) to use for perm testing.
strength <- "weak"     # Preservation criterion: strong or weak.
res <- 47              # Resolutions to be analyzed.
partition <- "6142226" # Which partition file to use as input?
contrast <- "6171865"  # Which network comparison file to use as input?

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

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

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

# Calculate power for approximate scale free fit.
sft <- silently({
	sapply(list(wtDat,koDat),function(x) 
	       pickSoftThreshold(x,
				 corFnc="bicor",
				 networkType="signed",
				 RsquaredCut=0.8)$powerEstimate)
})
names(sft) <- c("wt","ko")

# Load network partitions. Self-preservation enforced.
myfile <- list.files(rdatdir, pattern = partition, full.names = TRUE)
partitions <- readRDS(myfile)

# Load network comparison results.
myfile <- list.files(rdatdir,pattern=contrast,full.names=TRUE)
comparisons <- readRDS(myfile)

#------------------------------------------------------------------------------
# Perform permutation test.
#------------------------------------------------------------------------------

# WT and KO partitions.
wtPartition <- partitions[[r]][["wt"]]
koPartition <- partitions[[r]][["ko"]]

# Total number of modules.
nModules <- c(
"wt" = sum(names(table(wtPartition)) != 0),
"ko" = sum(names(table(koPartition)) != 0)
)

# Split into modules.
wtModules <- split(wtPartition, wtPartition)
koModules <- split(koPartition, koPartition)

# Input for NetRep:
data_list <- list(wt = wtDat, ko = koDat)
correlation_list <- list(wt = wtAdjm, ko = koAdjm)
network_list <- list(wt = abs(wtAdjm^sft["wt"]), ko = abs(koAdjm^sft["ko"]))
module_list <- list(wt = koPartition, ko = wtPartition)

# Hypothesis for self-preservation.
h0 <- list(
wt = c(discovery = "wt", test = "wt"),
ko = c(discovery = "ko", test = "ko")
)

# Perform permutation testing.
preservation <- lapply(h0, function(x) {
NetRep::modulePreservation(
network = network_list,
data = data_list,
correlation = correlation_list,
moduleAssignments = module_list,
modules = NULL,
backgroundLabel = "0",
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

# Identify preserved and divergent modules.
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
    return(v)
  } # ENDS function

  # Collect strong or weak changes...
  module_changes <- check_modules(preservation)
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
  pdKO <- sum(koProts == "divergent") / length(koProts)
  ppKO <- sum(koProts == "preserved") / length(koProts)
  total_divergent <- (sum(wtProts == "divergent") + sum(koProts == "divergent"))/(2*length(wtProts))
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
  # Total divergent.
  message(paste0(
    "... Total percentage of proteins assigned to divergent modules: ",
    round(100*total_divergent), " (%).","\n"
	  ))
  # Return.
  output[[r]] <- list("wtPartition" = wtPartition, "wtProts" = wtProts, 
		      "koPartition" = koPartition, "koProts" = koProts)
} # ENDS LOOP.
