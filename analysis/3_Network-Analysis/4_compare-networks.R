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

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change:
stats <- c(1:7)      # Which of the seven module statistics to use.
strength <- "strong" # Preservation criterion strong = all, or weak = any sig stats.
res <- 44            # Resolutions to analyze.

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
wtDat <- t(readRDS(list.files(rdatdir, pattern = "WT_cleanDat", 
			      full.names = TRUE)))
koDat <- t(readRDS(list.files(rdatdir, pattern = "KO_cleanDat", 
			      full.names = TRUE)))

# Load adjmatrices.
wtAdjm <- t(readRDS(list.files(rdatdir, pattern = "WT_Adjm.RData", 
			       full.names = TRUE)))
koAdjm <- t(readRDS(list.files(rdatdir, pattern = "KO_Adjm.RData", 
			       full.names = TRUE)))

# Load network partitions. Self-preservation enforced.
myfile <- list.files(rdatdir, pattern = "preservation", full.names = TRUE)
partitions <- readRDS(myfile)

#------------------------------------------------------------------------------
# Loop through all resolutions and perform permutation test.
#------------------------------------------------------------------------------

# All module statistics.
module_stats <- c("avg.weight","coherence","cor.cor",
		  "cor.degree","cor.contrib","avg.cor","avg.contrib")[stats]

# Status report:
nres <- length(res)
message(paste("Analyzing all resolutions in:", nres))
message(paste("Module statistic used to evaluate module preservation:",
	      module_stats))
message(paste("Criterion for module preservation/divergence:", 
	      c("strong"="all","weak"="any")[strength],"\n"))

# LOOP TO ANALYZE ALL RESOLUTIONS:
output <- list()

# for (r in seq_along(1:nres)) {
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
	  obs <- x$observed[,stats]
	  nulls <- apply(x$nulls, 2, function(x) apply(x, 1, mean))[,stats]
	  q <- apply(x$p.values,2,function(x) p.adjust(x,"bonferroni"))[,stats]
       	  q[is.na(q)] <- 1
	  # If testing more than one statistic.
	  if (length(stats)>1) {
		  fx <- c("strong"="all","weak"="any")[strength]
		  sig <- apply(q<0.05,1,eval(fx))
		  greater <- apply(obs>nulls,1,eval(fx))
		  less <- apply(obs<nulls,1,eval(fx))
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
  names(module_changes) <- c("ko","wt")
  # Calculate percent NS, divergent, preserved.
  wtProts <- wtPartition
  wtProts[wtProts == "0"] <- "ns"
  wtProts[wtProts %in% names(wtModules)[module_changes$wt == "preserved"]] <- "preserved"
  wtProts[wtProts %in% names(wtModules)[module_changes$wt == "divergent"]] <- "divergent"
  wtProts[wtProts %in% names(wtModules)[module_changes$wt == "ns"]] <- "ns"
  koProts <- koPartition
  koProts[koProts == "0"] <- "ns"
  koProts[koProts %in% names(koModules)[module_changes$ko == "preserved"]] <- "preserved"
  koProts[koProts %in% names(koModules)[module_changes$ko == "divergent"]] <- "divergent"
  koProts[koProts %in% names(koModules)[module_changes$ko == "ns"]] <- "ns"
  n_divergent <- sum(wtProts == "divergent") + sum(koProts == "divergent")
  n_preserved <- sum(wtProts == "preserved") + sum(koProts == "preserved")
  n_ns <- sum(wtProts == "ns") + sum(koProts == "ns")
  percent_divergent <- n_divergent / (n_divergent + n_preserved + n_ns)
  percent_preserved <- n_preserved / (n_divergent + n_preserved + n_ns)
  percent_ns <- n_ns / (n_divergent + n_preserved + n_ns)
  # Status report.
  message(paste("... Total number of WT modules:", nModules["wt"]))
  message(paste(
    "... ... Number of WT modules preserved in KO graph:",
    sum(module_changes$ko == "preserved")
  ))
  message(paste(
    "... ... Number of WT modules divergent in KO graph:",
    sum(module_changes$ko == "divergent")
  ))
  message(paste("... Total number of KO modules:", nModules["ko"]))
  message(paste(
    "... ... Number of KO modules preserved in WT graph:",
    sum(module_changes$wt == "preserved")
  ))
  message(paste(
    "... ... Number of KO modules divergent in WT graph:",
    sum(module_changes$wt == "divergent")
  ))
  message(paste("... Percent preserved:", percent_preserved))
  message(paste("... Percent divergent:", percent_divergent))
  message(paste("... ... .. Percent ns:", percent_ns, "\n"))
  # Return resolution, total number of modules, and module changes.
  output[[r]] <- list(
    "resolution" = r, "nModules" = nModules, 
    "wtChanges" = module_changes$wt,
    "koChanges" = module_changes$ko,
    "wtPartition" = wtPartition,
    "wtProteins" = wtProts,
    "koPartition" = koPartition,
    "koProteins" = koProts,
    "percent_preserved" = percent_preserved,
    "percent_divergent" = percent_divergent,
    "percent_ns" = percent_ns
  )
} # ENDS LOOP.

# Save output to file.
output_name <- paste0(job, "_Network_Comparisons.RData")
myfile <- file.path(rdatdir, output_name)
saveRDS(output, myfile)
