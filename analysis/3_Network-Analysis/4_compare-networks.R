#!/usr/bin/env Rscript

#' ---
#' title:
#' description: Compare modules in WT and KO networks with permutation test.
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

## Permutation Statistics:
# 1. avg.weight
# 2. coherence
# 3. cor.cor - Don't use!
# 4. cor.degree - Don't use!
# 5. cor.contrib - Don't use!
# 6. avg.cor
# 7. avg.contrib

# User parameters to change:
stats <- c(1:7) # Which permutation statistics to use for perm testing.
strength <- "strong" # Preservation criterion: strong = all, weak = any sig stats.
res <- c(1:100) # c(29, 35, 36, 40, 41, 42, 44, 45, 48, 49, 55, 58, 66, 79)
cutoff <- 1 # Size cutoff for a module 1 = single protein.
net1 <- "Cortex"
net2 <- "Striatum"
partition1 <- "6142226" # Which partition file to use as input? Used self-pres enforced partition.
partition2 <- "6142226" # Which partition file to use as input? Used self-pres enforced partition.
save_results <- FALSE # Should permutation results be saved?

# Is this a slurm job?
slurm <- any(grepl("SLURM", names(Sys.getenv())))
if (slurm) {
  # SLURM job notes - sent to job_*.info
  nThreads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
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

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load expression data.
mydata <- paste0("3_",c(net1,net2),"_cleanDat.RData")
myfiles <- sapply(mydata,function(x) list.files(rdatdir,x,full.names=TRUE))
data <- lapply(myfiles,readRDS)
data <- lapply(data,t)
names(data) <- c(net1,net2)

# Load adjmatrices.
mynetworks <- paste0("3_",c(net1,net2),"_Adjm.RData")
myfiles <- sapply(mynetworks,function(x) list.files(rdatdir,x,full.names=TRUE))
adjm <- lapply(myfiles,function(x) as.matrix(readRDS(x)))
adjm <- lapply(adjm,function(x) { rownames(x) <- colnames(x); return(x) })
names(adjm) <- c(net1,net2)

# Create unsigned networks for NetRep.
networks <- lapply(adjm,abs)

# partitions files: 
# 10342568 = Striatum - self-preservation enforced.

# Load network partitions.
# If working with csv files:
mypartitions <- paste0(c(net1,net2),"_partitions.csv")
myfiles <- sapply(mypartitions,function(x) list.files(rdatdir,x,full.names=TRUE))
all_partitions <- lapply(myfiles,function(x) fread(x,drop=1,skip=1))
all_partitions <- lapply(c(1,2),function(x) { 
			     colnames(all_partitions[[x]]) <- colnames(adjm[[x]])
			     return(all_partitions[[x]])
		      })
names(all_partitions) <- c(net1,net2)

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
  partitions <- lapply(all_partitions, function(x) as.integer(x[r, ] + 1))
  prots <- colnames(all_partitions[[1]])
  partitions <- lapply(partitions,function(x) { names(x) <- prots; return(x) })
  # Remove small modules.
  partitions <- lapply(partitions,filter_modules)
  # Total number of modules; ignore 0.
  nModules <- sapply(partitions,function(x) sum(names(table(x)) != 0))
  # Split partitions into modules.
  modules <- lapply(partitions,function(x) split(x,x))
  # Input for NetRep:
  data_list <- data
  correlation_list <- adjm
  network_list <- networks 
  module_list <- partitions # Zero index modules will be ignored.

  # Perform permutation testing.
  # Suppress warnings which arise from small modules;
  # NA p.vals are replaced with 1.
  H0 <- list(c(discovery = net1, test = net2),
	     c(discovery = net2, test = net1))
  names(H0) <- c(net1,net2) # Test cox in str and str in cox...
  suppressWarnings({
  preservation <- sapply(H0,function(x) {
					 NetRep::modulePreservation(
		                         network = network_list,
		                         data = data_list,
		                         correlation = correlation_list,
					 moduleAssignments = module_list,
					 modules = NULL,
				 	 backgroundLabel = "0",
					 discovery = x["discovery"],
					 test = x["test"],
					 selfPreservation = FALSE,
					 nThreads = nThreads,
					 # nPerm = 100000,  # determined by the function.
					 null = "overlap",
					 alternative = "two.sided", # c(greater,less,two.sided)
					 simplify = TRUE,
					 verbose = TRUE)
})
  }) # Ends sapply.

  names(preservation)
  class(preservation)
  length(preservation)

  # Identify preserved and divergent modules.
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
  pdKO <- sum(koProts == "divergent") / length(koProts)
  ppKO <- sum(koProts == "preserved") / length(koProts)
  total_divergent <- (sum(wtProts == "divergent") + sum(koProts == "divergent")) / (2 * length(wtProts))
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
    round(100 * total_divergent), " (%).", "\n"
  ))
  # Return.
  output[[r]] <- list(
    "wtPartition" = wtPartition, "wtProts" = wtProts,
    "koPartition" = koPartition, "koProts" = koProts
  )
} # ENDS LOOP.

# Save output to file.
output_name <- paste0(job, "_Network_Comparisons.RData")
myfile <- file.path(rdatdir, output_name)
saveRDS(output, myfile)
