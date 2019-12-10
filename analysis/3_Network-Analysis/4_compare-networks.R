#!/usr/bin/env Rscript

#' ---
#' title:
#' description: Compare modules in networks with permutation test.
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change:
stats <- c(1:7)          # Which permutation statistics to use for perm testing.
strength <- "strong"     # Preservation criterion: strong = all, weak = any sig stats.
res <- c(1:100)          # Which resolutions to analyze?
net1 <- "Cortex"         # Network 1.
net2 <- "Striatum"       # Network 2.
partition1 <- "10360847" # Partition file for first network. 
partition2 <- "10342568" # Partition file for second network.
#save_results <- FALSE    # Should permutation results be saved?

## Permutation Statistics:
# 1. avg.weight
# 2. coherence
# 3. cor.cor - Don't use!
# 4. cor.degree - Don't use!
# 5. cor.contrib - Don't use!
# 6. avg.cor
# 7. avg.contrib

## Partitions files:
# 10342568 = Striatum - self-preservation enforced.
# 10360847 = Cortex - self-preservation enforced.
# 6142226  = Combined - self-preservation enforced.

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
mydata <- paste0("3_", c(net1, net2), "_cleanDat.RData")
myfiles <- sapply(mydata, function(x) list.files(rdatdir, x, full.names = TRUE))
data <- lapply(myfiles, readRDS)
data <- lapply(data, t)
names(data) <- c(net1, net2)

# Load adjmatrices.
mynetworks <- paste0("3_", c(net1, net2), "_Adjm.RData")
myfiles <- sapply(mynetworks, function(x) list.files(rdatdir, x, full.names = TRUE))
adjm <- lapply(myfiles, function(x) as.matrix(readRDS(x)))
adjm <- lapply(adjm, function(x) {
  rownames(x) <- colnames(x)
  return(x)
})
names(adjm) <- c(net1, net2)

# Create unsigned networks for NetRep.
# Not weighted.
networks <- lapply(adjm, abs)

# Load network partitions.
mypartitions <- c(partition1,partition2)
myfiles <- sapply(mypartitions, function(x) list.files(rdatdir, x, full.names = TRUE))
all_partitions <- lapply(myfiles, readRDS)
names(all_partitions) <- c(net1, net2)

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
message(paste("Comparing", net1, "and", net2, "networks at", nres, "resolution(s)."))
message(paste0(
  "Module statistic(s) used to evaluate module preservation: ",
  module_stats
), ".")
message(paste0(
  "Criterion for module preservation/divergence: ",
  strength, ".", "\n"
))

# Loop to analyze all resolutions:
output <- list()
for (r in res) {
  # Status report.
  message(paste("Working on resolution:", r, "..."))
  # Extract from list.
  partitions <- list(all_partitions[[net1]][[r]],
		     all_partitions[[net2]][[r]])
  names(partitions) <- c(net1,net2)
  # Total number of modules; ignore 0.
  nModules <- sapply(partitions, function(x) sum(names(table(x)) != 0))
  # Input for NetRep:
  data_list <- data
  correlation_list <- adjm
  network_list <- networks
  module_list <- partitions # Zero index modules will be ignored.
  # Perform permutation testing.
  # Suppress warnings which arise from small modules;
  # NA p.vals are replaced with 1.
  H0 <- list(
    c(discovery = net1, test = net2),
    c(discovery = net2, test = net1)
  )
  names(H0) <- c(net1, net2) # Test cox in str and str in cox...
  suppressWarnings({
    preservation <- lapply(H0, function(x) {
      result <- NetRep::modulePreservation(
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
        verbose = FALSE
      )
      return(result)
    })
  }) # Ends lapply.
  # Identify preserved and divergent modules.
  module_changes <- lapply(preservation, check_modules)
  # Summarize number of preserved, divergent, ns modules.
  nPres <- sapply(module_changes, function(x) sum(x == "preserved"))
  nNS <- sapply(module_changes, function(x) sum(x == "ns"))
  nDiv <- sapply(module_changes, function(x) sum(x == "divergent"))
  ## Status messages:
  # Preservation of Network 1 in Network 2:
  message(paste0("...", "Summary of ", net1, " module preservation in ",net2,":"))
  message(paste("... ...", "N modules significantly preserved:", nPres[net1]))
  message(paste("... ...", "N modules significantly divergent:", nDiv[net1]))
  message(paste("... ...", "N modules NS (no sig. changes)   :", nNS[net1]))
  # Preservation of Network 2 in Network 1:
  message(paste0("...", "Summary of ", net2, " module preservation in ",net1,":"))
  message(paste("... ...", "N modules significantly preserved:", nPres[net2]))
  message(paste("... ...", "N modules significantly divergent:", nDiv[net2]))
  message(paste("... ...", "N modules NS (no sig. changes)   :", nNS[net2]))
  # Return output
  output[[r]] <- module_changes
} # ENDS LOOP.

# Save output to file.
output_name <- paste0(job, "_Network_Comparisons.RData")
myfile <- file.path(rdatdir, output_name)
saveRDS(output, myfile)
