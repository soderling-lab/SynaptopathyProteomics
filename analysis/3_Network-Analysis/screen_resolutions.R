#!/usr/bin/env Rscript
# Loop through all resolutions of WT and KO graphs, comapring modules with 
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
myfile <- list.files(rdatdir, pattern = "preservation", full.names = TRUE)
partitions <- readRDS(myfile)

# LOOP TO ANALYZE ALL RESOLUTIONS:
# Using partitions with self-resolution already enforced.
output <- list()
for (r in 1:100) {
	r_wt <- r_ko <- r
	# Extract from list.
	wtPartition <- partitions[[r_wt]][["wt"]]
	koPartition <- partitions[[r_ko]][["ko"]]
	# Split into modules.
	wtModules <- split(wtPartition, wtPartition)
	koModules <- split(koPartition, koPartition)
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
	module_list <- list(wt = koPartition, ko = wtPartition) # Switch!
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
	    backgroundLabel = 0,
	    discovery = x["discovery"],
	    test = x["test"],  
	    selfPreservation = TRUE, 
	    nThreads = 24,
	    # nPerm = 100000,  # determined by the function.
	    null = "overlap",
	    alternative = "two.sided", # c(greater,less,two.sided)
	    simplify = TRUE,
	    verbose = FALSE
	    )
	})
	output[[r]] <- preservation
	## Identify divergent modules.
	# Divergent modules are those whose observed correlation structure is significantly
	# less? than the null model.
	# Require all stats to be less than 0.05.
	 q <- lapply(preservation, function(x) p.adjust(apply(x$p.values, 1, max)))
	# Which modules have significant statistics?
	sigModules <- lapply(q, function(x) names(x)[x < 0.05])
	# Status.
	sigWT <- sigModules$wt
	sigKO <- sigModules$ko
	message(paste("Number of WT modules that exhibit divergence:", length(sigWT)))
	message(paste("Number of KO modules that exhibit divergence:", length(sigKO)))
} # END LOOP.

# Save output.
myfile <- file.path(rdatdir,"Module_Divergence_Preservation.RData")
saveRDS(output,myfile)
