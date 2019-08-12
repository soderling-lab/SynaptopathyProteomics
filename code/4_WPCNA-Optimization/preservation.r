#!/usr/bin/env Rscript

## Examine module self-preservation.

#-------------------------------------------------------------------------------
## # Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
	library(igraph)
	library(NetRep)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
tables <- paste(root,"tables",sep="/") 

# Load expression data and compute adjmatrix:
cleanDat <- readRDS("wtDat.Rds")
r <- WGCNA::bicor(cleanDat)
adjm <- r

# Read network partition info.
clufile <- paste(tables,"partitions.csv", sep = "/")
cluDat <- read.csv(clufile) + 1 # add 1 so that all module membership > 0 
partitions <- split(cluDat, seq(nrow(cluDat)))

# Create graph.
netfile <- paste(here,"wtAdjm.net",sep = "/")
g <- read_graph(netfile,format = "pajek") # Did I need this?

#-------------------------------------------------------------------------------
# ## Enforce module preservation.
#-------------------------------------------------------------------------------
# Remove modules that are not preserved (i.e. have insignificant module
# preservation statistics). The number of random permutations used to generate the
# null distributions is increased to 100,000 in order to stabilize the result with
# large number of modules. This computation is expensive and will take several
# minutes.


# Input for NetRep:
data_list        <- list(data = cleanDat) # The protein expression data.
correlation_list <- list(data = r) # The bicor correlation matrix.
network_list     <- list(data = adjm)  # The weighted, signed co-expresion network.

# Empty list for output of loop.
preserved_partitions <- list()

for (i in seq_along(partitions)){
	message(paste("Working on partition", i ,"..."))
	# Get partition
	module_labels <- as.vector(as.matrix(partitions[[i]])) 
	names(module_labels) <- colnames(cleanDat)
	nModules <- length(unique(module_labels))
	module_list <- list(data = module_labels)
	# Perform permutation test for self-preservation.
	preservation <- NetRep::modulePreservation(
						   network = network_list,
		  				   data = data_list,
				      		   correlation = correlation_list,
				  		   moduleAssignments = module_list,
						   modules = NULL,
						   backgroundLabel = "grey",
						   discovery = "data",
						   test = "data",
						   selfPreservation = TRUE,
						   nThreads = 8,
						   #nPerm = 100000, 
						   null = "overlap",
						   alternative = "greater",
						   simplify = TRUE,
						   verbose = FALSE
						   )
	# Get the maximum permutation test p-value.
	maxp <- apply(preservation$p.values, 1, function(x) max(x, na.rm = TRUE))
	# Modules removed if adjusted pvalue is greater than alpha = 0.05.
	alpha <- 0.05
	modules_out <- names(maxp)[maxp > alpha / nModules]
	nModules_out <- length(modules_out)
	if (length(nModules_out) > 0) {
		idx <- module_labels %in% modules_out
		module_labels[idx] <- 0
	}
	# Return module membership.
	preserved_partitions[[i]] <- module_labels
}

saveRDS(preserved_partitions, "preserved_partitions.Rds")

# ENDOFILE
#------------------------------------------------------------------------------
