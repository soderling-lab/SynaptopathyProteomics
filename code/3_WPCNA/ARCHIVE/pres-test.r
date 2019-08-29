#!/usr/bin/env Rscript

## Examine preservation of modules identified in the discovery dataset, either
#  WT or KO protein co-expression graph, in the the opposite test dataset.

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
	library(igraph)
	library(NetRep)
})

# Which type of analysis are we doing?
type <- 1
perform_permutations  <- FALSE # Should permutation test be carried out?
discovery <- c("wt","ko")[type]
test      <- c("ko","wt")[type]
msg <- paste("Analyzing preservation of", discovery, "modules in", test, "network!")
message(msg)

# Directories.
here <- getwd()
root <- dirname(dirname(here))
subdir <- c("WT_in_KO", "KO_in_WT")[type] # For output of the permutation test.
data <- paste(root,"data", "Preservation_Results", subdir, sep="/")
tables <- paste(root,"tables",sep="/") 

# Function to suppress unwanted output from WGCNA functions.
silently <- function(func, ...) {
	sink(tempfile())
	out <- func(...)
	sink(NULL)
	return(out)
}

#-------------------------------------------------------------------------------
## Parse the users input.
#-------------------------------------------------------------------------------

# Load expression data and compute adjmatrix:
wtDat <- readRDS("wtDat.Rds")
koDat <- readRDS("koDat.Rds")

wt_adjm <- silently(WGCNA::bicor, wtDat)
ko_adjm <- silently(WGCNA::bicor, koDat)

# Read network partition info.
files <- c("wt_preserved_partitions.Rds",
	   "ko_preserved_partitions.Rds")[type]
clufile <- paste(here, files, sep = "/")
partitions <- readRDS(clufile)  # len(WT partitions) == 110 

# Checks:
if (!all(colnames(wtDat) == colnames(koDat))) { stop("Input data don't match!") }
if (!all(colnames(wt_adjm) == colnames(ko_adjm))) { stop("Input data don't match!") }
if (!all(names(partitions[[1]]) %in% colnames(wtDat))) { stop("Input data don't match!") }

#-------------------------------------------------------------------------------
# ## Examine module preservation.
#-------------------------------------------------------------------------------

# Input for NetRep:
data_list        <- list(wt = wtDat,   ko = koDat)   # The protein expression data.
correlation_list <- list(wt = wt_adjm, ko = ko_adjm) # The bicor correlation matrix.
network_list     <- list(wt = wt_adjm, ko = ko_adjm) # The weighted, signed co-expresion network.

# Lets look at partition #110: max resolution.
i = 110

# Get partition
module_labels <- partitions[[i]] 
nModules <- length(unique(module_labels))
module_list <- list(wt = module_labels)

# Perform permutation test.
preservation <- NetRep::modulePreservation(
					   network = network_list,
					   data = data_list,
					   correlation = correlation_list,
					   moduleAssignments = module_list,
					   modules = NULL,
					   backgroundLabel = 0,
					   discovery = discovery,
					   test = test,
					   selfPreservation = FALSE,
					   nThreads = 8,
					   #nPerm = 100000, 
					   null = "overlap",
					   alternative = "less", #c(greater,less,two.sided)
					   simplify = TRUE,
					   verbose = TRUE
					   )

#------------------------------------------------------------------------------
## Examine the results.
#------------------------------------------------------------------------------

pvalues <- as.data.frame(preservation$p.values)
pvalues$min.p <- apply(pvalues,1,min)

module <- as.character(136)
idx <- c(1:73)[rownames(pvalues) == module]

preservation$nVarsPresent[idx] # Number of proteins in that module.

moduleList <- split(module_labels,module_labels)

moduleList[module]



# ENDOFILE
#------------------------------------------------------------------------------
