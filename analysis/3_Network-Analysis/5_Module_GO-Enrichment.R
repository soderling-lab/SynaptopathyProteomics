#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Parse command line arguments:
if (interactive()) {
	## If interactive:
	# User defined parameters (you only need to change these two):
	analysis_type = "Cortex" # Tissue type for analysis.
} else if (!interactive()) {
	## If not interactive, check that only 1 arg is passed.
	args <- commandArgs(trailingOnly=TRUE)
	if (length(args) == 1) { 
		analysis_type = commandArgs(trailingOnly=TRUE)[1]
		start <- Sys.time()
		message(paste("Starting analysis at:", start))
		message(paste0("Analyzing ", analysis_type,"..."))
	} else { 
		stop("Specify either 'Cortex' or 'Striatum'.",call.=FALSE) 
	}
}

## Optional parameters:

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root)

# Global imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Functions.
suppressWarnings({ devtools::load_all() })

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#--------------------------------------------------------------------
## Load the data.
#--------------------------------------------------------------------

# Load the cortex and striatum data from root/data.
data(cortex_data)
data(striatum_data)

# Load the graph partitions:
data(cortex_partition)
data(striatum_partition)

# Grab the data for the tissue type we are analyzing.
data_list <- list("Cortex"=cortex_data,
		  "Striatum"=striatum_data)[[analysis_type]]
partition <- list("Cortex"=cortex_partition,
		  "Striatum"=striatum_partition)[[analysis_type]]

# Reset partition index for self-preserved modules.
partition <- reset_index(partition)

# Data matrix:
dm <- data_list$Data

# Load the sample meta data.
data(samples)

# Load adjacency matrix.
adjm <- data_list$Adjm

# Load network.
netw <- data_list$Netw

# Load gene identifier map.
data(gene_map)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

# Create list of modules.
modules <- split(partition,partition)
names(modules) <- paste0("M",names(modules))

# Remove M0.
modules <- modules[-which(names(modules) == "M0")]

#--------------------------------------------------------------------
## Save results.
#--------------------------------------------------------------------

# Done!
end <- Sys.time()
message(paste("\nCompleted analysis at:", end))
message(paste("Elapsed time:",
	      round(difftime(end,start,units="mins"),2),"minutes."))
