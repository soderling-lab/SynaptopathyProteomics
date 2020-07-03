#!/usr/bin/env Rscript

#' ---
#' title:
#' description: Preprocessing of TMT data.
#' authors: Tyler W Bradshaw
#' ---

## Input:

## Options:

#---------------------------------------------------------------------
## Misc function - getrd().
#---------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
	in_root <- function(h=here, dir=dpat) { 
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE))) 
		return(check)
	}
	# Loop to find root.
	while (!in_root(here)) { 
		here <- dirname(here) 
	}
	root <- here
	return(root)
}

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. Load custom functions and prepare
# the project directory for saving output files.

# Load the R env.
rootdir <- getrd()
renv::load(rootdir,quiet=TRUE)

# Load required packages.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Load functions in root/R.
suppressWarnings({ devtools::load_all() })

# Directories:
rdatdir <- file.path(rootdir, "rdata")

#---------------------------------------------------------------------
## Module-level statistical analysis
#---------------------------------------------------------------------

## 
# Load the partition.
data(cortex_partition)
partition = cortex_partition

# Percent clustered
modules = split(partition,partition)
pclustered <- sum(partition!=0)/length(partition)
pclustered

# module glm

