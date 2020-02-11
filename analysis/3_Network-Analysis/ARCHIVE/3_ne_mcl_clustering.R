#!/usr/bin/env Rscript

#' ---
#' title: 
#' description: Breaking down large, cohesive communities with NE + MCL.
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

## User parameters to change:
net = "Cortex" # Which network are we analyzing? 
partition_file = "Cortex_Surprise"
max_size = 500 # maximum allowable size of modules before apply MCL.
inflation = seq(1.2,5,0.2) # Inflation space to explore.

# Global options and imports.
suppressPackageStartupMessages({
	library(neten)
	library(igraph)
	library(parallel)
	library(doParallel)
})

# Number of threads for parallel processing.
nThreads <- detectCores() - 1

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatdir <- file.path(root, "rdata")
funcdir <- file.path(root, "R")

# Functions.
myfun <- list.files(funcdir, pattern="*.R",full.names = TRUE)
invisible(sapply(myfun, source))

# Load co-expression (adjacency) matrices.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_Adjm.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_Adjm.RData"))
adjm <- as.matrix(readRDS(myfiles[net]))
rownames(adjm) <- colnames(adjm)

# Perform network enhancement.
adjm <- neten(adjm)

# Load Leidenalg graph partitions from 2_la-clustering.
myfiles <- c("Cortex" = "147731383_Cortex_CPMVertexPartition_partitions.csv",
	    "Striatum" = "148436673_Striatum_CPMVertexPartition_partitions.csv",
	    "Cortex_Surprise" = "3_Cortex_SurpriseVertexPartition_partitions.csv")
partitions <- data.table::fread(file.path(rdatdir,myfiles[partition_file]), 
				header=TRUE,drop = 1)

# Collect all partitions in a list.
all_partitions <- lapply(seq(nrow(partitions)),function(x) {
  partition <- as.numeric(partitions[x,]) + 1
  names(partition) <- colnames(partitions) 
  return(partition)
})

#--------------------------------------------------------------------
## Loop to break down large communities with NE + MCL.
#--------------------------------------------------------------------

# Create a Co-expression graph.
g0 <- graph_from_adjacency_matrix(adjm,mode="undirected",weighted=TRUE)

# Find modules that are too big.
partition <- all_partitions[[1]]
too_big <- names(table(partition))[which(table(partition) > max_size)]
m <- split(partition,partition)
modules <- m[too_big]

