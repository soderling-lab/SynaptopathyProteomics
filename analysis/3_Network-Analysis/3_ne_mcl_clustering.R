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
net = "Striatum" # Which network are we analyzing? 
max_size = 500 # maximum allowable size of modules before apply MCL.
inflation = seq(1.2,5,0.2) # Inflation space to explore.
resolutions = seq(1,100) # Resolutions to analyze.

# Global options and imports.
suppressPackageStartupMessages({
	library(neten)
	library(igraph)
	library(parallel)
	library(doParallel)
})

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

# Load Leidenalg graph partitions from 2_la-clustering.
myfiles <- c("Cortex" = file.path(rdatdir,"147731383_Cortex_CPMVertexPartition_partitions.csv"),
	    "Striatum" = file.path(rdatdir,"148436673_Striatum_CPMVertexPartition_partitions.csv"))
partitions <- data.table::fread(myfiles[net], header=TRUE,drop = 1)

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

# Register parallel clusters.
workers <- makeCluster(c(rep("localhost",nThreads)),type="SOCK")
registerDoParallel(workers)

# Parallel execution:
results <- foreach(resolution=resolutions) %dopar% {
	cluster_mcl(resolution)
}

# Shut down nodes.
suppressWarnings(stopCluster(workers))

# Loop to combine mcl partitions at every resolution.
mcl_partitions <- lapply(results,combine_partitions)

# Loop to merge mcl partition into original partition at every resolution.
final_partitions <- list()
for (i in seq_along(mcl_partitions)) {
	p1 <- all_partitions[[i]]
	p2 <- mcl_partitions[[i]]
	p3 <- merge_partitions(p1,p2)
	final_partitions[[i]] <- p3
}

# Write to file.
# Keep row indices so that output matches that from La clustering script.
myfile <- file.path(rdatdir,paste0("3_",net,"_MCL_partitions.csv"))
df <- as.data.frame(do.call(rbind,final_partitions))
data.table::fwrite(df,myfile,row.names=TRUE)
