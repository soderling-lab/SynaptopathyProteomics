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

## Loop through resolutions, for each module get best mcl partition.
results <- list()

# Function to do mcl clustering for a given resolution.
cluster_mcl <- function(resolution) {
suppressPackageStartupMessages({
	library(neten)
	library(igraph)
	library(parallel)
	library(doParallel)
})
	#message(paste("\nWorking on resolution",resolution,"..."))
	part <- all_partitions[[resolution]]
	modules <- split(part,part)
	module_sizes <- table(part)
	too_big <- names(module_sizes)[module_sizes > max_size]
	#message(paste0("Breaking down (",length(too_big),") large communities."))
	mcl_partitions <- list()
	## Loop through big modules, bet best MCL partition.
	for (m in too_big){
		#message(paste("... Working on module",m,"..."))
		# Get subgraph.
		prots <- names(modules[[m]])
		g1 <- induced_subgraph(g0,vids=V(g0)[match(prots,names(V(g0)))])
		# Network Enhancement.
		g_ne <- neten(g1)
		## Loop to sample inflation space.
		output <- list()
		for (i in seq_along(inflation)){
			# MCL clustering.
			mcl_part <- clusterMCL(g_ne,weight="weight",inflation[i])
			# Modularity of re-weighted (NE) graph.
			q <- modularity(g_ne,membership=mcl_part[names(V(g_ne))],
					weights=abs(edge_attr(g_ne,'weight')))
			# Summary stats.
			k <- sum(names(table(mcl_part)) != "0")
			output[[i]] <- list("Partition"=mcl_part,
					    "Inflation"=inflation[i],
					    "Clusters"=k,
					    "Modularity"=q)
		} # Ends loop through inflation space.
		# Get ~best partition.
		Q <- sapply(output,function(x) x$Modularity)
		K <- sapply(output,function(x) x$Clusters)
		parts <- lapply(output,function(x) x$Partition)
		idx <- which(Q==max(Q))
		best_q <- Q[idx]
		best_i <- inflation[idx]
		best_k <- K[idx]
		# Status.
		#message(paste("... ... Inflation:", best_i))
		#message(paste("... ... Clusters :", best_k))
		#message(paste("... ... Modularity:",round(best_q,3)))
		# return best mcl partition.
		mcl_partitions[[m]] <- parts[[idx]]
	} # Ends loop through modules.
	# Collect best mcl partitions.
	return(mcl_partitions)
} # Ends function.

# Serial execution:
#t0 <- Sys.time()
#results <- foreach(resolution=seq(1,2)) %do% {
#	cluster_mcl(resolution)
#}
#t1 <- Sys.time()
#dt <- t1-t0
#message(paste("Time to execute in series:",dt))

# Register parallel clusters.
workers <- makeCluster(c(rep("localhost",nThreads)),type="SOCK")
registerDoParallel(workers)

# Parallel execution:
t0 <- Sys.time()
results <- foreach(resolution=seq(1,2)) %dopar% {
	cluster_mcl(resolution)
}
t1 <- Sys.time()
dt <- t1-t0
message(paste("Time to execute in parallel:",dt))

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
