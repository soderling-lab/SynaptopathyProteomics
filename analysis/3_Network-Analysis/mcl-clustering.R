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
max_size = 500
min_size = 5
inflation = seq(1.2,5,0.2) # Inflation space to explore

# Is this a slurm job?
slurm <- any(grepl("SLURM", names(Sys.getenv())))
if (slurm) {
  # SLURM job notes - sent to job_*.info
  nThreads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
  jobID <- as.integer(Sys.getenv("SLURM_JOBID"))
  info <- as.matrix(Sys.getenv())
  idx <- grepl("SLURM", rownames(info))
  myfile <- file.path("./out", paste0("job_", jobID, ".info"))
  write.table(info[idx, ], myfile, col.names = FALSE, quote = FALSE, sep = "\t")
} else {
  nThreads <- 8
  jobID <- ""
}

# Global options and imports.
suppressPackageStartupMessages({
	library(neten)
	library(igraph)
})

# Directories.
if (rstudioapi::isAvailable()) {
  setwd("D:/projects/SynaptopathyProteomics/analysis/3_Network-Analysis")
}
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


# Enforce consistent order between data, adjm, and partitions.
idx <- idy <- match(colnames(data),colnames(adjm))
adjm <- adjm[idx,idy]
check <- all(colnames(data) == colnames(adjm))
if (!check) { message("Problem: adjm and data columns dont match!") }
idy <- match(colnames(data),colnames(partitions))
partitions <- as.data.frame(partitions)[,idy]
check <- all(colnames(data) == colnames(partitions))
if (!check) { message("Problem: data and partition names dont match!") }

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

# Loop through resolutions, for each module get best mcl partition.
results <- list()
for (resolution in seq(1,100)) {
	part <- all_partitions[[resoltuion]]
	modules <- split(part,part)
	module_sizes <- table(part)
	too_big <- names(module_sizes)[module_sizes > max_size]
	mcl_partitions <- list()
	## Loop through big modules, bet best MCL partition.
	for (m in too_big){
		message(paste("Working on module",m,"..."))
		# Get subgraph.
		prots <- names(modules[[m]])
		g1 <- induced_subgraph(g0,vids=V(g0)[match(prots,names(V(g0)))])
		# Network Enhancement.
		g_ne <- neten(g1)
		## Loop to sample inflation space.
		message("Exploring MCL inflation space...")
		output <- list()
		for (i in seq_along(inflation)){
			# MCL clustering.
			mcl_part <- clusterMCL(g_ne,weight="weight",inflation[i])
			# Modularity of re-weighted (NE) graph.
			q <- modularity(g_ne,membership=mcl_part[names(V(g_ne))],weights=abs(edge_attr(g_ne,'weight')))
			# Summary stats.
			k <- sum(names(table(mcl_part)) != "0")
			# Status.
			message(paste("... Inflation:",inflation[i],";",
				      "Modularity:", round(q,3)))
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
		message(paste("Module:",m))
		message(paste("Best Inflation:",best_i))
		message(paste("Number of Clusters:",best_k))
		message(paste("Modularity:",best_q))
		# return best mcl partition.
		mcl_partitions[[m]] <- parts[[idx]]
	} # ends loop through modules.
	# Collect best mcl partitions.
	results[[resolution]] <- mcl_partitions
} # Ends loop through resolutions.