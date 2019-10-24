#!/usr/bin/env Rscript
# Examine module self-preservation.

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(readxl)
	library(igraph)
	library(reshape2)
	library(ggplot2)
	library(anRichment)
	library(TBmiscr)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
datadir <- file.path(root,"data")
tabsdir <- file.path(root,"tables")

# Load expression data
wtDat <- t(readRDS(file.path(datadir,"wtDat.Rds")))
koDat <- t(readRDS(file.path(datadir,"koDat.Rds")))

# Fix rownames.
colnames(wtDat) <- colnames(koDat) <- rownames(readRDS(file.path(datadir,"wtDat.Rds")))

# Compute adjmatrix:
wtAdjm <- silently(WGCNA::bicor, wtDat)
koAdjm <- silently(WGCNA::bicor, koDat)

# Load partitions.
wtParts <- data.table::fread(file.path(datadir,"WT_partitions.csv"),drop=1)
koParts <- data.table::fread(file.path(datadir,"KO_partitions.csv"),drop=1)

#-------------------------------------------------------------------------------
# Input for NetRep:
data_list        <- list(wt = wtDat, ko = koDat)   
correlation_list <- list(wt = wtAdjm, ko = koAdjm) 
network_list     <- list(wt = wtAdjm, ko = koAdjm) 

# Loop through partitions, evaluating self-preservation.
results <- list()

for (i in 1:100) {
	# status
	message(paste("working on partition",i,"..."))
	# Get partition.
	wtPartition <- as.integer(wtParts[i,])+1
	koPartition <- as.integer(koParts[i,])+1
	names(wtPartition) <- names(koPartition) <- colnames(wtAdjm)
	module_list <- list(wt = wtPartition, ko = koPartition)
	# Perform permutation test for module self-preservation.
	self = as.list(c("wt","ko"))
	selfPreservation <- lapply(self,function(x) {
					   NetRep::modulePreservation(
					   network = network_list,
					   data = data_list,
					   correlation = correlation_list,
					   moduleAssignments = module_list,
					   modules = NULL,
					   backgroundLabel = 0,
					   discovery = x,
					   test = x,
					   selfPreservation = TRUE,
					   nThreads = 8,
					   #nPerm = 100000, 
					   null = "overlap",
					   alternative = "greater",
					   simplify = TRUE,
					   verbose = TRUE)})
	# Function to get max pvalue.
	maxp <- function(preservation) {
	p <- apply(preservation$p.values,1,function(x) max(x,na.rm=TRUE))
	q <- p.adjust(p,"bonferroni")
	return(q)
	}
	q <- lapply(selfPreservation, maxp)
	# Modules with NS preservation stats. 
	out <- lapply(q,function(x)names(x)[x>0.05])
	# For NS modules, set module membership to 0.
	wtPartition[wtPartition %in% out[[1]]] <- 0
	koPartition[koPartition %in% out[[2]]] <- 0
	# Return results.
	results[[i]] <- list(wt = wtPartition, ko = koPartition)
}

# Save to Rdata.
saveRDS(results,file.path(datadir,"self_preservation_results.RDS")
