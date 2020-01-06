#!/usr/bin/env Rscript

# Compare Co-expression networks to GO functional similarity network.

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters. Which co-expression network to analyze?
net <- "Cortex" 


# Imports. 
suppressPackageStartupMessages({
	library(data.table)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatdir <- file.path(root,"rdata")
funcdir <- file.path(root,"R")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein identifier map.
prot_map <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load co-expression network partitions--self-preservation enforced.
mypart <- c(Cortex = "10773682",Striatum = "10781799")[net] # relaxed criterion
myfile <- list.files(rdatdir,pattern=mypart,full.names=TRUE)
partitions <- readRDS(myfile)

# Load GO semantic similarity adjm.
myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_RMS_Adjm.csv")
GOadjm <- fread(myfile,header=TRUE,drop=1)

# Load GO semantic similarity network partitions.
myfile <- file.path(rdatdir,"3_GO_partitions.csv")
GOparts_df <- fread(myfile,header=TRUE,drop=1)
colnames(GOparts_df) <- colnames(GOadjm)

# Split GOSemSim partition df into list of partitions.
# First, fix names.
entrez_map <- as.list(prot_map$ids)
names(entrez_map) <- prot_map$entrez
ids <- unlist(entrez_map[colnames(GOparts_df)])
GOpartitions <- lapply(c(1:100),function(x) {
	      part <- as.numeric(GOparts_df[x,])
	      names(part) <- ids
	      return(part)
})

# Compare co-expresion and GO similarity partitions.
partition_similarity <- vector(mode="numeric",length=100)
for (i in 1:100) {
p1 = partitions[[i]]
p2 = GOpartitions[[i]]
# Add missing genes. Assign to module 0.
missing_ids <- names(p1)[names(p1) %notin% names(p2)]
missing <- vector(mode="numeric",length(missing_ids))
names(missing) <- missing_ids
p2 <- c(p2,missing)
# Not sure why there is a duplicate...
p2 <- p2[!duplicated(names(p2))]
# Need to fix calculation for M0. All 0, then should be 0 similarity.
module_similarity[i] <- module_assignment_similarity(p1,p2)
}


