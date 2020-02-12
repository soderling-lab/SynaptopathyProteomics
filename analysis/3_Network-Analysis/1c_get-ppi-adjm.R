#!/usr/bin/env Rscript

#' ---
#' title: 1c_get-ppi-adjm.R
#' description: Generate PPI adjacency matrix.
#' authors: Tyler W Bradshaw
#' ---

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
	library(getPPIs)
	library(ggplot2)
	library(data.table)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatdir <- file.path(root, "rdata")
funcdir <- file.path(root, "R")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein identifier map.
prot_map <- readRDS(file.path(rdatdir,"2_Protein_ID_Map.RData"))

#---------------------------------------------------------------------
## Generate PPI graph.
#---------------------------------------------------------------------

# Load mouse interactome.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
orgs <- c(10090, 9606, 10116)
idx <- musInteractome$Interactor_A_Taxonomy %in% orgs
ppis <- subset(musInteractome, idx)

# Get entrez IDs for all proteins in co-expression networks.
entrez <- prot_map$entrez

# Build a ppi graph.
message("Building PPI graph...")
g <- buildNetwork(ppis, entrez, taxid = 10090)

# Louvain clustering.
#partition = cluster_louvain(g, weights = NULL)
#m <- partition$membership
#names(m) <- partition$names
#modules = split(m,m)
#module_sizes <- sapply(modules,length)

# Get ppi adjacency matrix.
PPIadjm <- as.matrix(as_adjacency_matrix(g))

# Map entrez back to protein ids.
ids <- prot_map$ids[match(colnames(PPIadjm),prot_map$entrez)]
colnames(PPIadjm) <- rownames(PPIadjm) <- ids

# One protein is mapped to the same entrez id.
missing <- prot_map$id[prot_map$ids %notin% colnames(PPIadjm)]
names(missing) <- prot_map$entrez[match(missing,prot_map$id)]
duplicated <- prot_map$id[which(prot_map$entrez ==  names(missing))]

# Add missing column and row. 
# Just copy column and row for other protein.
missing_col <- PPIadjm[,duplicated[duplicated %notin% missing]]
PPIadjm <- cbind(PPIadjm,missing_col)
colnames(PPIadjm)[ncol(PPIadjm)] <- missing
missing_row <- PPIadjm[duplicated[duplicated %notin% missing],]
PPIadjm <- rbind(PPIadjm,missing_row)
rownames(PPIadjm)[nrow(PPIadjm)] <- missing

# Evaluate scale free fit.
dc <- apply(PPIadjm,2,sum) # Node degree.
fit <- WGCNA::scaleFreeFitIndex(dc,nBreaks=10,removeFirst=FALSE)
r <- fit$Rsquared.SFT
message(paste("Scale free fit of PPI graph:",round(r,3)))

# Write to file.
myfile <- file.path(rdatdir,"3_PPI_Adjm.csv")
fwrite(as.data.table(PPIadjm),myfile,row.names=TRUE)

message("Done!\n")
