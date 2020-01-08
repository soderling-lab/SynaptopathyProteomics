#!/usr/bin/env Rscript

#' ---
#' title: 
#' description:
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
	library(getPPIs)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
datadir <- file.path(root, "rdata")
funcdir <- file.path(root, "R")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

#---------------------------------------------------------------------
## Generate PPI graph.
#---------------------------------------------------------------------

# Load mouse interactome.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
idx <- musInteractome$Interactor_A_Taxonomy %in% c(10090, 9606, 10116)
ppis <- subset(musInteractome, idx)

# Get entrez IDs for all proteins in co-expression network.
prots <- colnames(adjm)
entrez <- prot_map$entrez[match(prots, prot_map$ids)]

# Build a ppi graph.
g <- buildNetwork(ppis, entrez, taxid = 10090)

# Get ppi adjacency matrix.
PPIadjm <- as.matrix(as_adjacency_matrix(g))

# Write to file.
myfile <- file.path(rdatdir,"3_PPI_Adjm.csv")
data.table::fwrite(PPIadjm,myfile,row.names=TRUE)
