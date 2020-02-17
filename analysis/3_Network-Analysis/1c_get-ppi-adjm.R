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
subdir <- dirname(here)
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables",subdir)

# Functions.
suppressWarnings({ devtools::load_all() })

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

# Save PPIs data frame--this contains evidence information.
myfile <- file.path(rdatdir,"3_All_PPIs.RData")
saveRDS(ppis,myfile)

# Get entrez IDs for all proteins in co-expression networks.
entrez <- prot_map$entrez

# Build a ppi graph.
message("Building PPI graph...")
g <- buildNetwork(ppis, entrez, taxid = 10090)

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

# Save as Rdata.
myfile <- file.path(rdatdir,"3_PPI_Adjm.RData")
saveRDS(PPIadjm,myfile)

message("Done!\n")
