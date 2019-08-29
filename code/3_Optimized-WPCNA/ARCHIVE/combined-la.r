#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
## Build PPI and PCE graph to pass as layers to leidenalg.
#------------------------------------------------------------------------------

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
	library(WGCNA)
	library(igraph)
	library(reshape2)
})

here <- getwd()
root <- dirname(dirname(here))
tables <- paste(root,"tables",sep="/")

# Build complete Protein co-expression (PCE) graph.
file <- paste(here,"wtDat.RDS", sep="/")
exprDat <- readRDS(file)
adjm <- bicor(exprDat)

# Load ppi data.
file <- paste(tables,"3_Compiled_PPIs.csv", sep="/")
ppiDat <- read.csv(file)

# Load protein map.
file <- paste(here,"map.csv",sep="/")
prot_map <- read.csv(file)

# Check.
all(ppiDat$musEntrezA %in% prot_map$entrez)
all(ppiDat$musEntrezB %in% prot_map$entrez)

# Map PPI EntrezA and EntrezB to protein ID.
idx <- match(ppiDat$musEntrezA, prot_map$entrez)
ppiDat$protA <- prot_map$prot[idx]

idx <- match(ppiDat$musEntrezB, prot_map$entrez)
ppiDat$protB <- prot_map$prot[idx]

ppiDat$ida <- paste(ppiDat$protA, ppiDat$protB, sep = "_")
ppiDat$idb <- paste(ppiDat$protB, ppiDat$protA, sep = "_")

# Create edge list
edge_list <- melt(adjm)
edge_list$id <- paste(edge_list$Var1, edge_list$Var2,sep="_")
edge_list$ppi <- as.numeric(edge_list$id %in% ppiDat$ida | edge_list$id %in% ppiDat$idb)

# Create new ppi adjm.
ppi_adjm <- matrix(edge_list$ppi, nrow = 3022, ncol = 3022)
colnames(ppi_adjm) <- rownames(ppi_adjm) <- colnames(adjm)

# Write to file as adjm.
# Cannot have column or row names for proper import as python igraphs!:w
write.csv(adjm, "bicor_adjm.csv")
write.csv(ppi_adjm, "ppi_adjm.csv")

# ENDOFILE
#------------------------------------------------------------------------------
