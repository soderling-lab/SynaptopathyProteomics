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
	require(WGCNA)
	require(igraph)
	require(getPPIs)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatdir <- file.path(root, "rdata")
funcdir <- file.path(root, "R")

# Functions.
myfun <- list.files(funcdir, pattern = "silently.R", full.names = TRUE)
invisible(sapply(myfun, source))

# Load expression data. Transpose -> rows = samples; columns = genes.
wtDat <- t(readRDS(file.path(rdatdir, "3_WT_cleanDat.RData")))
koDat <- t(readRDS(file.path(rdatdir, "3_KO_cleanDat.RData")))

# Compute adjmatrix:
wtAdjm <- silently(WGCNA::bicor(wtDat))
koAdjm <- silently(WGCNA::bicor(koDat))

# Calculate power for approximate scale free fit.
sft <- silently({
	sapply(list(wtDat,koDat),function(x) 
	       pickSoftThreshold(x,
				 corFnc="bicor",
				 networkType="signed",
				 RsquaredCut=0.8)$powerEstimate)
})
names(sft) <- c("wt","ko")

# All proteins, map to entrez.
prots <- colnames(wtAdjm)

# Protein id map.
myfile <- list.files(rdatdir,"Map",full.names=TRUE)
protmap <- readRDS(myfile) 

# Get Entrez.
entrez <- protmap$entrez[match(prots,protmap$ids)]

# Build a network.
data("musInteractome")
g <- buildNetwork(musInteractome,entrez,species=c(10090,9606,10116),taxid=10090)

# Degree connectivity.
dc <- apply(as.matrix(as_adjacency_matrix(g)),2,sum)

# Evaluate scale free fit.
fit <- scaleFreeFitIndex(dc, nBreaks =10, removeFirst = FALSE)

# Load statistical results.
myfile <- list.files(rdatdir,pattern="GLM",full.names=TRUE)
data <- readRDS(myfile)


