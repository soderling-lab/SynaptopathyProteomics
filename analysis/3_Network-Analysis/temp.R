#!/usr/bin/env Rscript
# Examine module self-preservation.

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(WGCNA)
  library(NetRep)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
datadir <- file.path(root, "rdata")
funcdir <- file.path(root, "R")

# Functions.
myfun <- list.files(funcdir, pattern = "silently.R", full.names = TRUE)
invisible(sapply(myfun, source))

# Load expression data. Transpose -> rows = samples; columns = genes.
wtDat <- t(readRDS(file.path(datadir, "3_WT_cleanDat.RData")))
koDat <- t(readRDS(file.path(datadir, "3_KO_cleanDat.RData")))

# Compute adjmatrix:
wtAdjm <- silently(WGCNA::bicor(wtDat))
koAdjm <- silently(WGCNA::bicor(koDat))

# Compute TOM adjcacency matrices--this insures that all edges are positve.
wtTOM <- TOMsimilarity(wtAdjm, TOMType = "signed", verbose = 0)
koTOM <- TOMsimilarity(koAdjm, TOMType = "signed", verbose = 0)
rownames(wtTOM) <- colnames(wtTOM) <- colnames(wtAdjm)
rownames(koTOM) <- colnames(koTOM) <- colnames(koAdjm)

# Load partitions.
myfiles <- list.files(datadir, pattern = "*partitions.csv", full.names = TRUE)
koParts <- data.table::fread(myfiles[1], drop = 1, skip = 1)
wtParts <- data.table::fread(myfiles[2], drop = 1, skip = 1)
colnames(koParts) <- colnames(wtParts) <- colnames(wtAdjm)

x = list()
for (i in 1:100){
	wt = as.numeric(as.data.frame(wtParts)[i,])+1
	ko = as.numeric(as.data.frame(koParts)[i,])+1
	names(wt) <- names(ko) <- colnames(wtParts)
	x[[i]] = list(wt=wt,ko=ko)
}
myfile <- file.path(datadir,"3_la_partitions.RData")
saveRDS(x,myfile)
