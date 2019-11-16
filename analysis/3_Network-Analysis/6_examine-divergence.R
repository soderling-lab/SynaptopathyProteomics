#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
## Prepare the workspace.
#------------------------------------------------------------------------------

# Directories.
here <- getwd()
rootdir <- dirname(dirname(here))
datadir <- file.path(rootdir, "data")
rdatdir <- file.path(rootdir, "rdata")
tabsdir <- file.path(rootdir, "tables")
figsdir <- file.path(rootdir, "figures")
funcdir <- file.path(rootdir, "R")

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(org.Mm.eg.db)
  library(anRichment)
  library(readxl)
  library(getPPIs)
})

# Load required custom functions.
myfun <- list.files(funcdir, pattern = ".R", full.names = TRUE)
invisible(sapply(myfun, source))

# Store all plots in list.
all_plots <- list()

# ggplot theme
ggtheme()

# Load expression data.
wtDat <- t(readRDS(list.files(rdatdir, pattern = "WT_cleanDat", full.names = TRUE)))
koDat <- t(readRDS(list.files(rdatdir, pattern = "KO_cleanDat", full.names = TRUE)))

# Load adjmatrices.
wtAdjm <- t(readRDS(list.files(rdatdir, pattern = "WT_Adjm.RData", full.names = TRUE)))
koAdjm <- t(readRDS(list.files(rdatdir, pattern = "KO_Adjm.RData", full.names = TRUE)))

# Load preserved partitions of co-expression graph:
myfile <- list.files(rdatdir, pattern = "5716254", full.names = TRUE)
partitions <- readRDS(myfile)

# Which resolution?
res <- 23
wtParts <- partitions[[res]][["wt"]]
koParts <- partitions[[res]][["ko"]]

# Modules.
wtModules <- split(wtParts, wtParts)
koModules <- split(koParts, koParts)
nWT <- length(wtModules)
nKO <- length(koModules)

# Load protein map.
protmap <- readRDS(file.path(rdatdir, "2_Prot_Map.RData"))

# Get entrez genes for prots of interest.
getEntrez <- function(prots, protmap) {
  return(protmap$entrez[match(names(prots), protmap$ids)])
}
wtEntrez <- lapply(wtModules, function(x) getEntrez(x, protmap))
koEntrez <- lapply(koModules, function(x) getEntrez(x, protmap))

# Build ppi graph.
# Load mouse interactome data.
data(musInteractome)
wtGraphs <- lapply(wtEntrez, function(x) buildNetwork(musInteractome, x, taxid = 10090))
koGraphs <- lapply(koEntrez, function(x) buildNetwork(musInteractome, x, taxid = 10090))

#

