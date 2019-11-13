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
source_myfun <- function() {
  myfun <- list.files(funcdir, pattern = ".R", full.names = TRUE)
  invisible(sapply(myfun, source))
}
source_myfun()

# Store all plots in list.
all_plots <- list()
ggtheme()

# Load expression data.
wtDat <- t(readRDS(list.files(rdatdir, pattern = "WT_cleanDat", full.names = TRUE)))
koDat <- t(readRDS(list.files(rdatdir, pattern = "KO_cleanDat", full.names = TRUE)))

# Load adjmatrices.
wtAdjm <- t(readRDS(list.files(rdatdir, pattern = "WT_Adjm.RData", full.names = TRUE)))
koAdjm <- t(readRDS(list.files(rdatdir, pattern = "KO_Adjm.RData", full.names = TRUE)))

# Load preserved partitions of co-expression graph:
myfile <- list.files(rdatdir, pattern = "preservation", full.names = TRUE)
partitions <- readRDS(myfile)

# Look at resolution 44.
wtParts <- partitions[[44]][['wt']]
koParts <- partitions[[44]][['ko']]

# Modules.
wtModules <- split(wtParts,wtParts)
koModules <- split(koParts,koParts)

# Load the module changes.
data <- readRDS(list.files(datadir,pattern="5595360",full.names=TRUE))
subdat <- data[[44]]

# Modules with changes.
wtDiff <- subdat$Changes$ko
koDiff <- subdat$Changes$wt
names(wtDiff) <- names(wtModules)[-1]
names(koDiff) <- names(koModules)[-1]

# koDiff.
names(koDiff)[koDiff=="divergent"]

# Divergent KO modules (GOF).
prots <- list(names(koModules[["4"]])
	      ,names(koModules[["5"]]),
	      names(koModules[["11"]]))

# Load protein map.
protmap <- readRDS(file.path(rdatdir,"2_Prot_Map.RData"))

# Get entrez genes for prots of interest.
getEntrez <- function(prots,protmap){
	return(protmap$entrez[match(prots,protmap$ids)])
}
entrez <- lapply(prots,function(x) getEntrez(x,protmap))

# Build ppi graph.
# Load mouse interactome data.
data(musInteractome)
g <- lapply(entrez,function(x) buildNetwork(musInteractome,x,taxid=10090))


