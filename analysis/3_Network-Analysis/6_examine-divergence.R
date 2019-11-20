#!/usr/bin/env Rscript

#' ---
#' title:
#' description: 
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(WGCNA)
  library(NetRep)
  library(getPPIs)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein id map.
myfile <- list.files(rdatdir, "Map", full.names = TRUE)
protmap <- readRDS(myfile)


# Load GO results.
myfiles <- list.files(rdatdir,pattern="Module_GO_Results",full.names=TRUE)
koAllGO <- readRDS(myfiles[1])
wtAllGO <- readRDS(myfiles[2])

# Load expression data.
wtDat <- t(readRDS(list.files(rdatdir,
  pattern = "WT_cleanDat",
  full.names = TRUE
)))
koDat <- t(readRDS(list.files(rdatdir,
  pattern = "KO_cleanDat",
  full.names = TRUE
)))

# Load adjmatrices.
wtAdjm <- t(readRDS(list.files(rdatdir,
  pattern = "WT_Adjm.RData",
  full.names = TRUE
)))
koAdjm <- t(readRDS(list.files(rdatdir,
  pattern = "KO_Adjm.RData",
  full.names = TRUE
)))

# Calculate power for approximate scale free fit.
sft <- silently({
	sapply(list(wtDat,koDat),function(x) 
	       pickSoftThreshold(x,
				 corFnc="bicor",
				 networkType="signed",
				 RsquaredCut=0.8)$powerEstimate)
})
names(sft) <- c("wt","ko")

# Load network comparison results.
myfile <- list.files(rdatdir,pattern="6171865",full.names=TRUE)
comparisons <- readRDS(myfile)

#------------------------------------------------------------------------------
# Prepare ppi graph.
#------------------------------------------------------------------------------

# Load mouse interactome.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
idx <- musInteractome$Interactor_A_Taxonomy %in% c(10090, 9606, 10116)
ppis <- subset(musInteractome, idx)

# Build a graph with all proteins.
graph <- buildNetwork(ppis, entrez, taxid = 10090)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# Best resolution based on number of divergent modules.
# WT = 90
# KO = 92

# Best resolution based on percent proteins assigned to divergent modules:
# R = 47

# Best resolution based on total number of divergent modules.
# R = 92

# Best resolution based on sum GO p for all modules.
# WT = 52 # problem: no wt divergent.
# KO = 50

# Best resolution based on sum GO p for DIVERGENT modules.
# WT = 42 - single divergent modules, top go = icpm for both wt and ko.
# KO = 33 - module seems to be the SAME proteins...

# Best resolution based on GO semantic similarity of modules.
# WT = 97, 100 (modularity) - problem no sig KO GO!?
# KO = 99, 100 (modularity)

# Look at "optimal" resolution.
wtR <- 90
koR <- 92

# Examine WT.
data <- comparisons[[wtR]]
wtModules <- split(data$wtProts,data$wtPartition)
wtChanges <- sapply(wtModules,unique)
wtDivergent <- wtModules[names(wtChanges[wtChanges=="divergent"])]
message(paste("Number of WT divergent modules:",length(wtDivergent)))
wtGO <- wtAllGO[[wtR]]
names(wtGO) <- names(wtModules) # for ease, fix names.

# Sizes of divergent modules.
sapply(wtDivergent,length)

# Check top go.
lapply(wtGO[names(wtChanges[wtChanges=="divergent"])],function(x){
	       cbind(x$shortDataSetName[c(1:5)],x$FDR[c(1:5)])
})

# Check correlation coefficients.
wtSubAdjm <- sapply(wtDivergent,function(x){
	       idx <- idy <- match(names(x),rownames(wtAdjm))
	       wtSubAdjm <- wtAdjm[idx,idy]
	       idx <- idy <- match(names(x),rownames(koAdjm))
	       koSubAdjm <- koAdjm[idx,idy]
	       return(mean(koSubAdjm)-mean(wtSubAdjm))
})

x = wtSubAdjm[[1]]
wt = x[[1]]
ko = x[[2]]
mean(wt)
mean(ko)

sum(abs(wtSubAdjm))

sapply(wtSubAdjm,mean) # High average edge weight.

# Write to file.
write_excel(wtGO,"wtGO.xlsx")

#------------------------------------------------------------------------------
# Examine KO.
data <- comparisons[[koR]]
koModules <- split(data$koProts,data$koPartition)
koChanges <- sapply(koModules,unique)
koDivergent <- koModules[names(koChanges[koChanges=="divergent"])]
message(paste("Number of KO divergent modules:",length(koDivergent)))
names(koDivergent)

# Sizes of divergent modules.
sapply(koDivergent,length)

# Collect GO data.
koGO <- koAllGO[[koR]]
names(koGO) <- names(koModules) # for ease, fix names.

# Check top go.
lapply(koGO[names(koChanges[koChanges=="divergent"])],function(x){
	       cbind(x$shortDataSetName[c(1:5)],x$FDR[c(1:5)])
})

# Check correlation coefficients.
koSubAdjm <- lapply(koDivergent,function(x){
	       idx <- idy <- match(names(x),rownames(koAdjm))
	       subAdjm <- koAdjm[idx,idy]
})
sapply(koSubAdjm,mean) # High average edge weight.

# Write GO to file.
write_excel(koGO,"koGO.xlsx")

# But are they connected?
# Not really...
prots <- names(wtDivergent[[3]])
entrez <- protmap$entrez[match(prots, protmap$ids)]
subg <- induced_subgraph(graph,entrez)
length(V(subg))
length(E(subg))
