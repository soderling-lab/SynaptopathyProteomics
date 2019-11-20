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

# Load statistical results.
glmDat <- readRDS(file.path(rdatdir,"2_GLM_Results.RData"))
rownames(glmDat) <- protmap$ids[match(glmDat$Uniprot,protmap$uniprot)]

# Collect sig prots.
sig <- apply(glmDat[,c(2:ncol(glmDat))],1,function(x) any(x<0.05))
sigProts <- rownames(glmDat)[sig]

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
# Module sig prot enrichment.
#------------------------------------------------------------------------------

### SKIP ### 
freq <- length(unique(sigProts))/dim(protmap)[1]

ratio <- list()
for (i in 1:length(comparisons)){
data <- comparisons[[i]]
modules <- split(data$wtProts,data$wtPartition)
changes <- sapply(modules,unique)
keep <- names(changes)[changes == "divergent"]
if (length(keep)>0){
modules <- modules[names(modules)[names(modules) %in% keep]]
observed <- sapply(modules,function(x) sum(names(x) %in% sigProts))
expected <- sapply(modules, length) * freq
ratio[[i]] <- observed/expected
} else {
	ratio[[i]] <- 0
}
}

c(1:100)[sapply(ratio,max)==max(sapply(ratio,max))]

c(1:100)[sapply(ratio,median)==max(sapply(ratio,median))]
# Resolution with best enrichment of sigprots in WT modules = 88
# Resolution with best enrichment of sigprots in KO modules = 95

#------------------------------------------------------------------------------
# Prepare ppi graph.
#------------------------------------------------------------------------------

# Load mouse interactome.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
idx <- musInteractome$Interactor_A_Taxonomy %in% c(10090, 9606, 10116)
ppis <- subset(musInteractome, idx)

prots <- colnames(wtAdjm)
entrez <- protmap$entrez[match(prots,protmap$ids)]

# Build a graph with all proteins.
graph <- buildNetwork(ppis, entrez, taxid = 10090)

#------------------------------------------------------------------------------
# Examine ~optimal resolution.
#------------------------------------------------------------------------------

# Best resolution based on average go sem sim of divergent modules (sum).
# WT = 90
# KO = 89

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
koR <- 89

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
getAdjm <- function(module,adjm){
	idx <- idy <- match(names(module),rownames(adjm))
	       subAdjm <- adjm[idx,idy]
	       return(subAdjm)
}
subWT <- lapply(wtDivergent,function(x) getAdjm(x,wtAdjm))
subKO <- lapply(wtDivergent,function(x) getAdjm(x,koAdjm))

# Mean edge strength
sapply(subKO,mean)
sapply(subWT,mean)

# rewired proteins.
n <- 1
df <- melt(subKO[[n]]-subWT[[n]])
df$abs <- abs(df$value)
df$wt <- melt(subWT[[n]])$value
df$ko <- melt(subKO[[n]])$value
df <- df[order(df$abs,decreasing=TRUE),]
head(df)

# But are they connected?
prots <- names(wtDivergent[[3]])
entrez <- protmap$entrez[match(prots, protmap$ids)]
subg <- induced_subgraph(graph,entrez)
length(V(subg))
length(E(subg))
sum(prots %in% sigProts)

prots

#------------------------------------------------------------------------------
# Examine KO.
#------------------------------------------------------------------------------

data <- comparisons[[koR]]
koModules <- split(data$koProts,data$koPartition)
koChanges <- sapply(koModules,unique)
koDivergent <- koModules[names(koChanges[koChanges=="divergent"])]
message(paste("Number of KO divergent modules:",length(koDivergent)))
names(koDivergent)

# Sizes of divergent modules.
sapply(koDivergent,length)

# Check correlation coefficients.
subWT <- lapply(koDivergent,function(x) getAdjm(x,wtAdjm))
subKO <- lapply(koDivergent,function(x) getAdjm(x,koAdjm))

# Mean edge strengths.
sapply(subWT,mean)
sapply(subKO,mean)

# rewired proteins.
n <- 4
df <- melt(subKO[[n]]-subWT[[n]])
df$abs <- abs(df$value)
df$wt <- melt(subWT[[n]])$value
df$ko <- melt(subKO[[n]])$value
df <- df[order(df$abs,decreasing=TRUE),]
head(df)

# Collect GO data.
koGO <- koAllGO[[koR]]
names(koGO) <- names(koModules) # for ease, fix names.

# Check top go.
lapply(koGO[names(koChanges[koChanges=="divergent"])],function(x){
	       cbind(x$shortDataSetName[c(1:5)],x$FDR[c(1:5)])
})

# But are they connected?
# Not really...
prots <- names(koDivergent[[1]])
entrez <- protmap$entrez[match(prots, protmap$ids)]
subg <- induced_subgraph(graph,entrez)
length(V(subg))
length(E(subg))
sum(prots %in% sigProts)
prots


