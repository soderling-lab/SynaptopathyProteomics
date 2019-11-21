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
  library(igraph)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

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
myfile <- list.files(rdatdir,pattern="6490667",full.names=TRUE)
comparisons <- readRDS(myfile)

# Load module divergence df.
moduleChanges <- fread(list.files(rdatdir,pattern="Divergence.csv",full.names=TRUE))

# Fix column names.
moduleChanges <- setNames(moduleChanges,gsub(" ","_",colnames(moduleChanges)))

# Which resolutions have changes?
x <- moduleChanges %>% filter(ko_nDivergent == 1) %>% 
	dplyr::select(resolution)
#dput(as.numeric(x$resolution))
#c(29, 35, 36, 40, 41, 42, 44, 45, 48, 49, 55, 58, 66, 79)

# What proteins are these?
#getProts <- function(comparisons)
subComp <- comparisons[c(29, 35, 36, 40, 41, 42, 44, 45, 48, 49, 55, 58, 66, 79)]
namen <- rep(NA,length(subComp))

# Collect proteins in divergent ko modules.
modProts <- list()
for (i in 1:length(subComp)){
	prots = subComp[[i]]$koProts
	parts = subComp[[i]]$koPartition
	modules = split(prots,parts)
	changes = sapply(modules,unique)
	namen[i] <- names(changes)[changes=="divergent"]
	modProts[[i]] <-names(modules[[namen[i]]])
}
names(modProts) <- namen

# How do these groups of proteins relate to each other?
# All possible combinations...
x = expand.grid(1:14,1:14)
dm <- matrix(ncol=14,nrow=14)
for (i in 1:dim(x)[1]){
	idx <- x[i,1]
	idy <- x[i,2]
	int_ <- sum(modProts[[idx]] %in% modProts[[idy]])
	union_ <- length(unique(c(modProts[[idx]],modProts[[idy]])))
	ji <- int_/union_
	dm[idx,idy] <- ji
}

# Convert to distance matrix and then plot.


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

#c(29, 35, 36, 40, 41, 42, 44, 45, 48, 49, 55, 58, 66, 79)
res <- 58 # best GO = 29, biggest with with dysregulate module = 79
GO <- koAllGO
geno <- "ko"
###
#examinePartition <- function(comparisons,res,geno,GO){
allDat <- comparisons[[res]]
namen <- names(allDat)[grep(geno,names(allDat))]
data <- allDat[namen]
modules <- split(data[[2]],data[[1]])
nModules <- length(modules)
changes <- sapply(modules,unique)
divergent <- modules[names(changes[changes=="divergent"])]
goDat <- GO[[res]]
myfile <- file.path(tabsdir,paste0("3_",toupper(geno),"_R",res,"_Module_GO_enrichment.xlsx"))
write_excel(goDat,myfile)
# for ease, fix names.
names(goDat) <- names(modules) 
# Sizes of divergent modules.
nDivergent <- sapply(divergent,length) # Proteins.
message(paste("Divergent module name :",names(divergent)))
message(paste("Number of divergent modules:",length(divergent)))
message(paste("Number of proteins in divergent module:",length(unlist(divergent))))
#}

# Check correlation coefficients.
getAdjm <- function(module,adjm){
	idx <- idy <- match(names(module),rownames(adjm))
	       subAdjm <- adjm[idx,idy]
	       return(subAdjm)
}
subWT <- lapply(divergent,function(x) getAdjm(x,wtAdjm))
subKO <- lapply(divergent,function(x) getAdjm(x,koAdjm))
# Mean edge strength
sapply(subKO,mean)
sapply(subWT,mean)

# But are they connected?
prots <- names(divergent[[1]])
entrez <- protmap$entrez[match(prots, protmap$ids)]
subg <- induced_subgraph(graph,entrez)
length(V(subg))
length(E(subg))

# How many sig?
sum(prots %in% sigProts)

# From which genos?
subdat <- melt(subset(glmDat,rownames(glmDat) %in% prots),id.vars="Uniprot")
subdat %>% group_by(variable) %>% summarize(sum(value<0.05))
