#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# User params.
save = FALSE # Save or Load goSemSim graph.

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(dendextend)
  library(getPPIs)
  library(anRichment)
  library(org.Mm.eg.db)
  library(utils)
  library(GOSemSim)
  library(igraph)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Load protein map.
protmap <- readRDS(file.path(rdatdir, "2_Prot_Map.RData"))

# Load network comparison results.
myfile <- list.files(rdatdir, pattern = "6490667", full.names = TRUE)

comparisons <- readRDS(myfile)

# Load GO data.
# Use Molecular function ontology.
msGO <- godata("org.Mm.eg.db", ont = "MF",computeIC=FALSE)

#-------------------------------------------------------------------------------
# Create GO semantic similarity network.
#-------------------------------------------------------------------------------

# Map prots to entrez.
prots <- names(comparisons[[1]]$wtProts)
entrez <- protmap$entrez[match(prots, protmap$ids)]

if (save) {
# Calculate Go sem sim.
goSim <- mgeneSim(entrez,semData = msGO,measure = "Wang",verbose = TRUE)
# Save data.
saveRDS(goSim, file.path(rdatdir, "3_All_GO_SemSim.RData"))
} else {
# Load data.
goSim <- readRDS(file.path(rdatdir, "3_All_GO_SemSim.RData"))
}

# For ease, map rows and column names back to prot ids.
rownames(goSim) <- protmap$ids[match(rownames(goSim),protmap$entrez)]
colnames(goSim) <- protmap$ids[match(colnames(goSim),protmap$entrez)]

#-------------------------------------------------------------------------------
# Which partition is best? Evaluate how much divergent modules have changed.
#-------------------------------------------------------------------------------

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

# Loop through all comparisons...
pb <- txtProgressBar(min = 0, max = length(comparisons), style = 3)
delta <- rep(NA,100)

for (i in 1:length(comparisons)){
    setTxtProgressBar(pb, i)
	wtModules <- split(comparisons[[i]][["wtProts"]],comparisons[[i]][["wtPartition"]])
	koModules <- split(comparisons[[i]][["koProts"]],comparisons[[i]][["koPartition"]])
	wtModules <- wtModules[!names(wtModules) == "0"]
	koModules <- koModules[!names(koModules) == "0"]
	wtDelta <- sapply(wtModules,function(x) {
		       idx <- idy <- match(names(x),colnames(wtAdjm))
		       delta <- mean(koAdjm[idx,idy] - wtAdjm[idx,idy])
})
	koDelta <- sapply(koModules,function(x) {
		       idx <- idy <- match(names(x),colnames(koAdjm))
		       delta <- mean(wtAdjm[idx,idy] - koAdjm[idx,idy])
})
	delta[i] <- sum(abs(wtDelta),abs(koDelta))
	# Close pb.
	if (i == length(comparisons)){
		close(pb)
		message(" Complete!\n")
	}
}

#-------------------------------------------------------------------------------
# Which partition is best?
#-------------------------------------------------------------------------------

# 1. average edge weight of modules.
# Best partition has modules with maximum average edge strenght. 
# ~ Biological cohesiveness.

# Function to calculate average edge weight of modules.
getAvgWeight <- function(goSim,comparisons,myProts,myPartition){
# Make a graph.
g <- graph_from_adjacency_matrix(goSim,weighted=TRUE,"undirected")
avg_weight <- list()
# Loop through all comparisons, caluculate avg edge weight of modules.
  pb <- txtProgressBar(min = 0, max = length(comparisons), style = 3)
for (i in 1:length(comparisons)){
    setTxtProgressBar(pb, i)
modules <- split(comparisons[[i]][[myProts]],comparisons[[i]][[myPartition]])
	avg_weight[[i]] <- sapply(modules,function(x) {
		       v <- names(x)
		       keep <- v %in% names(V(g))
		       v <- v[keep]
	               subg <- induced_subgraph(g,v)
		       return(mean(edge_attr(subg, "weight")))
})
}
close(pb)
return(avg_weight)
}

# Calculate average intra-module go sem sim.
wtWeight <- getAvgWeight(goSim,comparisons,"wtProts","wtPartition")

koWeight <- getAvgWeight(goSim,comparisons,"koProts","koPartition")

# Which partition is best?
# Calculate average edge weight of all modules, exlude module 0.
wt <- sapply(wtWeight,function(x) mean(x[!names(x)==0]))
ko <- sapply(koWeight,function(x) mean(x[!names(x)==0]))

best_r <- c(1:100)[wt==max(wt)]
best_r # WT #97

best_r <- c(1:100)[ko==max(ko)]
best_r # KO #99

# which partition is best?
# Calculate average edge weight of divergent modules!
getChanges <- function(comparisons,myProts,myPartition){
	changes <- lapply(comparisons,function(x) {
				  sapply(split(x[[myProts]],x[[myPartition]]),unique)
})
	return(changes)
}

# Changes.
wtChanges <- getChanges(comparisons,"wtProts","wtPartition")
koChanges <- getChanges(comparisons,"koProts","koPartition")

out <- list()
for (i in 1:length(comparisons)){
	idx = names(wtChanges[[i]])[wtChanges[[i]]=="divergent"]
	if (length(idx)>0){
		wt <- wtWeight[[i]][idx]
	} else {
		wt <- 0
	}
	idy = names(koChanges[[i]])[koChanges[[i]]=="divergent"]
	if (length(idy)>0){
		ko <- koWeight[[i]][idy]
	} else {
		ko <- 0
	}
	out[[i]] <- list("wt"=wt,"ko"=ko)
}

fx <- "sum"
wt <- lapply(out,function(x) x$wt)
wtAvg <- sapply(wt,eval(fx))
c(1:100)[wtAvg==max(wtAvg)] # 90
ko <- lapply(out,function(x) x$ko)
koAvg <- sapply(ko,eval(fx))
c(1:100)[koAvg==max(koAvg)] #89

#-------------------------------------------------------------------------------
# Which partition is best?
#-------------------------------------------------------------------------------

# 2. Modularity? -- will this be biased towards large modules?
# Best partition maximizes modularity of GO sem sim graph.
# ~ Biological cohesiveness.
# Remove 0 index modules = nodes that are not-clustered. 

# Function to calculate modularity. 
# Moderately time consuming...
getModularity <- function(goSim,comparisons,myProts,myPartition){
	# Make a graph.
	g <- graph_from_adjacency_matrix(goSim,weighted=TRUE,"undirected")
	q <- rep(NA,length(comparisons))
	# Loop through all comparisons, caluculate avg edge weight of modules.
	pb <- txtProgressBar(min = 0, max = length(comparisons), style = 3)
	for (i in 1:length(comparisons)){
		# Init progress bar.
		setTxtProgressBar(pb, i)
		# Get partition.
		v <- comparisons[[i]][[myPartition]]
		# Remove prots from membership vector that are not in go sem graph.
		keep <- names(v) %in% names(V(g))
		v <- v[keep]
		# Remove 0 index modules from graph.
		out <- v==0
		v <- v[!out]
		subg <- induced_subgraph(g,v)
		# Sort membership vector.
		v <- v[match(names(v),names(V(subg)))]
		# Calculate modularity.
		q[i] <- modularity(subg, membership=v, weights = edge_attr(g,"weight"))
	}
	close(pb)
	return(q)
}

# Calculate modularity of GO sem sim graph.
# NOTE: Zero index modules are excluded!
wtQ <- getModularity(goSim,comparisons,"wtProts","wtPartition")

koQ <- getModularity(goSim,comparisons,"koProts","wtPartition")

# Best partitions:
r_best <- c(1:100)[wtQ==max(wtQ)]
r_best # WT # 100

r_best <- c(1:100)[koQ==max(koQ)]
r_best # KO # 100
