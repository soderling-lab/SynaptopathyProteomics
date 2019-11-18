#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(dendextend)
  library(getPPIs) 
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Load protein map.
protmap <- readRDS(file.path(rdatdir,"2_Prot_Map.RData"))

# Load mouse PPIs.
data(musInteractome)

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

# Load network partitions.
myfile <- list.files(rdatdir,pattern="1023746",full.names=TRUE)
partitions <- readRDS(myfile)

# Load network comparison results.
myfile <- list.files(rdatdir,pattern="6040093",full.names=TRUE)
comparisons <- readRDS(myfile)

# Combine output into df.
res <- data.frame(
		"nWT" =  unlist(lapply(output, function(x) sum(names(table(x$wtPartition))!="0"))),
		"nKO" = unlist(lapply(output, function(x) sum(names(table(x$wtPartition))!="0"))),
		"nPresWT" = unlist(lapply(output, function(x) sum(x$wtProts == "preserved"))),
		"nDivWT" = unlist(lapply(output, function(x) sum(x$wtProts == "divergent"))),
		"nPresKO" = unlist(lapply(output, function(x) sum(x$koProts == "preserved"))),
		"nDivKO" = unlist(lapply(output, function(x) sum(x$koProts == "divergent"))))
res$percentTotalDivergence <- (res$nDivWT + res$nDivKO)/(2*2918)
res$percentTotalPreservation <- (res$nPresWT+res$nPresKO)/(2*2918)
df <- list()
for (i in 1:length(output)){
	x <- output[[i]]
	df[[i]] <- data.frame(
			 k_WT_NS = sum(sapply(split(x$wtProts,x$wtPartition),unique)=="ns"),
			 k_WT_Di = sum(sapply(split(x$wtProts,x$wtPartition),unique)=="divergent"),
			 k_WT_Pr = sum(sapply(split(x$wtProts,x$wtPartition),unique)=="preserved"),
			 k_WT_NC = sum(sapply(split(x$wtProts,x$wtPartition),unique)=="not-clustered"),
			 k_KO_NS = sum(sapply(split(x$koProts,x$koPartition),unique)=="ns"),
			 k_KO_Di = sum(sapply(split(x$koProts,x$koPartition),unique)=="divergent"),
			 k_KO_Pr = sum(sapply(split(x$koProts,x$koPartition),unique)=="preserved"),
			 k_KO_NC = sum(sapply(split(x$koProts,x$koPartition),unique)=="not-clustered")
			 )
}
df <- do.call(rbind,df)
df <- cbind(res,df)

# Examine a partition.
res <- 37
dat <- output[[res]]

# Modules.
modules <- split(dat$koProts,dat$koPartition)
sapply(modules, unique)

# Divergent modules.
myprots <- names(modules[[4]])
getEntrez <- function(prots) {return(protmap$entrez[match(prots,protmap$ids)])}
g <- buildNetwork(hitpredict=musInteractome,getEntrez(myprots),taxid=10090)

idx <- idy <- match(myprots,colnames(wtAdjm))
subWT <- wtAdjm[idx,idy]
subKO <- koAdjm[idx,idy]

