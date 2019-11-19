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
  library(anRichment)
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
myfile <- list.files(rdatdir,pattern="6142226",full.names=TRUE)
partitions <- readRDS(myfile)

# Function to perform GO enrichment for all modules in a given partition.
getModuleGO <- function(partitions,geno,resolution,protmap,musGOcollection){
	part <- partitions[[resolution]]
	modules <- split(part[[geno]],part[[geno]])
	dm <- sapply(c(1:length(modules)), function(x) part[[geno]]==x)
	colnames(dm) <- paste0("R",resolution,"-M",names(modules))
	logic <- dm == TRUE
	for (i in 1:ncol(dm)){
		col_header <- colnames(dm)[i]
		dm[logic[,i],i] <- col_header
		dm[!logic[,i],i] <- "NA"
	}
	# Prots mapped to entrez.
	entrez <- protmap$entrez[match(rownames(dm),protmap$ids)]
	# Perform GO enrichment.
	GOenrichment <- enrichmentAnalysis(
					   classLabels = dm,
					   identifiers = entrez,
					   refCollection = musGOcollection,
					   useBackground = "given",
					   threshold = 0.05,
					   thresholdType = "Bonferroni",
					   getOverlapEntrez = TRUE,
					   getOverlapSymbols = TRUE,
					   ignoreLabels = "FALSE"
					   verbose = 1
					   )
	# Collect the results.
	GO_results <- list()
	for (r in 1:length(GOenrichment$setResults)) {
		GO_results[[r]] <- GOenrichment$setResults[[r]]$enrichmentTable
	}
	names(GO_results) <- colnames(dm)
	return(GO_results)
}

# Perform WT GO enrichment.
wtGO <- sapply(c(1:100),function(x) getModuleGO(partitions,"wt",x,protmap,musGOcollection))

# Perform WT GO enrichment.
koGO<- sapply(c(1:100),function(x) getModuleGO(partitions,"ko",x,protmap,musGOcollection))
