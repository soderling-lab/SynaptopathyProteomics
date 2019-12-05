#!/usr/bin/env Rscript

# Analyze combined network partitions. Which resolution to pick???

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(fgsea)
  library(getPPIs)
  library(purrr)
  library(WGCNA)
  library(org.Mm.eg.db)
  library(anRichment)
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

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Prot_Map.RData"))

# Load statistical results.
glm_results <- readRDS(file.path(rdatdir,"2_Combined_All_GLM_Results.RData"))

# Load expression data.
data <- t(readRDS(file.path(rdatdir,"3_Combined_cleanDat.RData")))

# Load correlation matrix.
combAdjm <- t(readRDS(file.path(rdatdir,"3_Combined_Adjm.RData")))

# Load network partitions.
#myfile <- list.files(rdatdir, pattern = "6142226", full.names = TRUE)
myfile <- list.files(rdatdir, pattern = "9109552", full.names = TRUE)
partitions <- readRDS(myfile)

#-------------------------------------------------------------------------------
## Collect GLM statistics in a list.
#-------------------------------------------------------------------------------

# Names of relevant columns.
colNames <- colnames(glm_results[[1]])[c(2,5:9)]
stats <- colNames[c(2:length(colNames))]

# Collect data.
subDat <- lapply(glm_results,function(x) x[,colNames])

# Combine into a single df.
df <- subDat  %>% reduce(left_join, by="Uniprot")

# Rename columns.
newNames <- paste(rep(names(glm_results),each=length(colNames)-1),
		  sapply(strsplit(colnames(df)[c(2:ncol(df))],"\\."),"[",1))
colnames(df)[c(2:ncol(df))] <- newNames

# Collect each statistic into a single df in a list.
glm_stats <- sapply(stats,function(x) df[,c(1,grep(x,colnames(df)))])

# Clean up data a little...
glm_stats <- lapply(glm_stats,function(x) { 
			    x <- x[order(x$Uniprot),]
			    idx <- match(x$Uniprot,protmap$uniprot)
			    rownames(x) <- protmap$ids[idx]
			    x$Uniprot <- NULL 
			    return(x)}
)

#-------------------------------------------------------------------------------
## Unpack the permutation results.
#-------------------------------------------------------------------------------

# Resolutions.
resolutions <- c(1:length(partitions))

# Collect combined partitions.
partitions <- lapply(partitions,function(x) x$combined)

# Number of modules. Subtract one for NS modules (not preserved).
nModules <- sapply(partitions,function(x) length(unique(x)))-1

# Collect modules from each partition.
modules <- lapply(partitions,function(x) split(x,x))

# Remove "0" modules.
filtModules <- lapply(modules, function(x) x[-c(1:length(x))[names(x) == "0"]])

# Percent NS.
percentNS <- sapply(partitions,function(x) sum(x==0)/length(x))

# Module summary expression profiles (eigen vectors).
MEs <- lapply(partitions,function(x) moduleEigengenes(data,x,impute=FALSE,
						      excludeGrey=FALSE,
						      softPower=1,
					              verbose=0)[[1]])


# Module coherence.
PVE <- sapply(resolutions,function(x) propVarExplained(data, 
						       partitions[[x]], 
						       MEs[[x]], 
						       corFnc = "bicor"))

# Fix module names.
newNames <- lapply(partitions,function(x) paste0("M",names(table(x))))
# Function to rename PVE and ME lists.
renameList <- function(myList,namesList) {
	myList <- lapply(c(1:length(myList)), function(x) {
				 names(myList[[x]]) <- namesList[[x]]
				 return(myList[[x]])
})
	return(myList)
} # Ends function.

MEs <- renameList(MEs,newNames)
PVE <- renameList(PVE,newNames)

# Remove M0.
MEs <- lapply(MEs,function(x) { x[,"M0"] <- NULL; return(x) })
PVE <- lapply(PVE,function(x) x[-1])

# Mean percent variance explained.
meanPVE <- sapply(PVE,function(x) mean(x))
medianPVE <- sapply(PVE,function(x) median(x))
maxPVE <- sapply(PVE,function(x) max(x))

#------------------------------------------------------------------------------
## Perform GO analysis of modules at every resolution.
#------------------------------------------------------------------------------

# Load mouse GO collection.
myfile <- list.files(rdatdir, "musGO", full.names = TRUE)
musGO <- readRDS(myfile)

# Function to perform GO enrichment for all modules in a given partition.
getModuleGO <- function(partitions, resolution, protmap, musGOcollection) {
  part <- partitions[[resolution]]
  modules <- split(part,part)
  dm <- sapply(names(modules), function(x) part == x)
  colnames(dm) <- paste0("R", resolution, "-M", names(modules))
  logic <- dm == TRUE
  for (i in 1:ncol(dm)) {
    col_header <- colnames(dm)[i]
    dm[logic[, i], i] <- col_header
    dm[!logic[, i], i] <- "FALSE"
  }
  # Prots mapped to entrez.
  entrez <- protmap$entrez[match(rownames(dm), protmap$ids)]
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
    ignoreLabels = "FALSE",
    verbose = 0
  )
  # Collect the results.
  GO_results <- list()
  for (r in 1:length(GOenrichment$setResults)) {
    GO_results[[r]] <- GOenrichment$setResults[[r]]$enrichmentTable
  }
  names(GO_results) <- colnames(dm)
  return(GO_results)
} # Ends function.

# Loop to perform GO enrichment for modules at every resolution. 
message(paste("Evaluating GO enrichment of WT modules at every resolution!", "\n"))
n <- length(partitions) # n resolutions.
results <- list()
for (i in seq_along(resolutions)) {
	# Initialize progress bar.
	if (i == 1) { 
		pb <- txtProgressBar(min = 0, max = n, style = 3)
	}
	# Perform GO analysis.
	results[[i]] <- getModuleGO(partitions,resolution=i,protmap,musGO)
	# Update progress bar.
	setTxtProgressBar(pb, i)
	if (i == n) {
		# Close pb, save.
		close(pb)
	        myfile <- file.path(rdatdir, "3_Module_GO_Results.RData")
	        saveRDS(results, myfile)
	        message("Done!")
	}
} # Ends loop.

# Examine results.

# Need to pick a resolution.... analyze module coherence...
