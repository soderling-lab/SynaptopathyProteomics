#!/usr/bin/env Rscript

# Examine the WGCNA partition.

# Load the expression data.

#------------------------------------------------------------------------------
## Prepare the workspace.
#------------------------------------------------------------------------------

# Directories.
here <- getwd()
rootdir <- dirname(dirname(here))
datadir <- file.path(rootdir,"data")
tabsdir <- file.path(rootdir,"tables")
figsdir <- file.path(rootdir, "figures")
funcdir <- file.path(rootdir,"functions")

# Global options and imports.
suppressPackageStartupMessages({
	library(data.table)
	library(WGCNA)
	library(ggplot2)
        library(grid)
	library(gridExtra)
})

# Load functions.
source(file.path(funcdir,"clean_fun.R"))

# Load the normalized expression data.
data <- readRDS(file.path(datadir,"wtDat.Rds"))

# Load the optimized parameters.
parameters <- fread(file.path(rootdir,"data","1_optimized_params.txt"))
colnames(parameters) <- c("Parameter","Value")

# Load the partition.
partition <- fread(file.path(rootdir,"data","partition1.txt"))
colnames(partition) <- c("Protein","Module")

# Check:
all(partition$Protein == colnames(data))

#------------------------------------------------------------------------------
## Examine nModules, parameters, other basics...
#------------------------------------------------------------------------------

# Examine number of modules, and numbers of proteins per module.
tab <- as.data.frame(table(partition$Module))
colnames(tab) <- c("Module","nProteins")
tab <- tableGrob(tab, rows = NULL)
ggsave(file.path(figsdir,"0_wgcna_table.tiff"), tab, 
		 width = grobsize(tab)[1], height = ggsize(tab)[2])

# Examine optimized parameters.
tab1 <- tableGrob(parameters, rows = NULL)
ggsave(file.path(figsdir,"1_wgcna_parameters.tiff"), tab1, 
		 width = grobsize(tab1)[1], height = grobsize(tab1)[2])

#------------------------------------------------------------------------------
## Examine heirarchical clustering.
#------------------------------------------------------------------------------

# Calculate TOM dissimilarity matrix, excluding grey.
idy <- partition$Module != "grey"
diss <- 1 - TOMsimilarityFromExpr(
		    datExpr     = data[,idy],
		    corType     = "bicor",
		    networkType = "signed",
		    power       = 13,
		    TOMType     = "signed",
		    TOMDenom    = "min",
		    verbose     = 1
				  )
#colnames(diss) <- rownames(diss) <- colnames(data[,idy])

# Heirarchical clustering.
hc <- hclust(as.dist(diss), method = "average")

# Generate plot.
tiff(file.path(figsdir,"2_dendrogram.tiff"))
plotDendroAndColors(hc, 
		    partition$Module[idy], 
		    "Dynamic Tree Cut", 
		    rowText = NULL,
		    dendroLabels = FALSE,
		    hang = 0.03, 
		    addGuide = TRUE,
		    guideHang = 0.05, 
		    main = "Heirarchical Clustering Dendrogram"
		    )
dev.off() 

#------------------------------------------------------------------------------
##
#------------------------------------------------------------------------------
