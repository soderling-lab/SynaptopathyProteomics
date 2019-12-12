#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# Best resolution may be resolution 100, but at this resolution there are no
# significant Striatum modules (KW p.adj[bonferroni] <0.05).
# Relaxing threshold for significance by using method = "BH" and alpha = 0.1,
# then there are some sig modules...

# User parameters to change:
net <- "Cortex" # Which network are we analyzing?
#mypart <- c(Cortex = "10360847",Striatum = "10342568")[net]
mypart <- c(Cortex = "10773682",Striatum = "10781799")[net] # relaxed criterion

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(WGCNA)
  library(fgsea)
  library(reactome.db)
})

# Directories.
if (rstudioapi::isAvailable()) {
	setwd("D:/projects/SynaptopathyProteomics/analysis/3_Network-Analysis")
}
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load expression data.
myfiles <- list(
  Cortex = file.path(rdatdir, "3_Cortex_cleanDat.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_cleanDat.RData")
)
data <- lapply(myfiles, readRDS)
# Data should be transposed: rows, proteins.
data <- lapply(data, t)[[net]]

# Load Sample info.
traits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load correlation (adjacency) matrices.
myfiles <- list(
  Cortex = file.path(rdatdir, "3_Cortex_Adjm.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_Adjm.RData")
)
adjm <- lapply(myfiles, readRDS)[[net]]

# Load network partitions-- self-preservation enforced.
# myfile <- list.files(rdatdir, pattern = "1023746", full.names = TRUE) # WT and KO
#relaxed_criterion <- c("10773682","10781799")
# myfile <- list.files(rdatdir, pattern= "Combined_Module",full.names=TRUE) # Combined network only
myfiles <- list(
  Cortex = list.files(rdatdir, pattern = "10360847", full.names = TRUE),
  Striatum = list.files(rdatdir, pattern = "10342568", full.names = TRUE)
)
partitions <- lapply(myfiles, readRDS)[[net]]

#------------------------------------------------------------------------------
## Perform GSE analysis.
#------------------------------------------------------------------------------

# Compile pathways from reactome.db.
message("Compiling mouse pathways from the Reactome database.")
suppressMessages({
pathways <- reactomePathways(protmap$entrez)
})

# Loop to perform GSE analysis.
message(paste("Analyzing",net,"modules for gene set enrichment.","\n"))
GSEresults <- list()
for (i in seq_along(partitions)){
# Get partition.
partition <- partitions[[i]]
message(paste("Working on resolution:", i,"..."))
# Get Modules.
modules <- split(partition, partition)
names(modules) <- paste0("M", names(modules))
# Number of modules.
nModules <- sum(names(modules) != "M0")
# Power does not influence MEs.
MEdata <- moduleEigengenes(data,
  colors = partition,
  softPower = 1, impute = FALSE
)
MEs <- as.matrix(MEdata$eigengenes)
# Calculate module membership (kME).
KMEdata <- signedKME(data, MEs, corFnc = "bicor")
# Loop through modules, perform GSEA.
moduleGSE <- list()
for (m in seq_along(modules)){
  # Calculate ranks as KME.
  # Another metric for ranks by be signed Fold change * -log10pvalue.
  prots <- names(modules[[m]])
  subKME <- subset(KMEdata,rownames(KMEdata) %in% prots)
  colnames(subKME) <- names(modules)
  ranks <- subKME[,m]
  names(ranks) <- prots
  # Map proteins to entrez.
  entrez <- protmap$entrez[match(names(ranks), protmap$ids)]
  names(ranks) <- entrez
  # Ranks must be sorted in increasing order.
  ranks <- ranks[order(ranks)]
  # Perform GSEA.
  suppressWarnings({
  GSEdata <- fgsea(pathways, ranks, minSize=5, maxSize=500, nperm = 1e+05)
  })
  # Sort by p-value.
  GSEdata <- GSEdata[order(GSEdata$pval), ]
  # Add column for gene symbols to results.
  genes <- lapply(GSEdata$leadingEdge, function(x) {
			  protmap$gene[match(x, protmap$entrez)]
			  })
  GSEdata$Symbols <- genes
  moduleGSE[[m]] <- GSEdata
  } # Ends inner loop.
  # Status.
  nSig <- sum(sapply(moduleGSE,function(x) any(x$padj<0.05)))
  message(paste("... Modules with any significant GSE:",nSig,
		"of",nModules,"(",round(100*nSig/nModules,2),"%)","\n"))
# Return GSE results for 
GSEresults[[i]] <- moduleGSE
}

# Save results.
myfile <- file.path(rdatdir,paste0("3_",net,"_Module_GSE_Results.RData")) 
saveRDS(GSEresults, myfile)

quit()

#------------------------------------------------------------------------------
## Examine GO results in order to define ~best resolution.
#------------------------------------------------------------------------------

# Examine biological enrichment of modules at every resolution.
# Summarize the biological significance of a resolution as the sum of 
# -log(GO pvalues) for all modules.
modSig <- lapply(GOresults, function(x) 
		 sapply(x, function(y) sum(-log(y$pValue))))

modSig <- lapply(GOresults, function(x) 
		 sapply(x, function(y) mean(-log(y$pValue))))

# The code above is confusing, this is what it does:
#x = results[[1]] # list of go enrichment for all modules at res 1.
#y = x[[1]] # go enrichment df of module 1.
#y$pValue
#-log(y$pValue)
#sum(-log(y$pValue))
#sapply(x,function(y) sum(-log(y$pValue)))
#lapply(results, function(x) sapply(x,function(y) sum(-log(y$pValue))))
#out = lapply(results, function(x) sapply(x,function(y) sum(-log(y$pValue))))

# Summarize every resolution.
resSum <- sapply(modSig, sum)
best_res <- c(1:length(resSum))[resSum == max(resSum)]
names(best_res) <- net
# Status report. 
message(paste("Best resolution based on GO enrichment:",best_res))

# Save results.
myfile <- file.path(rdatdir,paste0("3_",net,"_Best_Resolution.RData"))
saveRDS(best_res, myfile)
