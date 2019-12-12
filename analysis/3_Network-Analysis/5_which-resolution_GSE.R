#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# GSE significance seems to decline with increasing resolution...

# User parameters to change:
net <- "Cortex" # Which network are we analyzing?
# mypart <- c(Cortex = "10360848",Striatum = "10342568")[net] # all 7 stats.
mypart <- c(Cortex = "10773682", Striatum = "10781799")[net] # relaxed criterion

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
  setwd("D:/projects/SynaptopathyProteomics/analysis/4_Network-Analysis")
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
myfile <- file.path(rdatdir, paste0("3_", net, "_Adjm.RData"))
adjm <- readRDS(myfile)

# Load network partitions-- self-preservation enforced.
myfile <- list.files(rdatdir, pattern = mypart, full.names = TRUE)
partitions <- readRDS(myfile)

#------------------------------------------------------------------------------
## Perform GSE analysis.
#------------------------------------------------------------------------------

# Compile pathways from reactome.db.
message(paste("Preparing to analyze", net, "modules for gene set enrichment..."))
message(paste("... Compiling mouse pathways from the Reactome database.", "\n"))
suppressMessages({
  pathways <- reactomePathways(protmap$entrez)
})

# Loop to perform GSE analysis.
GSEresults <- list()
for (i in seq_along(partitions)) {
  # Get partition.
  partition <- partitions[[i]]
  message(paste("Working on resolution:", i, "..."))
  # Get Modules, remove M0.
  modules <- split(partition, partition)
  names(modules) <- paste0("M", names(modules))
  modules <- modules[names(modules) != "M0"]
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
  modSig <- list()
  modSig2 <- list()
  modSig3 <- list()
  for (m in seq_along(modules)) {
    # Calculate ranks as KME.
    # Another metric for ranks by be signed Fold change * -log10pvalue.
    prots <- names(modules[[m]])
    subKME <- subset(KMEdata, rownames(KMEdata) %in% prots)
    colnames(subKME) <- names(modules)
    ranks <- subKME[, m]
    names(ranks) <- prots
    # Map proteins to entrez.
    entrez <- protmap$entrez[match(names(ranks), protmap$ids)]
    names(ranks) <- entrez
    # Ranks must be sorted in increasing order.
    ranks <- ranks[order(ranks)]
    # Perform GSEA.
    suppressWarnings({
      GSEdata <- fgsea(pathways, ranks, minSize = 5, maxSize = 500, nperm = 1e+05)
    })
    # Sort by p-value.
    GSEdata <- GSEdata[order(GSEdata$pval), ]
    modSig[[m]] <- sum(-log(GSEdata$pval)) # Sum of -log(pvals) for the module.
    modSig2[[m]] <- sum(-log(GSEdata$pval) * GSEdata$ES) # Sum of the product...
    modSig3[[m]] <- sum(-log(GSEdata$pval) * GSEdata$NES) # Normalized enrichment score...
    moduleGSE[[m]] <- GSEdata # GSE table for the module.
  } # Ends inner loop.
  # Fix names. Sum module significance for the resolution.
  names(moduleGSE) <- names(modules)
  names(modSig) <- names(modules)
  names(modSig2) <- names(modules)
  names(modSig3) <- names(modules)
  resSig <- log2(sum(unlist(modSig)))
  resSig2 <- sum(unlist(modSig2))
  resSig3 <- sum(unlist(modSig3))
  # Status.
  nSig <- sum(sapply(moduleGSE, function(x) any(x$padj < 0.05)))
  message(paste0(
    "... Modules with any significant GSE: ", nSig,
    " of ", nModules, " (", round(100 * (nSig / nModules), 2), "%)"
  ))
  message(paste("... Resolution GSE significance score:", round(resSig, 4)))
  message(paste("... Alternative significance score #1:", round(resSig2, 4)))
  message(paste("... Alternative significance score #2:", round(resSig3, 4), "\n"))
  # Return GSE results for
  GSEresults[[i]] <- moduleGSE
}

# Save results.
myfile <- file.path(rdatdir, paste0("3_", net, "_Module_GSE_Results.RData"))
saveRDS(GSEresults, myfile)

#------------------------------------------------------------------------------
## Examine GO results in order to define ~best resolution.
#------------------------------------------------------------------------------

# Examine biological enrichment of modules at every resolution.
# Summarize the biological significance of a resolution as the sum of
# -log(GO pvalues) for all modules.
modSig <- lapply(GSEresults, function(x) {
  sapply(x, function(y) sum(-log(y$pval)))
})
names(modSig) <- paste0("R", seq_along(partitions))

# Insure that any list elements with length 0 are removed.
keep <- seq_along(modSig)[sapply(modSig, function(x) length(x) != 0)]
modSig <- modSig[keep]

## The code above is confusing, this is what it does step-by-step:
# x = results[[1]] # list of go enrichment for all modules at res 1.
# y = x[[1]] # go enrichment df of module 1.
# y$pValue
#-log(y$pValue)
# sum(-log(y$pValue))
# sapply(x,function(y) sum(-log(y$pValue)))
# lapply(results, function(x) sapply(x,function(y) sum(-log(y$pValue))))
# out = lapply(results, function(x) sapply(x,function(y) sum(-log(y$pValue))))

# Summarize every resolution.
resSum <- sapply(modSig, sum)
best_res <- names(resSum)[resSum == max(resSum)]
names(best_res) <- net

# Status report.
message(paste("Best resolution based on GSE:", best_res))

# Save results.
myfile <- file.path(rdatdir, paste0("3_", net, "_GSE_Best_Resolution.RData"))
saveRDS(best_res, myfile)
