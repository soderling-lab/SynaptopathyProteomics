#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change:
net <- "Cortex" # Which network are we analyzing?
# mypart <- c(Cortex = "10360848",Striatum = "10342568")[net] # all 7 stats.
mypart <- c(Cortex = "10773682", Striatum = "10781799")[net] # relaxed criterion

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(WGCNA)
  library(getPPIs)
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
# Data should be transposed: rows, proteins.
myfile <- file.path(rdatdir, paste0("3_", net, "_cleanDat.RData"))
data <- t(readRDS(myfile))

# Load Sample info.
traits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load correlation (adjacency) matrices.
myfile <- file.path(rdatdir, paste0("3_", net, "_Adjm.RData"))
adjm <- readRDS(myfile)

# Load network partitions-- self-preservation enforced.
myfile <- list.files(rdatdir, pattern = mypart, full.names = TRUE)
partitions <- readRDS(myfile)

#------------------------------------------------------------------------------
## Build gene set for GSEA.
#------------------------------------------------------------------------------

# Load GO dataset from BROAD Instititute. These lists are human entrez ids.
# We will map human genes to mouse using the getHomologs() function from
# GetPPIs.

# Broad datasets are from: http://software.broadinstitute.org/gsea/downloads.jsp
# Ge lab datasets are from: http://ge-lab.org/#/data
# Bader lab datasets are from: http://download.baderlab.org/EM_Genesets/current_release/
datasets <- c(
  All_curated = "c2.all.v7.0.entrez.gmt", # 1
  Hallmark = "h.all.v7.0.entrez.gmt", # 2
  All_canonical = "c2.cp.v7.0.entrez.gmt", # 3
  ALL_GO = "c5.all.v7.0.entrez.gmt", # 4
  GO_MF = "c5.mf.v7.0.entrez.gmt", # 5
  GO_BP = "c5.bp.v7.0.entrez.gmt", # 6
  GO_CC = "c5.cc.v7.0.entrez.gmt", # 7
  MSigDB = "msigdb.v7.0.entrez.gmt", # 8. This is all BROAD lists!
  Reactome = "c2.cp.reactome.v7.0.entrez.gmt", # 9 Human Reactome.
  TF_motif = "c3.tft.v7.0.entrez.gmt", # 10
  microRNA_motif = "c3.mir.v7.0.entrez.gmt", # 11
  All_motif = "c3.all.v7.0.entrez.gmt", # 12
  Ge_mus_curated = "mGSKB_Entrez.gmt", # 13
  Ge_mus_compiled = "Ge_Mus_GO_KEGG.gmt", # 14. Compiled Mouse GO/KEGG
  Bader = "Mouse_AllPathways_December_01_2019_entrezgene.gmt"
)

# Choose a dataset.
dataset <- datasets[3]
pathways <- gmtPathways(file.path(rdatdir, dataset))
# Remove paths with missing (NA) names.
keep <- c(1:length(pathways))[!is.na(names(pathways))]
pathways <- pathways[keep]
nPaths <- length(pathways)
# Get mouse homologs of human genes in pathways list.
hsEntrez <- unique(unlist(pathways))
msHomologs <- getHomologs(hsEntrez, taxid = "10090") # mouse taxid
names(msHomologs) <- hsEntrez
# Map to mouse, discard unmapped (NA) genes.
msCanonical <- lapply(pathways, function(x) purrr::discard(msHomologs[x], is.na))

# Load SynGO pathways.
myfile <- file.path(rdatdir, "SynGO_Pathways.RData")
syngo <- do.call(c, readRDS(myfile)) # Combine BP and CC.

# Combine Syngo and canonical pathways.
msPathways <- c(msCanonical, syngo, reactome)

# Check: What percentage of genes are in pathways list?
genes <- protmap$entrez
percentMapped <- sum(genes %in% unique(unlist(msCanonical))) / length(genes)
message(paste("Percent genes in pathways:", round(100 * percentMapped, 2)))

# Total number of pathways.
message(paste("Number of pathways:", length(msPathways)))

# Load SynGO pathways.
myfile <- file.path(rdatdir, "SynGO_Pathways.RData")
syngo <- do.call(c, readRDS(myfile)) # Combine BP and CC.

# Compile pathways from reactome.db.
message(paste("Preparing to analyze", net, "modules for gene set enrichment..."))
message(paste("... Compiling mouse pathways from the Reactome database.", "\n"))
suppressMessages({
  reactome <- reactomePathways(protmap$entrez)
})

# Combine SynGO and Reactome pathways.
pathways <- c(syngo, reactome)

# Number of pathways:
nPaths <- length(pathways)

# How many proteins in Synaptosome are mapped to pathways?
nProts <- length(protmap$entrez)
percent_mapped <- sum(protmap$entrez %in% unlist(pathways)) / nProts
message(paste(
  "Percent proteins mapped to compiled pathways:",
  round(100 * percent_mapped, 2), "(%)"
))

#------------------------------------------------------------------------------
## Perform GSE analysis.
#------------------------------------------------------------------------------

# Which pathways to use?
msPathways <- c(msCanonical,reactome,syngo)

# Loop to perform GSE analysis for all modules at every resolution...
GSEresults <- list()
for (i in seq_along(partitions)) {
  # Get partition.
  partition <- partitions[[i]]
  message(paste("Working on resolution:", i, "..."))
  # Get Modules.
  modules <- split(partition, partition)
  names(modules) <- paste0("M", names(modules))
  # Remove M0.
  modules <- modules[names(modules) != "M0"]
  # Number of modules.
  nModules <- sum(names(modules) != "M0")
  # Power does not influence MEs.
  MEdata <- moduleEigengenes(data,
    colors = partition,
    softPower = 1, impute = FALSE,
    excludeGrey = TRUE
  )
  MEs <- as.matrix(MEdata$eigengenes)
  # Calculate module membership (kME).
  KMEdata <- signedKME(data, MEs, corFnc = "bicor", outputColumnName = "")
  colnames(KMEdata) <- names(modules)
  # Clean up KME table... Add column with entrez ids.
  KMEdata <- add_column(KMEdata, Protein = rownames(KMEdata), .before = 1)
  Entrez <- protmap$entrez[match(rownames(KMEdata), protmap$ids)]
  KMEdata <- add_column(KMEdata, Entrez, .after = 1)
  # Add column with module assignment.
  KMEdata <- add_column(KMEdata, Module = partition[rownames(KMEdata)], .after = 2)
  # Empty lists for collection GSE results.
  moduleGSE <- list()
  modSig <- list()
  # Loop through modules, perform GSEA.
  for (m in seq_along(modules)) {
    # Calculate ranks as KME.
    prots <- rownames(KMEdata)
    ranks <- KMEdata[, names(modules)[m]]
    names(ranks) <- prots
    # Map proteins to entrez.
    entrez <- protmap$entrez[match(names(ranks), protmap$ids)]
    names(ranks) <- entrez
    # Ranks must be sorted in increasing order.
    ranks <- ranks[order(ranks)]
    # Perform GSEA.
    suppressWarnings({
      GSEdata <- fgsea(msPathways, ranks, minSize = 5, maxSize = 500, nperm = 1e+05)
    })
    # Sort by p-value.
    GSEdata <- GSEdata[order(GSEdata$pval), ]
    modSig[[m]] <- sum(-log(GSEdata$pval) * GSEdata$NES) # -log(Pval) * NES - normalized ES.
    moduleGSE[[m]] <- GSEdata # GSE table for the module.
  } # Ends inner loop.
  # Fix names. Sum module significance for the resolution.
  names(moduleGSE) <- names(modules)
  names(modSig) <- names(modules)
  sigScore <- log2(sum(unlist(modSig)))
  # Status.
  message(paste("... Resolution GSE significance score:", round(resSig, 4)))
  # Return GSE results for
  GSEresults[[i]] <- list("moduleGSE" = moduleGSE,"KMEdata" = KMEdata, "Score" = sigScore)
}

# Does the result make sense?
m = split(partitions[[1]],partitions[[1]])
prots = names(m[[1]])
gseDat <- moduleGSE[[1]]
subdat <- subset(gseDat, gseDat$ES > 0 & gseDat$padj < 0.05)
subdat <- subdat[order(subdat$pval, decreasing = TRUE), ]
moduleGenes <- protmap$entrez[match(prots, protmap$ids)]
subdat$PercentInModule <- sapply(subdat$leadingEdge,function(x) sum(x %in% moduleGenes)/length(x))
fwrite(subdat,"foo.csv",row.names=TRUE)

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
