#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change:
net <- "Cortex" # Which network are we analyzing?
# mypart <- c(Cortex = "10360848",Striatum = "10342568")[net] # Preservation enforced using all 7 stats.
mypart <- c(Cortex = "10773682", Striatum = "10781799")[net] # Relaxed preservation criterion.

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

# Load functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load expression data.
# Data should be transposed: rows, proteins.
myfile <- file.path(rdatdir, paste0("3_", net, "_cleanDat.RData"))
data <- t(readRDS(myfile))

# Load network partitions-- self-preservation enforced.
myfile <- list.files(rdatdir, pattern = mypart, full.names = TRUE)
partitions <- readRDS(myfile)

#------------------------------------------------------------------------------
## Build gene set for GSEA.
#------------------------------------------------------------------------------

# Load GSE datasets from BROAD Instititute. These lists are human entrez ids.
# We will map human genes to mouse homologs using the getHomologs() function 
# from my GetPPIs library.

## 1. Load Broad Institute gene lists:
# Broad datasets are from: http://software.broadinstitute.org/gsea/downloads.jsp
message(paste("Preparing to analyze", net, "modules for gene set enrichment..."))
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
  All_motif = "c3.all.v7.0.entrez.gmt" # 12
)

# Choose a Broad dataset.
message(paste("... Mapping curated human pathways from the Broad",
	     "Institute to mouse."))
dataset <- datasets[1] # All BROAD curated pathways.
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

## 2. Load SynGO data:
message(paste("... Compiling mouse SynGO pathways."))
myfile <- file.path(rdatdir, "SynGO_Pathways.RData")
msSyngo <- do.call(c, readRDS(myfile)) # Combine BP and CC.

## 3. Compile mouse pathways from reactome.db.
message(paste("... Compiling mouse Reactome pathways.", "\n"))
suppressMessages({
  msReactome <- fgsea::reactomePathways(protmap$entrez)
})

# Combine SynGO and Reactome pathways.
msPathways <- c(msCanonical, msSyngo, msReactome)

# Check: What percentage of genes are in pathways list?
genes <- protmap$entrez
percentMapped <- sum(genes %in% unique(unlist(msPathways))) / length(genes)
message(paste("... Percent genes in compiled pathways:", round(100 * percentMapped, 2)))

# Total number of pathways.
message(paste("... Total number of pathways:", length(msPathways)))

#------------------------------------------------------------------------------
## Perform GSE analysis.
#------------------------------------------------------------------------------

## Outstanding questions?
# Which pathways to use? --> Synprot cononical pathways and msReactome.
# Should NS Proteins be removed? --> Don't calculate GSE for M0.
# Should GSE be performed on module genes or all genes? --> All genes.
# Should GSE of proteins outside module be penalized? --> Use KME to rank genes.

# Loop to perform GSE analysis for all modules at every resolution...
GSEresults <- list()
subsetProts <- FALSE
message("Perfomring GSE analysis for all modules at every resolution...")
for (i in c(1,20,40,60,80,100)) {
#for (i in seq_along(partitions)) {
  # Get partition.
  partition <- partitions[[i]]
  message(paste("Working on resolution:", i, "..."))
  # Get modules for given partition. Exclude M0.
  modules <- split(partition, partition)
  names(modules) <- paste0("M", names(modules))
  modules <- modules[names(modules)[!names(modules) == "M0"]]
  # Number of modules.
  nModules <- length(modules)
  message(paste("... Number of modules:",nModules))
  # Calculate module eigengenes (MEs).
  MEdata <- moduleEigengenes(data,
    colors = partition,
    softPower = 1, impute = FALSE,
    excludeGrey = TRUE
  )
  MEs <- as.matrix(MEdata$eigengenes)
  # Calculate signed module membership (kME).
  KMEdata <- signedKME(data, MEs, corFnc = "bicor", outputColumnName = "")
  colnames(KMEdata) <- names(modules)
  # Clean-up kME table: Add protein id, entrez id, and module assignment cols.
  KMEdata <- add_column(KMEdata, Protein = rownames(KMEdata), .before = 1)
  Entrez <- protmap$entrez[match(rownames(KMEdata), protmap$ids)]
  KMEdata <- add_column(KMEdata, Entrez, .after = 1)
  KMEdata <- add_column(KMEdata, Module = partition[rownames(KMEdata)], .after = 2)
  # Empty lists for collection GSE results.
  moduleGSE <- list()
  modSignificance <- list()
  # Loop through modules, perform GSEA.
  for (m in seq_along(modules)) {
    # Calculate ranks as signed module membership (KME).
	  if (subsetProts) {
		  prots <- names(modules[[m]]) # just module proteins...
	  } else {
		  prots <- rownames(KMEdata) # all proteins...
	  }
    ranks <- KMEdata[prots, names(modules)[m]]
    names(ranks) <- prots
    # Map proteins to entrez.
    entrez <- protmap$entrez[match(names(ranks), protmap$ids)]
    names(ranks) <- entrez
    # Ranks must be sorted in increasing order.
    ranks <- ranks[order(ranks)]
    # Perform GSEA.
    suppressWarnings({
      GSEdata <- fgsea(msPathways, ranks,
        minSize = 5,
        maxSize = 500, nperm = 1e+05
      )
    })
    # Calculate module significance score as:
    # -log(Pval) * abs(NES) | NES = normalized enrichment score
    modSignificance[[m]] <- sum(-log(GSEdata$pval) * abs(GSEdata$NES))
    moduleGSE[[m]] <- GSEdata # GSE table for the module.
  } # Ends loop.
  # Fix names. Sum module significance for the resolution.
  names(moduleGSE) <- names(modules)
  # Number of modules with significant enrichment.
  nSigMods <- sum(sapply(moduleGSE,function(x) any(x$padj<0.05 & x$NES>0)))
  message(paste("... Number of modules with any significant enrichment:",nSigMods))
  # Summarize resolution as the log2 of sum module significance scores.
  sigScore <- log2(sum(unlist(modSignificance)))
  message(paste("... Summary GSE significance score:", round(sigScore, 4),"\n"))
  # Return GSE results for
  GSEresults[[i]] <- list(
    "moduleGSE" = moduleGSE,
    "KMEdata" = KMEdata,
    "Score" = sigScore
  )
  # Save results.
  if (i == 100) {
    myfile <- file.path(rdatdir, paste0("3_", net, "_Module_GSE_Results.RData"))
    saveRDS(GSEresults, myfile)
  }
}

# Sometimes proteins KME is high, but it is not assigned to that module... 
# how to handle this senario... 
output = GSEresults
x = output[[60]]
write_excel(x$moduleGSE,"temp.xlsx")
fwrite(x$KMEdata,"kme.csv")

#------------------------------------------------------------------------------
## Examine GSE results in order to define ~best resolution.
#------------------------------------------------------------------------------

# Status report.
#message(paste("Best resolution based on GSE:", best_res))

# Save results.
#myfile <- file.path(rdatdir, paste0("3_", net, "_GSE_Best_Resolution.RData"))
#saveRDS(best_res, myfile)
