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

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(getPPIs)
  library(purrr)
  library(WGCNA)
  library(org.Mm.eg.db)
  library(anRichment)
  library(getPPIs)
  library(DescTools)
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

# Load GLM stats.
myfile <- file.path(rdatdir, "2_GLM_Stats.RData")
glm_stats <- readRDS(myfile)

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
##
#------------------------------------------------------------------------------

i <- 1

# Get partition.
partition <- partitions[[i]]
message(paste("Working on resolution:", i,"..."))
# Get Modules.
modules <- split(partition, partition)
names(modules) <- paste0("M", names(modules))
# Number of modules.
nModules <- sum(names(modules) != "M0")
message(paste0("... Number of modules: ", nModules))
# Percent not clustered.
percentNC <- sum(partition == 0) / length(partition)
message(paste("... Percent of proteins not clustered:", round(100 * percentNC, 2), "(%)"))
# Calculate Module eigengenes.
# Power does not influence MEs.
MEdata <- moduleEigengenes(data,
  colors = partition,
  softPower = 1, impute = FALSE
)
MEs <- as.matrix(MEdata$eigengenes)
# Get Percent Variance explained (PVE)
PVE <- MEdata$varExplained
names(PVE) <- names(modules)
meanPVE <- mean(as.numeric(PVE[names(PVE) != "M0"]))
message(paste("... Mean module coherence (PVE):", round(100 * meanPVE, 2), "(%)"))
# Create list of MEs.
ME_list <- split(MEs, rep(1:ncol(MEs), each = nrow(MEs)))
names(ME_list) <- names(modules)
# Calculate module membership (kME).
KMEdata <- signedKME(data, MEs, corFnc = "bicor")

library(fgsea)
data(examplePathways)

m = 1

# Loop through modules, perform GSEA.
moduleGO <- list()

for (m in seq_along(modules)){

  # Calculate ranks as KME.
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
  # Another metric for ranks by be signed Fold change * -log10pvalue.

  library(reactome.db)

  my_pathways <- reactomePathways(names(exampleRanks))

  # Perform GSEA.
  suppressWarnings({
  GSEdata <- fgsea(examplePathways, ranks, minSize=5, maxSize=15, nperm = 1e+05)
  })

  # Sort by p-value.
  GSEdata <- GSEdata[order(GSEdata$pval), ]


# plot the most significantly enriched pathway
  library(ggplot2)
  plotEnrichment(examplePathways[[head(GSEdata[order(pval), ], 1)$pathway]],
		                exampleRanks) + labs(title=head(GSEdata[order(pval), ], 1)$pathway)

  topPathwaysUp <- GSEdata[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- GSEdata[ES < 0][head(order(pval), n=10), pathway]

  # Add column for gene symbols to results.
  genes <- lapply(GSEdata$leadingEdge, function(x) {
			  protmap$gene[match(x, protmap$entrez)]
			  })
  GSEdata$Symbols <- genes
  moduleGO[[m]] <- GSEdata
  } # Ends inner loop.
