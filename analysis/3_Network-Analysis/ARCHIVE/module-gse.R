#!/usr/bin/env Rscript

#' ---
#' title: 
#' description: 
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
## Prepare the workspace.
#-------------------------------------------------------------------------------

# Imports 
suppressPackageStartupMessages({
  library(dplyr)
  library(WGCNA)
  library(getPPIs)
  library(fgsea)
})


# Directories.
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

# Load statistical results.
#glm_results <- readRDS(file.path(rdatdir, "2_Combined_All_GLM_Results.RData"))

# Load GLM stats.
#myfile <- file.path(rdatdir, "2_GLM_Stats.RData")
#glm_stats <- readRDS(myfile)

# Load expression data.
cleanDat <- readRDS(file.path(rdatdir, "3_Cortex_cleanDat.RData"))
data <- t(cleanDat)
colnames(data) <- rownames(cleanDat)

# Load Sample info.
#traits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load network partitions-- self-preservation enforced.
myfile <- list.files(rdatdir,pattern="10360847",full.names=TRUE) # Cortex
#myfile <- list.files(rdatdir, pattern="10342568",full.names=TRUE) # Striatum
partitions <- readRDS(myfile)

#------------------------------------------------------------------------------
## Prepare a GSE dataset.
#------------------------------------------------------------------------------

# Load GO dataset from BRAD Instititute. We will map human genes to mouse.
# From: http://software.broadinstitute.org/gsea/downloads.jsp
# Ge lab data from: http://ge-lab.org/#/data
# Baderlab data: http://download.baderlab.org/EM_Genesets/current_release/
  datasets <- c(
    All_curated = "c2.all.v7.0.entrez.gmt", # 1
    Hallmark = "h.all.v7.0.entrez.gmt", # 2
    All_canonical = "c2.cp.v7.0.entrez.gmt", # 3
    ALL_GO = "c5.all.v7.0.entrez.gmt", # 4
    GO_MF = "c5.mf.v7.0.entrez.gmt", # 5
    GO_BP = "c5.bp.v7.0.entrez.gmt", # 6
    GO_CC = "c5.cc.v7.0.entrez.gmt", # 7
    MSigDB = "msigdb.v7.0.entrez.gmt", # 8. This is all lists!
    Reactome = "c2.cp.reactome.v7.0.entrez.gmt", # 9
    TF_motif = "c3.tft.v7.0.entrez.gmt", # 10
    microRNA_motif = "c3.mir.v7.0.entrez.gmt", # 11
    All_motif = "c3.all.v7.0.entrez.gmt", # 12
    Ge_mus_curated = "mGSKB_Entrez.gmt", # 13
    Ge_mus_compiled = "Ge_Mus_GO_KEGG.gmt", # 14. Compiled GO/KEGG
    Bader = "Mouse_AllPathways_December_01_2019_entrezgene.gmt"
  )

  # Choose a dataset.
  dataset <- datasets[13]
  pathways <- gmtPathways(file.path(rdatdir, dataset))
  filter <- TRUE # Should genes that are not in synaptic proteome be removed?
  map2mouse <- FALSE # Map human genes to mouse?

  # Clean up pathways...
  keep <- c(1:length(pathways))[!is.na(names(pathways))]
  pathways <- pathways[keep]
  nPaths <- length(pathways)

  # Get mouse homologs of human genes in pathways list.
  if (map2mouse) {
    hsEntrez <- unique(unlist(pathways))
    msHomologs <- getHomologs(hsEntrez, taxid = "10090") # mouse taxid
    names(msHomologs) <- hsEntrez
    # Map to mouse, discard unmapped (NA) genes.
    msPathways <- lapply(pathways, function(x) purrr::discard(msHomologs[x], is.na))
  } else {
    msPathways <- pathways
  }

  # Filter pathways; keep genes in data; remove empty pathways;
  genes <- protmap$entrez[match(colnames(data), protmap$ids)]
  filter_pathways <- function(pathways, genes) {
    pathways <- lapply(pathways, function(x) purrr::discard(x, x %notin% genes))
    keep <- names(pathways)[!sapply(pathways, function(x) length(x) == 0)]
    pathways <- pathways[keep]
    return(pathways)
  }
  if (filter) {
    msPathways <- filter_pathways(msPathways, genes)
  }

  # Check: What percentage of genes are in pathways list?
  percentMapped <- sum(genes %in% unique(unlist(msPathways))) / length(genes)
  message(paste("Percent genes in pathways:", round(100 * percentMapped, 2)))
  # Total number of pathways.
  message(paste("Number of pathways:", length(msPathways)))

#-------------------------------------------------------------------------------
## Perform module gene set enrichment analysis (GSEA)
#-------------------------------------------------------------------------------

results <- list()
for (i in 50:50){
	message(paste("Working on partition:",i))
	partition <- partitions[[i]]
	# Get Modules.
	modules <- split(partition, partition)
	names(modules) <- paste0("M", names(modules))
	# Number of modules.
	nModules <- sum(names(modules) != "M0")
	message(paste("... Number of modules:",nModules))
	# Calculate Module eigengenes.
	MEdata <- moduleEigengenes(data, 
				  colors = partition, 
				  softPower = 1, 
				  impute = FALSE)
	MEs <- as.matrix(MEdata$eigengenes)
	# Calculate module membership (kME).
	KMEdata <- signedKME(data, MEs, corFnc = "bicor")
# Loop through modules, perform GSEA.
moduleGO <- list()
for (m in seq_along(modules)){
  # Calculate ranks.
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
  GSEdata <- fgsea(msPathways, ranks, minSize=5, maxSize=15, nperm = 1e+05)
  })
  # Sort by p-value.
  GSEdata <- GSEdata[order(GSEdata$pval), ]
  # Add column for gene symbols to results.
  genes <- lapply(GSEdata$leadingEdge, function(x) {
			  protmap$gene[match(x, protmap$entrez)]
			  })
  GSEdata$Symbols <- genes
  moduleGO[[m]] <- GSEdata
  } # Ends inner loop.
names(moduleGO) <- names(modules)
moduleGO <- moduleGO[seq_along(moduleGO)[names(moduleGO)!="M0"]]
  # Status report.
  nSig <- sum(sapply(moduleGO,function(x) any(x$padj<0.0)))
  message(paste("... Number of modules with any significant enrichment:",nSig))
  # Return results.
  results[[i]] <- moduleGO
} # Ends outer loop.

moduleGO <- results[[1]]
x = moduleGO[[1]]
data.table::fwrite(x,"foo.csv")
