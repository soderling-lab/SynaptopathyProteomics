#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

## User parameters to change:
net = "Striatum" # Which network are we analyzing? 
overwrite_figsdir = TRUE
do_DBD_enrichment = FALSE
do_GO_enrichment = FALSE
do_module_analysis = FALSE

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
  library(igraph)
  library(vegan)
  library(ggplot2)
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
netsdir <- file.path(root, "networks")

# Create directory for figure output.
script_name <- "4_network-analysis"
figsdir <- file.path(root,"figs",Sys.Date(),script_name)
if(overwrite_figsdir) {
  message(paste("Warning, overwriting files in:\n",figsdir))
  unlink(figsdir,recursive = TRUE)
  dir.create(figsdir,recursive = TRUE)
} else if (!dir.exists(figsdir)) {
  dir.create(figsdir,recursive = TRUE)
}

# Functions.
source_funcdir <- function(funcdir){
  myfun <- list.files(funcdir, pattern="*.R",full.names = TRUE)
  invisible(sapply(myfun, source))
}
source_funcdir(funcdir)

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load GLM stats.
myfile <- file.path(rdatdir, "2_GLM_Stats.RData")
glm_stats <- readRDS(myfile)

# Proteins with any significant change.
sigProts <- apply(glm_stats$FDR,1,function(x) any(x<0.05))
sigProts <- names(sigProts)[sigProts]

# SigProts by genotype.
fdr_df <- glm_stats$FDR
fdr_df$Protein <- rownames(fdr_df)
sig_df <- melt(fdr_df,id="Protein") %>% 
	group_by(variable) %>% filter(value < 0.05) %>% group_split()
sigProts_geno <- sapply(sig_df,function(x) x$Protein)
names(sigProts_geno) <- gsub(" ","_", gsub(" FDR", "", colnames(fdr_df)[1:8]))

# Load expression data.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_cleanDat.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_cleanDat.RData")
)
data <- t(readRDS(myfiles[net])) # Data should be transposed: rows, proteins.

# Load Sample info.
traits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load co-expression (adjacency) matrices.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_Adjm.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_Adjm.RData")
)
adjm <- as.matrix(readRDS(myfiles[net]))
rownames(adjm) <- colnames(adjm)

# Load GO semantic similarity graph.
myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_RMS_Adjm.csv")
adjm_go <- fread(myfile,drop=1)

# Load network partitions-- self-preservation enforced.
ids <- c("Cortex"="14942508","Striatum"="14940918")
myfile <- list.files(rdatdir, pattern = ids[net], full.names = TRUE)
partitions <- readRDS(myfile) 

# Reset partition index.
partitions <- lapply(partitions, reset_index)

# Load theme for plots.
ggtheme()

#------------------------------------------------------------------------------
## Collect all modules in a named list.
#------------------------------------------------------------------------------

# Name partitions.
named_parts <- partitions
names(named_parts) <- paste0("R",c(1:length(partitions)))

# Name modules.
modules_list <- lapply(named_parts,function(x) split(x,x))
modules_list <- sapply(modules_list,function(x) {
  names(x) <- paste0("M",names(x))
  return(x) })

# All named modules.
all_modules <- unlist(modules_list,recursive=FALSE)
all_modules <- sapply(all_modules,names)

#--------------------------------------------------------------------
## Create Synaptosome PPI graph.
#--------------------------------------------------------------------

# Load all ppis mapped to mouse genes.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
idx <- musInteractome$Interactor_A_Taxonomy %in% c(10090, 9606, 10116)
ppis <- subset(musInteractome, idx)

# Get entrez IDs for all proteins in data.
prots <- colnames(data)
entrez <- protmap$entrez[match(prots, protmap$ids)]

# Build a ppi graph with all proteins.
g1 <- buildNetwork(ppis, entrez, taxid = 10090)

# Remove self-connections and redundant edges.
g1 <- simplify(g1)

# Set vertex attribute as protein identifiers.
ids <- protmap$ids[match(names(V(g1)),protmap$entrez)]
g1 <- set_vertex_attr(g1,"name",value = ids)

# Add ppi edge attribute.
g1 <- set_edge_attr(g1,"ppi",value=TRUE)

#--------------------------------------------------------------------
## Plot topology of ppi graph.
#--------------------------------------------------------------------

# Check topology of PPI graph.
adjm_ppi <- as_adjacency_matrix(g1)

# node degree is column sum.
connectivity <- apply(adjm_ppi,2,sum) 

# Scatter plot.
p1 <- ggplotScaleFreeFit(connectivity)
myfile <- prefix_file({
  file.path(figsdir,paste0(net,"_ScaleFreeFit.tiff"))
  })
ggsave(myfile,p1,width=3,height=3)

# Histogram.
p2 <- ggplotHistConnectivity(connectivity)
myfile <- prefix_file({
  file.path(figsdir,paste0(net,"_ScaleFreeHist.tiff"))
  })
ggsave(myfile,p2,width=3,height=3)



# Break down a large module with MCL.
module = rep_dbd_modules[1]
prots=all_modules[[module]]

subg <- induced_subgraph(g0,vids=V(g0)[match(prots,names(V(g0)))])

# Scatter plot.
connectivity <- apply(adjm_ppi,2,sum) 
p1 <- ggplotScaleFreeFit(connectivity)

# Threshold the graph.
subg <- set_edge_attr(subg,'weight',value = edge_attr(subg,'weight')^11)

# Run MCL.
partition <- clusterMCL(subg,weight="weight",inflation=2.5)

# By thresholding the graph we are able to resolve its structure.
table(partition)


