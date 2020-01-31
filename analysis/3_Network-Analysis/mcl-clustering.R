#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

## User parameters to change:
net = "Cortex" # Which network are we analyzing? 
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
script_name <- "mcl-clustering"
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
# Get a large module.
#--------------------------------------------------------------------

max_size <- 500
moi <- all_modules[which(sapply(all_modules,length)>max_size)]
moi <- names(moi[-grep("M0",names(moi))])
module = moi[1]

#--------------------------------------------------------------------
## Check the PPI topology of large module.
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
graph <- buildNetwork(ppis, entrez, taxid = 10090)

# Remove self-connections and redundant edges.
graph <- simplify(graph)

# Set vertex attribute as protein identifiers.
ids <- protmap$ids[match(names(V(graph)),protmap$entrez)]
graph <- set_vertex_attr(graph,"name",value = ids)

# Add ppi edge attribute.
graph <- set_edge_attr(graph,"ppi",value=TRUE)

prots=all_modules[[module]]
prots <- prots[-which(prots %notin% names(V(graph)))]
subg <- induced_subgraph(graph,vids=V(graph)[match(prots,names(V(graph)))])

# Scatter plot.
connectivity <- apply(as.matrix(as_adjacency_matrix(subg)),2,sum)
p1 <- ggplotScaleFreeFit(connectivity)

#--------------------------------------------------------------------
## Create Synaptsome co-expression graph.
#--------------------------------------------------------------------

# Create co-expression graph.
graph <- graph_from_adjacency_matrix(adjm,mode="undirected",weighted=TRUE)
graph <- simplify(graph)

# Subset
prots=all_modules[[module]]
subg <- induced_subgraph(graph,vids=V(graph)[match(prots,names(V(graph)))])

# Get best soft-threshold.
idy <- colnames(data) %in% prots
subdat <- data[,idy]
sft <- pickSoftThreshold(subdat,corFnc="bicor",
			 networkType="signed",
			 RsquaredCut=0.8)
sft_power <- sft$powerEstimate
stats <- subset(sft$fitIndices,sft$fitIndices$Power==sft$powerEstimate)

# Threshold the co-expresion graph.
subg <- set_edge_attr(subg,'weight',value = edge_attr(subg,'weight')^sft_power)

# Run MCL.
partitions <- list()
inflation <- seq(1.2,5,by=0.1)
Q <- vector("numeric",length(inflation))
for (i in seq_along(inflation)){
	partition <- clusterMCL(subg,weight="weight",inflation[i])
	q <- modularity(subg,membership=partition[names(V(subg))],
			weights=abs(edge_attr(subg,'weight')))
	# Return modularity and partition.
	Q[i] <- q
	partitions[[i]] <- partition
	# Status.
	message(paste0("Inflation: ", inflation[i],
		       "; Modularity: ",round(q,3)))
}

