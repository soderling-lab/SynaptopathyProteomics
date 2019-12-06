#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
  require(WGCNA)
  require(igraph)
  require(getPPIs)
  require(data.table)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatdir <- file.path(root, "rdata")
funcdir <- file.path(root, "R")

# Load functions.
myfun <- list.files(funcdir, pattern = "silently.R", full.names = TRUE)
invisible(sapply(myfun, source))

# Load statistical results.
myfile <- list.files(rdatdir, pattern = "GLM", full.names = TRUE)
data <- readRDS(myfile)

# Load mouse ppi interaction data.
data("musInteractome")

# Load expression data. Transpose -> rows = samples; columns = genes.
wtDat <- t(readRDS(file.path(rdatdir, "3_WT_cleanDat.RData")))
koDat <- t(readRDS(file.path(rdatdir, "3_KO_cleanDat.RData")))

# Compute adjmatrix:
wtAdjm <- silently(WGCNA::bicor(wtDat))
koAdjm <- silently(WGCNA::bicor(koDat))

# Load protein id map.
myfile <- list.files(rdatdir, "Map", full.names = TRUE)
protmap <- readRDS(myfile)

# Given expression data, calculate power for ~scale free fit:
sft <- silently({
  sapply(list(wtDat, koDat), function(x) {
    pickSoftThreshold(x,
      corFnc = "bicor",
      networkType = "signed",
      RsquaredCut = 0.8
    )$powerEstimate
  })
})
names(sft) <- c("wt", "ko")

# Map all proteins to entrez.
prots <- colnames(wtAdjm)
entrez <- protmap$entrez[match(prots, protmap$ids)]

# Keep ppi data from mouse,human, and rat.
idx <- musInteractome$Interactor_A_Taxonomy %in% c(10090, 9606, 10116)
ppis <- subset(musInteractome, idx)

# Build a graph with all proteins.
graph <- buildNetwork(ppis, entrez, taxid = 10090)

# Evaluate scale free fit.
dc <- apply(as.matrix(as_adjacency_matrix(g)), 2, sum) # calc degree connectivity.
fit <- scaleFreeFitIndex(dc, nBreaks = 10, removeFirst = FALSE)

#------------------------------------------------------------------------------

# Melt stats, tidy data easier to work with.
df <- melt(data, id = "Uniprot")

# Collect in a list.
sig_prots <- df %>%
  filter(value < 0.05) %>%
  group_by(variable) %>%
  group_split()
names(sig_prots) <- levels(df$variable)

# Get Entrez for all sig genes.
sig_entrez <- lapply(sig_prots, function(x) protmap$entrez[match(x$Uniprot, protmap$uniprot)])

# Number of seeds:
n_seeds <- sapply(sig_entrez, length)

# Build community graphs.
graphs <- lapply(sig_entrez, function(x) getCommunity(graph, x, k = 2))

# Number of nodes per community.
n_nodes <- sapply(graphs, function(x) length(V(x)))

# Is there overlap?
