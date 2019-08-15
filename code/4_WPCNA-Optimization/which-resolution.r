#!/usr/bin/env Rscript

# Which parition maximizes modularity of the ppi graph?
here <- getwd()
root <- dirname(dirname(here))
data <- paste(root,"data",sep="/")
tables <- paste(root,"tables",sep="/")

# Load compiled PPIs
file <- paste(tables,"3_Compiled_PPIs.csv",sep="/")
sif <- read.csv(file)

# Load gene gene map.
file <- paste(here,"map.csv", sep="/")
map <- read.csv(file) 
map$X <- NULL

# Add protein identifiers to SIF.
sif$protA <- map$prots[match(sif$musEntrezA, map$entrez)]
sif$protB <- map$prots[match(sif$musEntrezB, map$entrez)]

# Make igraph object.
library(igraph)
edges <- sif[,c("protA","protB")]
g <- graph_from_data_frame(d = edges, vertices = map, directed = FALSE)

# Coerce to simple graph--remove duplicate edges and self-loops.
g <- simplify(g)
is.simple(g)

# Number of nodes and edges.
length(V(g))
length(E(g))

# Load partitions.
file <- "wt_preserved_partitions.Rds"
profile <- readRDS(file)

## Calculate modularity of a given partition.

# Insure that partitions are in the same order as the vertices of the graph.
all(unlist(lapply(profile, function(x) all(names(x) == names(V(g))))))

# Calculate modularity. Add 1 to membership vector such that all are non-Zero.
q <- lapply(profile, function(x) modularity(g, membership = x + 1))

Q <- unlist(q)

#------------------------------------------------------------------------------
## Vizualize heirarchy of MEs.

# Calculate MEs for each partition.
# Use WTdata.

here <- getwd()
root <- dirname(dirname(here))
fun <- paste(root,"functions",sep="/")

source(paste(fun,"clean_fun.R",sep="/"))

# Load expression data.
wtDat <- readRDS("wtDat.Rds")
profile <- readRDS("wt_preserved_partitions.Rds")

# Calculate MEs for every partition.
# Exclude grey.
library(WGCNA)
MEdat <- lapply(profile, function(x) moduleEigengenes(wtDat, x, impute = FALSE, excludeGrey = TRUE))

# collect MEs.
MEs <- lapply(MEdat, function(x) x$eigengenes)
names(MEs) <- paste0("P",c(1:length(MEs)))

# Build dm of ME vectors.
MEdm <- do.call(cbind,MEs)

# Subset
MEdm <- MEdm[,grepl("P1\\.|P110\\.",colnames(MEdm))]

# Calculate bicor matrix.
MEcor <- silently(bicor,MEdm)

# Heirarchical clustering of the ME correlation matrix.
hc <- hclust(as.dist(1 - MEcor))

# Plot
library(ggdendro)
ggdendrogram(hc)

# Load expression data.



