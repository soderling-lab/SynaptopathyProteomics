#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
library(dendextend)
here <- getwd()
rootdir <- dirname(dirname(here))
datadir <- file.path(rootdir,"data")
myfile <- file.path(datadir,"wt_preserved_partitions.Rds")
data <- readRDS(myfile)

# Generate a list of all contrasts.
contrasts <- expand.grid(seq_along(data),seq_along(data))
colnames(contrasts) <- c("a1","a2")
contrasts_list <- apply(contrasts,1,function(x) list(a1=x[1],a2=x[2]))

# Loop through contrasts list, comparing partitions in data.
# This will take several minutes! 
# FIXME: Its really doing >2x the necessary work as self comparisons are 1 and matrix is symmetric.
fm <- lapply(contrasts_list, function(x) FM_index_R(data[[x$a1]],data[[x$a2]]))

# Extract the FM similarity statistic and convert this into a matrix. Labels are (P)artition(Number)
dm <- matrix(unlist(fm), nrow=length(data),ncol=length(data))
rownames(dm) <- colnames(dm) <- paste0("P",seq(dim(dm)[1]))

# Convert to distance matrix and hclust.
hc <- hclust(as.dist(1-dm), method = "ward.D2")

# Examine tree.
library(ggdendro)
ggdendrogram(hc)

# Generate groups.
k <- cutree(hc,k=4)
groups <- split(k,k)

# Get representative paritition from each group, its medoid.
# The medoid is the partition which is most similar (closest) to all others in its group.
# The distance between the medoid and all other partitions in its group should be minimized.
# Loop to get the medoid of each group:
medoid <- list()
for (i in 1:length(groups)){
	v <- names(groups[[i]])
	idx <- idy <- colnames(dm) %in% v
	subdm <- 1 - dm[idx,idy]
	diag(subdm) <- NA
	col_sums <- apply(subdm,2,function(x) sum(x,na.rm=TRUE))
	medoid[[i]] <- names(col_sums[col_sums == min(col_sums)])
}

# Examine these partitions....
unlist(medoid)















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



