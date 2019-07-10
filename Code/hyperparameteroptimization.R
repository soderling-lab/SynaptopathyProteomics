#------------------------------------------------------------------------------
## Prepare the workspace.
#------------------------------------------------------------------------------

rm(list = ls())
if (.Device != "null device" ) dev.off()
cat("\f")
options(stringsAsFactors = FALSE)

# Insure no other packages are loaded.
library(magrittr)
library(JGmisc)
detachAllPackages(keep = NULL)

# Set the working directory.
#rootdir <- "D:/projects/Synaptopathy-Proteomics"
rootdir <- "C:/Users/User/Documents/Tyler/Synaptopathy-Proteomics"
setwd(rootdir)

# Set any other directories.
Rdatadir <- paste(rootdir, "RData", sep = "/")

# Load TAMPOR cleanDat from file:
datafile <- paste(Rdatadir, "Combined","TAMPOR_data_outliersRemoved.Rds", sep = "/")
data <- log2(readRDS(datafile))

# Get a random subset of 1000 genes. Obscure labels.
rand_idx <- sample(dim(data)[1], 1000, replace = FALSE)
data <- data[rand_idx,]
rownames(data) <- paste0("protein_", c(1:nrow(data)))
colnames(data) <- paste0("sample_", c(1:ncol(data)))
saveRDS(data)

# Load dependencies:
suppressPackageStartupMessages({
  library(ParBayesianOptimization)
  library(WGCNA)
  library(igraph)
  library(doParallel)
  library(ggdendro)
})

# Allow parallel calculations:
allowWGCNAThreads()
nThreads <- 9
clusterLocal <- makeCluster(c(rep("localhost", nThreads)), type = "SOCK")
registerDoParallel(clusterLocal)

#------------------------------------------------------------------------------
## Define a function to be optimized.
#------------------------------------------------------------------------------

score_network_hclust <- function(k){
  # Calculate the signed adjaceny matrix of data.
  # Perform heirarchical clustering and divide into k groups.
  r <- 1 - bicor(t(data))
  adjm <- ((1 + r) / 2)^12 # Signed, weighted adjacency matrix. 
  diss <- 1 - adjm
  diag(diss) <- 0
  hc <- hclust(as.dist(diss), method = "average")
  membership <- cutree(hc, k)
  g <- graph_from_adjacency_matrix(adjm, mode = "undirected", weighted = TRUE)
  q <- modularity(g, membership, weights = edge_attr(g, "weight"))
  return(list(Score = q, dendro = hc))
}

# Before proceeding, test the function:
res1 <- score_network_hclust(k = 2)
res2 <- score_network_hclust(k = 20)

# Examine the results.
res1$Score
res2$Score
plot(res1$dendro, labels = FALSE)

#------------------------------------------------------------------------------
## Perform Bayesian optimization.
#------------------------------------------------------------------------------

# What is the optimal number of heirarchical groups to divide the network into?

# BayesianOptimization:
results <- BayesianOptimization(FUN = score_network_hclust, 
                                bounds = list(k = c(2L,3022L)), # maximum number of clusters is number of nodes...
                                parallel = TRUE, 
                                bulkNew = 9,
                                packages = c("WGCNA","igraph"), 
                                export = "data",
                                initPoints = 5, 
                                nIters = 100,
                                gsPoints = 100, 
                                verbose = 1)

results$BestPars
# Iteration   k      Score elapsedSecs
# 1:         0 486 0.01670862     12 secs
# 2:         1 486 0.01670862     23 secs

#------------------------------------------------------------------------------
## Are these really the best parameters?
score_network(k = 146)
score_network(k = 1000)
score_network(k = 3022)

#------------------------------------------------------------------------------

