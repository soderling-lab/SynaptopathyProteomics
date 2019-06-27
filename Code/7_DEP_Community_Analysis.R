#' ---
#' title: Building DEP communities.
#' author: Tyler W Bradshaw
#' urlcolor: blue
#' header-includes:
#' - \usepackage{float}
#' - \floatplacement{figure}{H}
#' output:
#'    pdf_document:
#'      fig_caption: true
#'      toc: true
#'      number_sections: false
#'      highlight: tango
#' ---

#-------------------------------------------------------------------------------
#' ## Prepare the workspace.
#-------------------------------------------------------------------------------
#+ eval = TRUE, echo = FALSE, error = FALSE

# Use ctl+alt+T to execute a code chunk.
# Use ctl+shift+W to close all tabs.

# Run this chunk before doing anything!
rm(list = ls())
dev.off()
cat("\014") # alternative is cat("\f")
options(stringsAsFactors = FALSE)

# If you have not cleared the workspace of all loaded packages, you may
# incounter problems. To remove all packages, you can call the following:
library(magrittr)
library(JGmisc)
detachAllPackages(keep = NULL)

#  Load required packages.
suppressPackageStartupMessages({
  library(JGmisc)
  library(readxl)
  library(knitr)
  library(readr)
  library(dplyr)
  library(reshape2)
  library(DEP)
  library(tibble)
  library(SummarizedExperiment)
  library(ggplot2)
  library(hexbin)
  library(vsn)
  library(BurStMisc)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(edgeR)
  library(openxlsx)
  library(stringr)
  library(imp4p)
  library(Cairo)
  library(pryr)
  library(qvalue)
  library(gridExtra)
  library(cowplot)
  library(WGCNA)
  library(impute)
  library(ggrepel)
  library(sva)
  library(anRichment)
  library(ggdendro)
  library(flashClust)
  library(purrr)
  library(ggpubr)
  library(doParallel)
  library(NMF)
  library(FSA)
  library(plyr)
  library(RColorBrewer)
  library(gtable)
  library(grid)
  library(ggplotify)
  library(igraph)
  library(RCy3)
  library(DescTools)
  library(TBmiscr)
  library(GOSemSim)
})

# Define version of the code.
CodeVersion <- "DEP_Community_Analysis"

# Define tisue type: cortex = 1; striatum = 2; combined = 3.
type <- 3
tissue <- c("Cortex", "Striatum", "Combined")[type]

# Set the working directory.
rootdir <- "D:/Documents/R/Synaptopathy-Proteomics"
setwd(rootdir)

# Set any other directories.
functiondir <- paste(rootdir, "Code", sep = "/")
datadir <- paste(rootdir, "Input", sep = "/")
Rdatadir <- paste(rootdir, "RData", sep = "/")

# Create code-version specific directories for figures and tables.
# Creat otuput direcotry for figures.
outputfigs <- paste(rootdir, "Figures", tissue, sep = "/")
outputfigsdir <- paste(outputfigs, CodeVersion, sep = "/")
if (!file.exists(outputfigsdir)) {
  dir.create(file.path(outputfigsdir))
} else {
  print("This directory already exists. Warning: Some files may be overwritten when running this script.")
}
# Create output directory for tables.
outputtabs <- paste(rootdir, "Tables", tissue, sep = "/")
outputtabsdir <- paste(outputtabs, CodeVersion, sep = "/")
if (!file.exists(outputtabsdir)) {
  dir.create(file.path(outputtabsdir))
} else {
  print("This directory already exists. Warning: Some files may be overwritten when running this script.")
}
# Create output directory for reports.
outputreports <- paste(rootdir, "Reports", tissue, sep = "/")
outputrepsdir <- paste(outputreports, CodeVersion, sep = "/")
if (!file.exists(outputrepsdir)) {
  dir.create(file.path(outputrepsdir))
} else {
  print("This directory already exists. Warning: Some files may be overwritten when running this script.")
}

# Load required custom functions.
my_functions <- paste(functiondir, "0_TMT_Preprocess_Functions.R", sep = "/")
source(my_functions)

# Define prefix for output figures and tables.
outputMatName <- paste(tissue, "_Network_Analysis_", sep = "")

# Globally set ggplots theme.
ggplot2::theme_set(theme_gray())

#-------------------------------------------------------------------------------
#' ## Load the WGNCA and TMT data.
#-------------------------------------------------------------------------------

# Load TAMPOR cleanDat from file:
datafile <- paste(Rdatadir, tissue, "TAMPOR_data_outliersRemoved.Rds", sep = "/")
cleanDat <- readRDS(datafile)
cleanDat <- log2(cleanDat)
cleanDat[1:5, 1:5]
dim(cleanDat)

# Load WGCNA network and meta data.
file <- paste(Rdatadir, "Network_and_metaModules.Rds", sep = "/")
data <- readRDS(file)
net <- data$net
meta <- data$meta

# Map 3 un-mapped entrez IDs by hand.
meta[is.na(meta$entrez), ]
not_mapped <- list(
  "Ndufb1|P0DN34" = 102631912,
  "F8a1|Q00558" = 14070,
  "Pc|Q05920" = 18563
)
for (i in 1:length(not_mapped)) {
  idx <- meta$protein == names(not_mapped)[i]
  meta$entrez[idx] <- not_mapped[[i]]
}
# Check.
sum(is.na(meta$entrez))

# Add gene names and uniprot with entrez for consistent mapping.
meta$gene <- mapIds(org.Mm.eg.db,
  keys = meta$entrez, column = "SYMBOL",
  keytype = "ENTREZID", multiVals = "first"
)

# Add uniprot.
meta$uniprot <- sapply(strsplit(meta$protein, "\\|"), "[", 2)

# Load TAMPOR statistical results.
file <- paste(outputtabs, "Final_TAMPOR",
  "Combined_TMT_Analysis_TAMPOR_GLM_Results.xlsx",
  sep = "/"
)
results <- lapply(as.list(c(1:8)), function(x) read_excel(file, x))
names(results) <- excel_sheets(file)

# Load DEP PPI communities.
file <- paste0(Rdatadir, "/", "DEP_Communities.Rds")
DEP_communities <- readRDS(file)

#-------------------------------------------------------------------------------
#' ## Build df of significantly dysregulated proteins.
#-------------------------------------------------------------------------------

# Create a dictionary like object for mapping entrez to gene|uniprot.
entrez2protein <- as.list(meta$protein)
names(entrez2protein) <- meta$entrez

protein2entrez <- as.list(meta$entrez)
names(protein2entrez) <- meta$protein

uniprot2entrez <- as.list(meta$entrez)
names(uniprot2entrez) <- meta$uniprot

gene2entrez <- as.list(meta$entrez)
names(gene2entrez) <- meta$gene

uniprot2gene <- as.list(meta$gene)
names(uniprot2gene) <- meta$uniprot

# Build a df with statistical results.
stats <- lapply(results, function(x)
  data.frame(
    Uniprot = x$Uniprot,
    Entrez = x$Entrez,
    FDR = x$FDR
  ))
names(stats) <- names(results)
df <- stats %>% reduce(left_join, by = c("Uniprot", "Entrez"))
colnames(df)[c(3:ncol(df))] <- names(stats)

# Proteins with any significant change.
df$sigProt <- apply(df, 1, function(x) any(as.numeric(x[c(3:ncol(df))]) < 0.05))
stats_df <- df

# Insure that all Uniprot have been mapped to entrez.
stats_df$Entrez <- unlist(uniprot2entrez[stats_df$Uniprot])

# Check.
sum(is.na(stats_df$Entrez)) == 0

#-------------------------------------------------------------------------------
#' ## Combine sigprots from cortex and striatum.
#-------------------------------------------------------------------------------

# Combine tissues for each genotype.
genos <- c("Shank2","Shank3","Syngap1","Ube3a")
u <- list()
for (geno in genos){
  idx <- grep(geno,names(sigProts))
  u[[geno]] <- union(sigProts[[idx[1]]],sigProts[[idx[2]]])
}
sigProts <- u

# Map sigProts to Entrez.
sigEntrez <- lapply(sigProts,function(x) unlist(uniprot2entrez[x]))

#-------------------------------------------------------------------------------
#' Load the PPI graph.
#-------------------------------------------------------------------------------

# Load the data.
dir <- "D:/Documents/R/Synaptopathy-Proteomics/Tables/Network"
file <- paste(dir,"SIF.xlsx",sep="/")
sif <- read_excel(file, sheet = 1)

# Send network to cytoscape?
send_to_cytoscape = FALSE 

# Create a data frame with all node attributes. 
nodes <- data.frame(Entrez = unlist(meta$entrez[!is.na(meta$entrez)]),
                    Uniprot = unlist(meta$uniprot[!is.na(meta$entrez)]),
                    Symbol = unlist(meta$gene[!is.na(meta$entrez)]),
                    Module = unlist(meta$module[!is.na(meta$entrez)]),
                    MetaModule = unlist(meta$metaModule[!is.na(meta$entrez)]))
nodes$MetaModule <- paste0("MM",nodes$MetaModule)
colors <- unlist(lapply(as.list(meta$module),function(x) col2hex(x)))
colors <- colors[!is.na(meta$entrez)]
nodes$ModulColor <- colors

# Compile EdgeR GLM stats. Build a df with logfc and p-values.
stats <- lapply(results, function(x)
  data.frame(Uniprot = x$Uniprot, 
             Entrez = x$Entrez, 
             logFC = x$logFC,
             PValue = x$PValue,
             FDR = x$FDR))
names(stats) <- names(results)
df <- stats %>% reduce(left_join, by = c("Uniprot","Entrez"))
namen <- sapply(strsplit(colnames(df)[c(3:ncol(df))],"\\."),"[",1)
colnames(df)[c(3:ncol(df))] <- paste(namen,rep(names(stats),each=3))

# Combine with node attributes. 
# use uniprot to combine as this is the MOST stable identifier for all proteins.
nodes <- merge(nodes, df, by = c("Uniprot"))

# Clean-up the df.
nodes$Entrez.y <- NULL
names(nodes)[2] <- "Entrez"
nodes <- nodes[,c(2,1,3:ncol(nodes))] # reorder so that entrez is first col
dim(nodes)

# Add labels for Shank2, Shank3, Syngap1, and Ube3a sigProts.
groups <- paste(c("Shank2","Shank3","Syngap1","Ube3a"), collapse="|")
cols <- grepl(groups,colnames(nodes)) & grepl("FDR",colnames(nodes))
dat <- nodes[,cols]
genos <- c("Shank2","Shank3","Syngap1","Ube3a")
out <- list()
for (i in seq_along(genos)){
  geno <- genos[i]
  idy <- grep(geno,colnames(dat))
  sub <- dat[,idy]
  out[[geno]] <- apply(sub,1,function(x) any(x<0.05))
}
dm <- do.call(cbind,out)
colnames(dm) <- paste(colnames(dm),"sigProt", sep="_")

# Any significant change.
nodes$sigProt <- apply(dat, 1, function(x) any(x<0.05))

# Determine node color by mixing genotype colors. 
# Custom colors:
colvec <- c("#FFF200", "#00A2E8", "#22B14C", "#A349A4")

node_colors <- list()
for (i in 1:nrow(dm)){
  x <- dm[i,]
  if (sum(x)>1) {
    node_colors[[i]] <- do.call(mixcolors, as.list(colvec[x]))
  } else if (sum(x)==1) { 
    node_colors[[i]] <- col2hex(colvec[x])
  } else if (sum(x)==0) {
    node_colors[[i]] <- col2hex("gray")
  }
}
nodes$NodeColor <- as.character(do.call(rbind, node_colors))

# Rownames need to be same and node names for loading into cytoscape. 
rownames(nodes) <- nodes$Entrez

# Make igraph object. 
g <- graph_from_data_frame(d=sif, vertices=nodes, directed=FALSE)

# Coerce to simple graph--remove duplicate edges and self-loops.
g <- simplify(g)
is.simple(g)

# Number of nodes and edges. 
length(V(g)) # All but three unmapped genes. 
length(E(g))

#-------------------------------------------------------------------------------
#' ## Find KNN in co-expression space.
#-------------------------------------------------------------------------------

# Calculate protein co-expression (correlation) matrix.
# Row and column names are entrez ids. 
r <- bicor(t(cleanDat))
diag(r) <- NA
colnames(r) <- unlist(protein2entrez[colnames(r)])
rownames(r) <- unlist(protein2entrez[rownames(r)])

# Define a function to get KNN neighbors.
get_knn <- function(x, k) {
  x.sorted <- sort(x, decreasing = TRUE)
  x.knn <- x.sorted[c(1:k)]
  out <- list(names(x.knn))
  names(out) <- names(x)[1]
  return(out)
}

# Get KNN for all nodes.
knn_list <- apply(r, 1, function(x) get_knn(x, k=3))

#-------------------------------------------------------------------------------
#' ## Build subraphs for DEPs communities. 
#-------------------------------------------------------------------------------

# Generate PPIs graphs using DEPs from each tissue:genotype contrast as seed 
# nodes. Add nodes with 2+ connections to these seed nodes for biological 
# context. Do not consider connections to 1433* chaperone proteins. Also add
# seed nearest neighbors in co-expression space (k=3).

# Defaults for analysis. 
degree_to_stay    = 2     # 2 degrees to seed nodes.

# Get DEP seeds neighborhoods. 
subg <- lapply(sigEntrez, function(x) make_ego_graph(g, nodes = x, mode = "all"))

# Combine subgraphs, union. 
uniong <- lapply(subg, function(x) do.call(igraph::union,x))

# Calculate distances to all nodes. 
dist <- lapply(uniong, function(x) distances(x, mode = "all", algorithm = "unweighted"))

# Keep distances to seed nodes. 
f <- function(x,y) { x <- x[,colnames(x) %in% y]; return(x) }
dist <- mapply(f,dist,sigEntrez)

# Only consider direct connections (distance == 1).
f <- function(x) { x[x!=1] <- 0; return(x) }
dist <- lapply(dist, function(x) f(x))

# Exclude Ywha* genes (14-3-3 proteins).
out <- as.character(meta$entrez[grep("Ywha*", meta$gene)])
f <- function(x) { x[rownames(x) %in% out,] <- 0; return(x)}
dist <- lapply(dist, function(x) f(x))

# Calculate degree to seed nodes (row sum).
deg <- lapply(dist, function(x) apply(x,1,function(y) sum(y)))

# We will keep nodes that have at least 2 connections with seed nodes.
keep <- lapply(deg, function(x) names(x)[x >= degree_to_stay])

# Combine with seed nodes.
community_nodes <- mapply(c,keep,sigEntrez)

# Add KNN.
knn_nodes <- lapply(sigEntrez, function(x) as.character(unlist(knn_list[x])))

# Combine.
combined_nodes <- lapply(mapply(c,community_nodes,knn_nodes), function(x) unique(x))

# Return output. 
community_results <- list("seed_nodes" = sigEntrez,
                          "community_nodes" = community_nodes,
                          "knn_nodes" = knn_nodes, 
                          "combined_nodes" = combined_nodes)

# Reorganize output by genotype. 
out <- list()
for (i in 1:4){
  out[[i]] <- sapply(community_results,"[",i)
}
names(out) <- c("Shank2","Shank3","Syngap1","Ube3a")
community_results <- out

# Summarize number of nodes.
unlist(lapply(community_results$Shank2, function(x) length(x)))
unlist(lapply(community_results$Shank3, function(x) length(x)))
unlist(lapply(community_results$Syngap1, function(x) length(x)))
unlist(lapply(community_results$Ube3a, function(x) length(x)))

#-------------------------------------------------------------------------------
#' ## Send DEP communities to Cytoscape!
#-------------------------------------------------------------------------------

send_to_cytoscape = FALSE

################################################################################
# Cytoscape should be open before proceeding!
################################################################################

# Custom colors:
colors <- as.list(c("#FFF200", "#00A2E8", "#22B14C", "#A349A4"))
names(colors) <- c("Shank2","Shank3","Syngap1","Ube3a")

for (i in 1:length(sigEntrez)){
  
  geno <- names(sigEntrez)[i]
  print(paste0("Working on ", geno, " subgraph", " (i=",i,")", "..."))
  
  # Get subset of nodes (v) in modules overlap.
  # We will use these to seed a network.
  v <- sigEntrez[[i]] # number of significant proteins for this exp.
  
  # Insure that all nodes are in the network.
  #table(v %in% vertex_attr(g, "name"))
  v <- v[v %in% vertex_attr(g, "name")]
  print(paste0("Seed Nodes: ", length(v)))
  seeds <- v
  
  # Create list of subgraphs (subg) for every seed node.
  subg <- make_ego_graph(g, 
                         order = 1, 
                         nodes = v,
                         mode = "all", 
                         mindist = 0)
  
  # Combine subgraphs, union. 
  uniong <- do.call(igraph::union,subg)
  
  # Calculate distances from seed nodes to all else.
  dist <- as.data.frame(
    distances(uniong,
              v = V(uniong), 
              to = v, 
              mode = "all",
              weights = NULL, 
              algorithm = "unweighted")
  )
  
  # Only consider direct connections to seed nodes (distance == 1).
  dist[dist!=1] <- 0
  
  # Exclude Ywha* genes (14-3-3 proteins).
  out <- as.character(unlist(gene2entrez[grep("Ywha*",names(gene2entrez))]))
  dist[rownames(dist) %in% out,] <- 0
  
  # Calculate degree to seed nodes (sum).
  dist$SeedDegree <- apply(dist[,-ncol(dist)],1,function(x) sum(x))
  
  # We will keep nodes that have at least 2 connections with seed nodes.
  # degree_to_stay = 2
  keep <- dist$SeedDegree>=degree_to_stay
  dist <- dist[keep,]
  keepers <- unique(c(v,rownames(dist)))
  print(paste("Additional nodes:", length(keepers) - length(seeds)))
  subg <- induced_subgraph(g, keepers)
  
  # Build df of node attributes. 
  df <- data.frame(sigProt = names(V(subg)) %in% sigEntrez[[geno]])
  rownames(df) <- names(V(subg))
  df$sigProt <- as.factor(df$sigProt)
  
  # How many nodes. 
  print(paste0("Total Nodes: ", length(V(subg))))
  
  # result
  community_results[[i]] <- list(seeds= seeds, 
                                 subg = subg, 
                                 nodes = names(V(subg)),
                                 nodes.entrez = keepers)
  
  # Send to cytoscape with RCy3!
  if (send_to_cytoscape == TRUE){
    cytoscapePing()
    quiet(RCy3::createNetworkFromIgraph(subg,names(sigEntrez)[i]))
    
    # Load node attribute table in cytoscape.
    loadTableData(df)
    
    # Create custom syle to customize appearance. 
    #geno <- strsplit(names(sigEntrez)[i], "\\.")[[1]][3]
    #colvec <- as.character(c(col2hex("gray"), unlist(colors[geno])))
    defaults <- list(NODE_SHAPE        = "Ellipse",
                     NODE_SIZE         = 55,
                     EDGE_WIDTH        = 2.0,
                     EDGE_TRANSPARENCY = 120)
    nodeLabels <- mapVisualProperty('node label','Symbol','p')
    #nodeFills <- mapVisualProperty('node fill color','sigProt','d',c(FALSE,TRUE), colvec) # why does this not work???
    #setNodeColorMapping("NodeColor", mapping.type = 'p')
    setNodeColorBypass(sigEntrez[[geno]], colors[[geno]], network = getNetworkSuid(geno))
    setNodeSizeBypass(sigEntrez[[geno]], new.sizes = 75, network = getNetworkSuid(geno))
    #edgeWidth <- mapVisualProperty('edge width','weight','p')
    createVisualStyle(style.name = geno, defaults, list(nodeLabels))
    lockNodeDimensions(TRUE, style.name = geno)
    setVisualStyle(style.name = geno)
    # Apply perfuse force directed layout. 
    layoutNetwork(layout.name = "force-directed")
    setNodeColorDefault(col2hex("grey"), style.name = geno)
  }
}

#-------------------------------------------------------------------------------
#' ## Generate randomly seeded subgraphs...
#-------------------------------------------------------------------------------

n_iter <- 1000
n_seeds <- lapply(DEP_communities, function(x) length(x$seeds))
output <- list()

# Loop to generate randomly seeded graphs. 
for (i in 1:n_iter){
  print(paste("Generating randomly seeded communities, iteration",i,"..."))
  # Randomly Sample nodes.
  rand_seeds <- lapply(n_seeds, function(x) sample(meta$entrez,x, replace = FALSE))
  # Get random seeds' neighborhoods. 
  subg <- lapply(rand_seeds, function(x) make_ego_graph(g, nodes = x, mode = "all"))
  # Combine subgraphs, union. 
  uniong <- lapply(subg, function(x) do.call(igraph::union,x))
  # Calculate distances to all nodes. 
  dist <- lapply(uniong, function(x) distances(x, mode = "all", algorithm = "unweighted"))
  # Keep distances to seed nodes. 
  f <- function(x,y) { x <- x[,colnames(x) %in% y]; return(x) }
  dist <- mapply(f,dist,rand_seeds)
  # Only consider direct connections (distance == 1).
  f <- function(x) { x[x!=1] <- 0; return(x) }
  dist <- lapply(dist, function(x) f(x))
  # Exclude Ywha* genes (14-3-3 proteins).
  out <- as.character(meta$entrez[grep("Ywha*", meta$gene)])
  f <- function(x) { x[rownames(x) %in% out,] <- 0; return(x)}
  dist <- lapply(dist, function(x) f(x))
  # Calculate degree to seed nodes (row sum).
  deg <- lapply(dist, function(x) apply(x,1,function(y) sum(y)))
  # We will keep nodes that have at least 2 connections with seed nodes.
  keep <- lapply(deg, function(x) names(x)[x >= 2])
  # Combine with seed nodes.
  community_nodes <- mapply(c,keep,rand_seeds)
  # Add KNN.
  knn_nodes <- lapply(rand_seeds, function(x) as.character(unlist(knn_list[x])))
  # Combine.
  combined_nodes <- lapply(mapply(c,community_nodes,knn_nodes), function(x) unique(x))
  # Return output. 
  out <- list("seed_nodes" = rand_seeds,
              "community_nodes" = community_nodes,
              "knn_nodes" = knn_nodes, 
              "combined_nodes" = combined_nodes)
  output[[i]] <- out
  names(output) <- paste("iter",i,sep="_")
  # Save result.
  if (i == c(1,seq(0,n_iter,50))){
    print("Saving progress to .RDS!")
    
    file <- paste(Rdatadir,"/","Random_Communities.RDS")
    saveRDS(output, file)
  }
  
}



###################################################
# Check...
unlist(lapply(output$iter_1$seed_nodes,function(x) length(x)))
unlist(lapply(output$iter_1$community_nodes,function(x) length(x)))
unlist(lapply(output$iter_1$knn_nodes,function(x) length(x)))
unlist(lapply(output$iter_1$combined_nodes,function(x) length(x)))
###################################################

#-------------------------------------------------------------------------------
#' ## Examine protein overlap between communities.
#-------------------------------------------------------------------------------

# Fixme: change color of plot. 

# Examine overlap in ppi networks...
nodes <- sapply(DEP_communities, "[", 3)
names(nodes) <- sapply(strsplit(names(nodes), "\\."), "[", 1)

# Build a matrix showing overlap.
col_names <- names(nodes)
row_names <- names(nodes)

# All possible combinations.
contrasts <- expand.grid(col_names, row_names)
colnames(contrasts) <- c("ConditionA", "ConditionB")
contrasts$ConditionA <- as.vector(contrasts$ConditionA)
contrasts$ConditionB <- as.vector(contrasts$ConditionB)

# Loop to calculate intersection for all contrasts.
int <- list()
A <- list()
B <- list()
for (i in 1:dim(contrasts)[1]) {
  a <- unlist(nodes[contrasts$ConditionA[i]])
  b <- unlist(nodes[contrasts$ConditionB[i]])
  A[[i]] <- length(a)
  B[[i]] <- length(b)
  int[[i]] <- intersect(a, b)
}

contrasts$Name <- paste(contrasts$ConditionA,
  contrasts$ConditionB,
  sep = "_U_"
)
names(int) <- contrasts$Name

# Add intersection to contrasts.
contrasts$Intersection <- lapply(int, function(x) length(x))
contrasts$A <- unlist(A)
contrasts$B <- unlist(B)
contrasts$C <- contrasts$A + contrasts$B

# Calculate percent intersection.
int <- list()
for (i in 1:dim(contrasts)[1]) {
  a <- unlist(nodes[contrasts$ConditionA[i]])
  b <- unlist(nodes[contrasts$ConditionB[i]])
  int[[i]] <- length(intersect(a, b)) / length(unique(c(a, b)))
}

# Add to contrasts dm
contrasts$Percent <- unlist(int)

# Make overlap matrix.
n <- length(nodes)
dm <- matrix(unlist(contrasts$Intersection), nrow = n, ncol = n)
rownames(dm) <- colnames(dm) <- row_names

# heirarchical clustering of dissimilarity matrix calculated as 1-percent_overlap.
dm <- matrix(contrasts$Percent, nrow = n, ncol = n)
rownames(dm) <- colnames(dm) <- row_names
diss <- 1 - dm
hc <- hclust(as.dist(diss), method = "average")
dendro <- ggdendrogram(hc, rotate = TRUE, labels = FALSE)
# Strip labels.
dendro <- dendro + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
dendro

# Remove upper tri and melt.
dm[lower.tri(dm)] <- NA
df <- melt(dm, na.rm = TRUE)

# Add intersection.
idx <- match(paste(df$Var1, df$Var2), paste(contrasts$ConditionA, contrasts$ConditionB))
df$intersection <- unlist(contrasts$Intersection[idx])

# Generate plot.
# Order df based on dendrogram.
levels(df$Var1) <- hc$labels[hc$order]
levels(df$Var2) <- hc$labels[hc$order]

plot <- ggplot(df, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "black") +
  geom_text(aes(Var2, Var1, label = intersection), color = "black", size = 4) +
  scale_fill_gradient2(name = "Percent Overlap") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.7),
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    barwidth = 7, barheight = 1,
    title.position = "top", title.hjust = 0.5
  )) +
  coord_fixed()

#plot <- plot + theme(legend.position = "none")
plot

# Save heatmap and dendrogram.
file <- paste0(outputfigsdir, "/", outputMatName, "DEP_Community_Overlap_matirx.eps")
ggsave(file, plot, width = 3, height = 3, units = "in", dpi = 300)

file <- paste0(outputfigsdir,"/",outputMatName,"DEP_Community_Overlap_dendro.eps")
ggsave(file,dendro, width = 3, height = 3, units = "in", dpi = 300)


#-------------------------------------------------------------------------------
#' Evaluate topology of the protein subgraphs.
#-------------------------------------------------------------------------------

# Do DEP communities exhibit a scale free topology?

# Load the data.
dir <- "D:/Documents/R/Synaptopathy-Proteomics/Tables/Network"
file <- paste(dir, "SIF.xlsx", sep = "/")
sif <- read_excel(file, sheet = 1)

# Create a data frame with all node attributes.
nodes <- data.frame(
  Entrez = unlist(meta$entrez[!is.na(meta$entrez)]),
  Uniprot = unlist(meta$uniprot[!is.na(meta$entrez)]),
  Symbol = unlist(meta$gene[!is.na(meta$entrez)])
)
# Make igraph object.
g <- graph_from_data_frame(d = sif, vertices = nodes, directed = FALSE)

# Coerce to simple graph--remove duplicate edges and self-loops.
g <- simplify(g)
is.simple(g)

# Number of nodes and edges.
length(V(g)) # All but three unmapped genes.
length(E(g))

# Number of connected components.
connected_components <- components(g)
connected_components$csize[1] # The largest connected component.

# Subset the graph.
prots <- sapply(DEP_communities, "[", 3)
subgraphs <- lapply(prots, function(x) induced_subgraph(g, x))

# Calculate node connectivity (degree).
connectivity <- lapply(subgraphs, function(x) degree(x, loops = FALSE))

# Evaluate scale free fit with WGCNA function scaleFreePlot()
plot_data <- lapply(connectivity, function(x)
  ggplotScaleFreePlot(x, nBreaks = 10))
plots <- sapply(plot_data, "[", 1)
names(plots) <- sapply(strsplit(names(plots), "\\."), "[", 1)

# Check the fit.
# Shank2 does not exhibit scale free fit!
plot_grid(plots$Shank2, plots$Shank3, plots$Syngap1, plots$Ube3a)

# Save figures.
for (i in 1:length(plots)) {
  namen <- names(plots)[i]
  file <- paste0(outputfigsdir, "/", namen, "_Community_ScaleFreeFit.eps")
  ggsave(file, plot, width = 3, height = 2.5, units = "in")
}

#-------------------------------------------------------------------------------
#' ## Start WGCNA. Choosing a soft thresholding power, Beta.
#-------------------------------------------------------------------------------
#+ eval = FALSE

##############################################################################
# Pick a community!
##############################################################################
# Subset cleanDat based on prots of interest.
n <- 2
group <- c("Shank2", "Shank3", "Syngap1", "Ube3a")[n]
prots <- out[[n]]$proteins
subg_name <- group
subDat <- subset(cleanDat, rownames(cleanDat) %in% prots)
dim(subDat)

##############################################################################

# Estimate powers?
estimatePower <- TRUE

# Load TAMPOR cleanDat from file: #3022 rows.
datafile <- paste(Rdatadir, tissue, "TAMPOR_data_outliersRemoved.Rds", sep = "/")
cleanDat <- readRDS(datafile)
cleanDat <- log2(cleanDat)
cleanDat[1:5, 1:5]
dim(cleanDat) # 1 outlier removed.

# Strip module, metamodule and modulColor attributes from subg.
attr <- names(vertex_attr(subg, index = V(subg)))[c(4, 5, 6)]
for (i in 1:length(attr)) {
  subg <- delete_vertex_attr(subg, attr[i])
}

# Check.
length(names(vertex_attr(subg, index = V(subg)))) == 35

# Load combined sample info.
traitsfile <- paste(Rdatadir, tissue, "Combined_Cortex_Striatum_traits.Rds", sep = "/")
sample_info <- readRDS(traitsfile)
sample_info[1:5, 1:5]
dim(sample_info)

# Allow parallel WGCNA calculations:
allowWGCNAThreads()
nThreads <- 9
clusterLocal <- makeCluster(c(rep("localhost", nThreads)), type = "SOCK")
registerDoParallel(clusterLocal)

## Determine soft power, beta.
# Vector of powers to test:
powers <- seq(4, 20, by = 1.0)

# Soft Power selection
if (estimatePower == TRUE) {
  sft <- pickSoftThreshold(t(subDat),
    powerVector = powers, corFnc = "bicor",
    blockSize = 15000, verbose = 3, networkType = "signed"
  )

  # Create table.
  mytable <- round(sft$fitIndices, 2)
  mytable$truncated.R.sq <- NULL
  table <- tableGrob(mytable, rows = NULL)
  grid.arrange(table)

  # Save table as tiff.
  # file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology_Table.tiff")
  # ggsave(file,table, width = 3, height = 2.5, units = "in")

  # Figure. ggplotScaleFreeFit() generates three plots.
  plots <- ggplotScaleFreeFit(sft)
  plots$Grid

  # Save as tiff.
  # file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology_1.tiff")
  # ggsave(file,plots$ScaleFreeFit, width = 3, height = 2.5, units = "in")
}

# Choose a scale free power (beta).
beta <- sft$fitIndices$Power[sft$fitIndices$SFT.R.sq > 0.8][1]
cat(paste(
  "Power (beta): ",
  beta, "\nScale Free fit: ", fit$stats$scaleFreeRsquared, "\n"
))

# Examine scale free fit...
# Calculate node connectivity in the Weighted co-expression network.
connectivity <- softConnectivity(
  datExpr = t(subDat),
  corFnc = "bicor",
  weights = NULL,
  type = "signed",
  power = beta,
  blockSize = 15000,
  minNSamples = NULL,
  verbose = 0,
  indent = 0
)

# Examine fit.
fit <- ggplotScaleFreePlot(connectivity)
fit$stats
plot <- fit$ggplot
plot

# Save fig.
file <- paste0(outputfigsdir, "/", outputMatName, group, "_WGCNA_ScaleFreeFit.eps")
ggsave(file, plot, width = 3, height = 2.5, units = "in")

collectGarbage()

#-------------------------------------------------------------------------------
#' ## Build the WPCNA Network.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Load blockwiseModules Parameters.
file <- paste(Rdatadir, list.files(Rdatadir, pattern = group), sep = "/")
file <- file[1]
if (!grepl(group, file) | !(length(file) == 1)) {
  print("Warning! Choose the correct sample params file.")
}

# Extract optimal parameters.
sampled_params <- do.call(rbind, sapply(readRDS(file), "[", 1))
rownames(sampled_params) <- paste0("params_", (1:nrow(sampled_params)))
sampled_params$iter <- c(1:nrow(sampled_params))
sampled_params$costRank <- rank(sampled_params$costGrey) / nrow(sampled_params)
sampled_params$modRank <- rank(sampled_params$q2) / nrow(sampled_params)
sampled_params$score <- sampled_params$modRank - sampled_params$costRank
sampled_params <- sampled_params[order(sampled_params$score, decreasing = TRUE), ]

# Choose network building parameters.
params_iter <- sampled_params$iter[1] # top ranked params
params <- sampled_params[params_iter, ]
params[, c(1:7)]

# Network Parameters table.
mytable <- params[, c(1:7)]
table <- tableGrob(mytable, rows = NULL)
# Save
file <- paste0(outputfigsdir, "/", outputMatName, group, "_Network_Parameters.eps")
ggsave(file, table, width = 9, height = 1)

# Network parameters:
beta
networkType <- "signed"
corType <- "bicor"
enforceMMS <- TRUE # Should minimal modules size be inforced? Done after network building.

# Other key parameters for optimization:
minModSize <- params$minModSize # Minimum module size. seq(1,50,by=1)
deepSplit <- params$deepSplit # Sensitivity for module splitting [0-4]. 4 is most sensitive. Increasing results in more modules.
mergeCutHeight <- params$mergeCutHeight # Cut height for module detection. Was 0.07. Increasing results in more modules.
reassignThresh <- params$reassignThresh # pvalue threshold for reassigning nodes to modules.
minKMEtoStay <- params$minKMEtoStay # minimum module connectivity score for assigning to a module.
minCoreKMESize <- params$minCoreKMESize # minimim number of genes in a modules with minKMEtoStay.
pamStage <- params$pamStage # partitioning about medioids
pamStage <- FALSE

# Insure that minCoreKME is not less than minModSize.
if (minCoreKMESize < minModSize) {
  minCoreKMESize <- minModSize
  print(paste0("Warning: minCoreMKESize set to minModsize", " (", minModSize, ")"))
}

# Defaults:
maxBlockSize <- 12000 # maximum block size for module detection.
detectCutHeight <- 0.995 # dendrogram cut height for module detection.

# Call blockwiseModules to build WGCNA network.
# Setting saveTOM = FALSE will really slow thing down.
net <- blockwiseModules(t(subDat),
  power = beta,
  deepSplit = deepSplit,
  minModuleSize = minModSize,
  mergeCutHeight = mergeCutHeight,
  TOMDenom = "mean",
  detectCutHeight = detectCutHeight,
  corType = corType,
  networkType = networkType,
  pamStage = pamStage,
  pamRespectsDendro = TRUE,
  reassignThresh = reassignThresh,
  minCoreKMESize = minCoreKMESize,
  minKMEtoStay = minKMEtoStay,
  verbose = 3,
  saveTOMs = FALSE,
  maxBlockSize = maxBlockSize
)

# Check the number of modules.
nModules_original <- length(unique(net$colors))
nModules_original

collectGarbage()

#-------------------------------------------------------------------------------
#' ## Enforce module preservation.
#-------------------------------------------------------------------------------
#' Remove modules that are not preserved (i.e. have insignificant module
#' preservation statistics). The number of random permutations used to generate the
#' null distributions is increased to 100,000 in order to stabilize the result with
#' large number of modules. This computation is expensive and will take several
#' minutes.
#'
#+ eval = FALSE

# Input for NetRep:
r <- bicor(t(subDat)) # Data has already be log-transformed.
adjm <- ((1 + r) / 2)^beta # Signed network.

data_list <- list(data = t(subDat)) # The protein expression data.
correlation_list <- list(data = r) # The bicor correlation matrix.
network_list <- list(data = adjm) # The weighted, signed co-expresion network.
module_labels <- net$colors # Module labels.
names(module_labels) <- rownames(subDat)

# Try self-preservation test.
preservation <- NetRep::modulePreservation(
  network = network_list,
  data = data_list,
  correlation = correlation_list,
  moduleAssignments = list(data = module_labels),
  modules = NULL,
  backgroundLabel = "grey",
  discovery = "data",
  test = "data",
  selfPreservation = TRUE,
  nThreads - 1,
  # nPerm = 100000, # Increase nPerm to 100,000 in order to stabilize the result with large number of modules.
  null = "overlap",
  alternative = "greater",
  simplify = TRUE,
  verbose = TRUE
)

# Collect stats.
preservation <- preservation[c("observed", "p.values")]

# Get the maximum permutation test p-value.
maxPval <- apply(preservation$p.values, 1, function(x) max(x, na.rm = TRUE))

# Modules removed if adjusted pvalue is greater than alpha = 0.05.
alpha <- 0.05
modules_out <- names(maxPval)[maxPval > alpha / nModules_original]
nModules_out <- length(modules_out)
names(net$MEs)[grepl(paste(modules_out, collapse = "|"), names(net$MEs))] <- "MEgrey"

# Set non-significant modules to grey.
net$colors[net$colors %in% modules_out] <- "grey"

# Total number of modules, excluding grey:
nModules <- length(unique(net$colors)) - 1
params$nModules <- nModules
nModules

# Check percent grey.
percent_grey <- round(100 * sum(net$colors == "grey") / length(net$colors), 2)
print(paste("Percent grey nodes =", percent_grey))
params$PercentGrayNodes <- percent_grey

# Table of key network stats.
mytable <- data.frame(
  nNodes = sum(net$colors != "grey"),
  PercentGrey = round(params$PercentGrayNodes, 2),
  nModules = nModules,
  MedianCoherence = round(params$medianModCoherence, 3)
)
# NetworkModularity = round(params$q2,3))
table <- tableGrob(mytable, rows = NULL)
grid.arrange(table)

# Save.
file <- paste0(outputfigsdir, "/", outputMatName, group, "_Key_Network_Stats.eps")
ggsave(file, table, width = 6, height = 1)

#-------------------------------------------------------------------------------
#' ## Enforce min module size. Recalculate MEs.
#-------------------------------------------------------------------------------

# Fraction of un-assigned proteins.
print(paste("Percent grey = ", round(100 * table(net$colors)["grey"] / length(net$goodGenes), 2),
  " (n=", table(net$colors)["grey"], ")",
  sep = ""
))

# Module Summary (excluding grey)
nModules <- length(table(net$colors)) - 1
print(paste("Total number of modules =", nModules, "(without grey)"))
modules <- cbind(colnames(as.matrix(table(net$colors))), table(net$colors))
orderedModules <- cbind(Mnum = paste("M", seq(1:nModules), sep = ""), Color = labels2colors(c(1:nModules)))
modules <- modules[match(as.character(orderedModules[, 2]), rownames(modules)), ]
moduleColors <- as.data.frame(cbind(orderedModules, Size = modules))

# If necessary, enforce minModSize. Set these to grey.
# Loop to enforce minModSize:
if (enforceMMS) {
  removedModules <- orderedModules[which(modules < minModSize), "Color"]
  print(paste(
    "There are", length(removedModules), "modules that contain less than",
    minModSize, "proteins."
  ))
  for (i in removedModules) {
    net$colors[net$colors == i] <- "grey"
  }
  for (i in removedModules) {
    net$MEs[, paste0("ME", i)] <- NULL
  }
  nModules <- length(table(net$colors)) - 1
  modules <- cbind(colnames(as.matrix(table(net$colors))), table(net$colors))
  orderedModules <- cbind(Mnum = paste("M", seq(1:nModules), sep = ""), Color = labels2colors(c(1:nModules)))
  modules <- modules[match(as.character(orderedModules[, 2]), rownames(modules)), ]
  moduleColors <- as.data.frame(cbind(orderedModules, Size = modules))
}

# Recalculate MEs.
MEs <- moduleEigengenes(t(subDat),
  colors = net$colors,
  impute = FALSE,
  nPC = 1,
  align = "along average",
  excludeGrey = FALSE,
  grey = if (is.numeric(colors)) 0 else "grey",
  subHubs = TRUE,
  trapErrors = FALSE,
  returnValidOnly = FALSE,
  softPower = beta,
  scale = TRUE,
  verbose = 0, indent = 0
)
net$MEs <- MEs$eigengenes

# module assignments.
moduleColors <- data.frame(
  protein = rownames(subDat),
  module = net$colors
)

# Check number of modules.
nModules_original
nModules
nModules_out

#-------------------------------------------------------------------------------
#' ## Examine WPCNA network.
#-------------------------------------------------------------------------------

# Fixme: which geno has poor dendro? # Syngap1 poor dendro in that some colors
# appear to be split up... inspected kme and it looks okay...

# Ube3a some colors split up a little bit...

# Generate dissimilarity matrix (1-TOM(adj)).
# Exclude grey proteins.
idx <- net$colors != "grey"
diss <- 1 - TOMsimilarityFromExpr(
  datExpr = t(subDat[idx, ]),
  corType = "bicor",
  networkType = "signed",
  power = beta,
  TOMType = "signed",
  TOMDenom = "mean",
  verbose = 0
)

# Perform hierarchical clustering.
GeneNames <- rownames(subDat)
colnames(diss) <- rownames(diss) <- GeneNames[idx]
hier <- hclust(as.dist(diss), method = "average") # average, single, complete, median, centroid

# Generate dendrogram with module color labels.
file <- paste0(outputfigsdir, "/", outputMatName, group, "_Network_Dendrogram.eps")
#setEPS()
#postscript(file)
plotDendroAndColors(hier,
  net$colors[idx],
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)

#dev.off()

#-------------------------------------------------------------------------------
#' ## Calculate Modularity of the WPCNA partition.
#-------------------------------------------------------------------------------

# Calculate the adjacency network.
r <- bicor(t(subDat))
adjm <- ((1 + r) / 2)^beta # signed.

# Create igraph object.
graph <- graph_from_adjacency_matrix(
  adjmatrix = adjm,
  mode = c("undirected"),
  weighted = TRUE,
  diag = FALSE
)

# Calculate modularity, q.
membership <- as.numeric(as.factor(net$colors))
q1 <- modularity(graph, membership, weights = edge_attr(graph, "weight"))
q1

# Without "grey" nodes.
v <- rownames(subDat)[!net$colors == "grey"]
subg <- induced_subgraph(graph, v)
membership <- as.numeric(as.factor(net$colors))
membership <- membership[!net$colors == "grey"]
q2 <- modularity(subg, membership, weights = edge_attr(subg, "weight"))
q2

#-------------------------------------------------------------------------------
#' ## Calculate module membership (kME).
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Module membership (kME) is the strength of association (correlation) between a
# proteins expression vector and that of module EigenProteins.

colors <- as.list(net$colors)
names(colors) <- rownames(subDat)

# The protein-wise correlation with module Eigenproteins.
kMEdat <- signedKME(t(subDat), net$MEs, corFnc = "bicor")
colnames(kMEdat) <- gsub("kME", "", colnames(kMEdat))

# Clean up a little.
idx <- match(rownames(kMEdat), rownames(subDat))
kMEdat <- add_column(kMEdat, Protein = rownames(kMEdat), .before = 1)
kMEdat <- add_column(kMEdat, Module = net$colors[idx], .after = 1)
rownames(kMEdat) <- NULL

# Write to excel with custom function write.excel().
# file <- paste0(outputtabsdir,"/",tissue,"_WGCNA_Analysis_ModuleMembership_kME.xlsx")
# write.excel(kMEdat,file)

#-------------------------------------------------------------------------------
#' ## Prepare Numeric metadata for WGCNA analysis.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Insure numericMeta rows and cleanDat columns match
numericMeta <- sample_info
numericMeta <- numericMeta[match(colnames(cleanDat), rownames(numericMeta)), ]
dim(numericMeta)

# Add column for SexType
numericMeta$SexType <- numericMeta$Sex

# Add column for TissueType
numericMeta$TissueType <- numericMeta$Tissue

# Add column for Tissue.Sample.Model
numericMeta$Tissue.Sample.Model <- paste(numericMeta$TissueType,
  numericMeta$Sample.Model,
  sep = "."
)

# Pool WT within a tissue.
idx <- grepl("Cortex.WT", numericMeta$Tissue.Sample.Model)
numericMeta$Tissue.Sample.Model[idx] <- "Cortex.WT"
idy <- grepl("Striatum.WT", numericMeta$Tissue.Sample.Model)
numericMeta$Tissue.Sample.Model[idy] <- "Striatum.WT"
unique(numericMeta$Tissue.Sample.Model)

# Pull out as many numeric traits for correlation later as we can.
numericMeta$Group <- tempVec <- as.vector(as.data.frame(do.call(
  rbind,
  strsplit(colnames(cleanDat), "\\.")
))[, 1])
numericMeta$Group <- as.numeric(as.factor(numericMeta$Group))

# Tissue
numericMeta$Tissue <- as.numeric(as.factor(numericMeta$Tissue)) - 1

# Genotype groupings (genetic background).
numericMeta$Model <- gsub(" ", "", numericMeta$Model)
numericMeta$Syngap1 <- as.numeric(numericMeta$Model == "Syngap1")
numericMeta$Ube3a <- as.numeric(numericMeta$Model == "Ube3a")
numericMeta$Shank2 <- as.numeric(numericMeta$Model == "Shank2")
numericMeta$Shank3 <- as.numeric(numericMeta$Model == "Shank3")

# Sex
numericMeta$Sex <- as.numeric(as.factor(numericMeta$Sex)) - 1

# Make Syngap1 WT v KO column.
temp_df <- as.data.frame(do.call(rbind, strsplit(numericMeta$Sample.Model, "\\.")))
numericMeta$Syngap1_KO <- NA
numericMeta$Syngap1_KO[which(temp_df$V1 == "WT")] <- 0
numericMeta$Syngap1_KO[which(temp_df$V1 == "HET" & temp_df$V2 == "Syngap1")] <- 1

# Make Ube3a WT v HET column.
temp_df <- as.data.frame(do.call(rbind, strsplit(numericMeta$Sample.Model, "\\.")))
numericMeta$Ube3a_KO <- NA
numericMeta$Ube3a_KO[which(temp_df$V1 == "WT")] <- 0
numericMeta$Ube3a_KO[which(temp_df$V1 == "KO" & temp_df$V2 == "Ube3a")] <- 1

# Make Shank2 WT v KO column.
temp_df <- as.data.frame(do.call(rbind, strsplit(numericMeta$Sample.Model, "\\.")))
numericMeta$Shank2_KO <- NA
numericMeta$Shank2_KO[which(temp_df$V1 == "WT")] <- 0
numericMeta$Shank2_KO[which(temp_df$V1 == "KO" & temp_df$V2 == "Shank2")] <- 1

# Make Shank3 WT v KO column.
temp_df <- as.data.frame(do.call(rbind, strsplit(numericMeta$Sample.Model, "\\.")))
numericMeta$Shank3_KO <- NA
numericMeta$Shank3_KO[which(temp_df$V1 == "WT")] <- 0
numericMeta$Shank3_KO[which(temp_df$V1 == "KO" & temp_df$V2 == "Shank3")] <- 1

# Control vs. disease.
numericMeta$SampleType <- gsub(" ", "", numericMeta$SampleType)
numericMeta$ASD <- NA
numericMeta$ASD[numericMeta$SampleType == "WT"] <- 0
numericMeta$ASD[numericMeta$SampleType == "KO"] <- 1
numericMeta$ASD[numericMeta$SampleType == "HET"] <- 1

# Tissue:Genotype groupings (genetic background).
f1 <- as.factor(numericMeta$TissueType)
f2 <- as.factor(numericMeta$Model)
mod <- model.matrix(~ 0 + f1:f2)
colnames(mod) <- apply(expand.grid(levels(f1), levels(f2)), 1, paste, collapse = ".")
numericMeta <- cbind(numericMeta, mod)

# Tissue:Sex
f1 <- as.factor(numericMeta$TissueType)
f2 <- as.factor(numericMeta$SexType)
mod <- model.matrix(~ 0 + f1:f2)
colnames(mod) <- apply(expand.grid(levels(f1), levels(f2)), 1, paste, collapse = ".")
numericMeta <- cbind(numericMeta, mod)

# Generate Tissue Specific Sample.Model contrasts
g <- paste(numericMeta$TissueType, numericMeta$Sample.Model, sep = ".")
allContrasts <- combn(unique(g), 2)
# Keep if tissue and model are equal.
keepers <- function(x) {
  y <- do.call(rbind, strsplit(x, "\\."))
  logic <- y[1, 1] == y[2, 1] & y[1, 3] == y[2, 3]
  return(logic)
}
keep <- apply(allContrasts, 2, function(x) keepers(x))
contrasts <- allContrasts[, keep]
# Pool WT within a tissue.
contrasts[grepl("Cortex.WT", contrasts)] <- "Cortex.WT"
contrasts[grepl("Striatum.WT", contrasts)] <- "Striatum.WT"
# Coerce to list.
contrasts_list <- apply(contrasts, 2, as.list)
# lapply through contrasts list and generate model vector.
mod_list <- lapply(contrasts_list, function(x) match(numericMeta$Tissue.Sample.Model, x) - 1)
mod <- do.call(cbind, mod_list)
# Add names
colnames(mod) <- sapply(contrasts_list, "[", 1)
# Fix Cortex.KO.Shank3 column name.
colnames(mod)[grep("Cortex.WT", colnames(mod))] <- "Cortex.KO.Shank3"

# Add to numericMeta
numericMeta <- cbind(numericMeta, mod)

# Tissue specific disease status (Control v ASD).
g <- numericMeta$SampleType
g[grepl("KO|HET", g)] <- "ASD"
g[grepl("WT", g)] <- "Control"
f1 <- paste(numericMeta$TissueType, g, sep = ".")
allContrasts <- list(
  c("Cortex.ASD", "Cortex.Control"),
  c("Striatum.ASD", "Striatum.Control")
)
mod_list <- lapply(allContrasts, function(x) match(f1, x) - 1)
mod <- do.call(cbind, mod_list)
colnames(mod) <- colnames(mod) <- sapply(allContrasts, "[", 1)
numericMeta <- cbind(numericMeta, mod)

# Sex specific disease status (control v ASD).
g <- numericMeta$SampleType
g[grepl("KO|HET", g)] <- "ASD"
g[grepl("WT", g)] <- "Control"
numericMeta$M.ASD <- as.numeric(numericMeta$SexType == "M" & g == "ASD")
numericMeta$F.ASD <- as.numeric(numericMeta$SexType == "F" & g == "ASD")

# Determine numerical indices. The columns of numericMeta with numerical data.
# Warnings OK; This determines which traits are numeric and if forced to numeric values,
# non-NA values do not sum to 0.
numericIndices <- unique(c(
  which(!is.na(apply(numericMeta, 2, function(x) sum(as.numeric(x))))),
  which(!(apply(numericMeta, 2, function(x) sum(as.numeric(x), na.rm = T))) == 0)
))

#-------------------------------------------------------------------------------
#' ## VerboseBoxplots - vizualize module expression across traits.
#-------------------------------------------------------------------------------

# Fixme: refine this chunk to focus on pairwise comparisons. 
# Fixme: remove unnecesary code.
# Fixme: problem with ggplot build / class NULL

prot2color <- as.list(net$colors)
names(prot2color) <- rownames(subDat)

color2prots <- split(rownames(subDat), net$colors)

# Should plots be saved?
saveplots <- FALSE
savegroups <- FALSE

# Calculate Module EigenProteins.
MEList <- moduleEigengenes(t(subDat), colors = net$colors)
MEs <- orderMEs(MEList$eigengenes)
colnames(MEs) <- gsub("ME", "", colnames(MEs))
rownames(MEs) <- rownames(numericMeta)

# Insure traits are in matching order.
traits <- sample_info
idx <- match(rownames(MEs), rownames(traits))
traits <- traits[idx, ]
all(rownames(traits) == rownames(MEs))

# Define groups, the biological groups of interest.
groups <- paste(numericMeta$TissueType, numericMeta$Sample.Model, sep = ".")
groups[grepl("Cortex.WT", groups)] <- "WT.Cortex"
groups[grepl("Striatum.WT", groups)] <- "WT.Striatum"
unique(groups)

# Calculate Kruskal-Wallis pvalues for all modules (columns of MEs df).
KWtest <- apply(MEs, 2, function(x) kruskal.test(x, as.factor(groups)))

##########################  Test ###########################
# Test just a single genotype... Not looking for convergence!

idx <- groups %in% c("WT.Cortex", "WT.Striatum", unique(groups[grepl(group,groups)]))
g <- as.factor(groups[idx])
x <- MEs$turquoise[idx]
x <- MEs$blue[idx]
x <- MEs$yellow[idx]
x <- MEs$brown[idx]
x <- MEs$green[idx]
#x <- MEs$brown[idx]
#x <- MEs$grey[idx]
kw <- kruskal.test(x, g)
kw
dunnTest(x, g)

##########################  Test ###########################

# Extract pvalues from the list of KW tests.
# The pvalue is the 3rd element of each list.
KW_results <- as.data.frame(do.call(rbind, sapply(KWtest, "[", 3)))

# Clean up the result. Drop grey!
colnames(KW_results) <- "p.value"
rownames(KW_results) <- gsub(".p.value", "", rownames(KW_results))
KW_results <- add_column(KW_results, Module = rownames(KW_results), .before = 1)
rownames(KW_results) <- NULL
KW_results <- KW_results[!KW_results$Module == "grey", ]

## Which comparisons are significant?
# Bonferroni correction for nModules.
KW_results$p.adj <- p.adjust(as.numeric(KW_results$p.value), method = "bonferroni")
# Benjamini-Hochberg FDR correction.
KW_results$FDR <- p.adjust(as.numeric(KW_results$p.value), method = "BH")
KW_results <- KW_results[order(KW_results$p.value), ]

# nsig modules.
KWsigModules <- KW_results$Module[KW_results$FDR < 0.05]
print(paste("nModules with KW FDR < 0.05 =", length(KWsigModules)))

# Split Module EigenProtein (MEs) dm into a list of column vectors for lapply.
ME_list <- split(as.matrix(MEs), rep(1:ncol(MEs), each = nrow(MEs)))
names(ME_list) <- colnames(MEs)

# Add vector of groups
ME_list <- lapply(ME_list, function(x) data.frame(x = x, groups = groups))

# Define levels for order of bars in plot.
# levels <- c("M.WT","M.ASD","F.WT","F.ASD")
levels <- c(
  "WT.Cortex", "WT.Striatum",
  "Cortex.KO.Shank2", "Striatum.KO.Shank2",
  "Cortex.KO.Shank3", "Striatum.KO.Shank3",
  "Cortex.HET.Syngap1", "Striatum.HET.Syngap1",
  "Cortex.KO.Ube3a", "Striatum.KO.Ube3a"
)

# Generate contrasts matrix for comparisons of interest.
contrasts <- makePairwiseContrasts(list("M.WT", "F.WT"), list("M.ASD", "F.ASD"))
g1 <- list("WT.Cortex", "WT.Striatum")
g2 <- list(
  c("Cortex.KO.Shank2", "Cortex.KO.Shank3", "Cortex.HET.Syngap1", "Cortex.KO.Ube3a"),
  c("Striatum.KO.Shank2", "Striatum.KO.Shank3", "Striatum.HET.Syngap1", "Striatum.KO.Ube3a")
)
contrasts <- makePairwiseContrasts(g1, g2)

# Loop through ME_list and generate verboseBoxPlot.
# lapply wont work here because the name is not preserved when you call lapply()...
# method is the p.adj method for the Dunn's test p-value.
plot_data <- list()
for (i in 1:dim(MEs)[2]) {
  x <- ME_list[[i]]$x
  g <- ME_list[[i]]$groups
  color <- names(ME_list)[[i]]
  plot <- ggplotVerboseBoxplot(x, g,
    levels,
    contrasts,
    color,
    stats = TRUE,
    method = "dunn",
    correction_factor = 1
  )
  plot_data[[i]] <- plot
  names(plot_data)[[i]] <- color
}

# Extract plots.
plot_list <- sapply(plot_data, "[", 1)

# Store as VBplots
vbplots <- plot_list

# Extract post-hoc test stats.
Dtest_stats <- sapply(plot_data, "[", 3)
names(Dtest_stats) <- sapply(strsplit(names(Dtest_stats), "\\."), "[", 1)

# Loop to add module column.
for (i in 1:length(Dtest_stats)) {
  df <- Dtest_stats[[i]]
  namen <- names(Dtest_stats)[i]
  df <- add_column(df, Module = namen, .before = 1)
  Dtest_stats[[i]] <- df
}

# Number of significant post-hoc tests.
nsig <- do.call(rbind, lapply(Dtest_stats, function(x) sum(x$P.adj < 0.05)))
idx <- match(KW_results$Module, rownames(nsig))
KW_results$DunnettTestNsig <- nsig[idx, ]

# Convert Dtest stats into data frame. By casting long data into wide.
# dcast should work for multiple value.var, but it isnt...
df <- do.call(
  rbind,
  lapply(
    Dtest_stats,
    function(x) reshape2::dcast(x, Module ~ Comparison, value.var = c("Z"))
  )
)
df2 <- do.call(
  rbind,
  lapply(
    Dtest_stats,
    function(x) reshape2::dcast(x, Module ~ Comparison, value.var = c("P.unadj"))
  )
)

# Fix column names.
colnames(df)[2:ncol(df)] <- paste(colnames(df)[2:ncol(df)], "Z")
colnames(df2)[2:ncol(df2)] <- paste(colnames(df2)[2:ncol(df2)], "P.value")

# Bind as single df. Drop grey.
data <- cbind(df, df2)
data <- data[!data$Module == "grey", ]
dim(data)

# Add KW P.value and FDR.
idx <- match(data$Module, KW_results$Module)
data <- add_column(data, KW.P.value = KW_results$p.value[idx], .after = 1)
data <- add_column(data, KW.FDR = KW_results$FDR[idx], .after = 2)
kw_data <- data

## Save plots...
# Fix names. 
names(vbplots) <- sapply(strsplit(names(vbplots),"\\."),"[", 1)

# Edit plots such that x-axis labels are red if post-hoc test is significant.
for (i in 1:dim(KW_results)[1]) {
  plot <- vbplots[[KW_results$Module[i]]]
  stats <- Dtest_stats[[KW_results$Module[i]]]
  b <- ggplot_build(plot)
  foo <- b$data[[3]]
  sig <- foo$label[order(foo$x)]
  a <- rep("black", 8)
  a[grepl("\\*", sig)] <- "red"
  a <- c(rep("black", 2), a)
  plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
  file <- paste0(outputfigsdir, "/", outputMatName, "_", group, "_", KW_results$Module[i], "_verboseBoxplot", ".tiff")
  ggsave(file, plot, width = 3.25, height = 2.5, units = "in")
}

#-------------------------------------------------------------------------------
#' ## Identify top module Hub proteins.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Should plots be saved?
saveplots <- FALSE

# Modify traits$Sample.Model for grouping all WT's together.
traits <- sample_info
traits_temp <- traits
# Column groupings are determined by traits$Sample.Model.
traits_temp$Sample.Model <- paste(traits_temp$Tissue,
                                  traits_temp$Sample.Model,
                                  sep = "."
)
traits_temp$Sample.Model[grepl("Cortex.WT", traits_temp$Sample.Model)] <- "Cortex.WT"
traits_temp$Sample.Model[grepl("Striatum.WT", traits_temp$Sample.Model)] <- "Striatum.WT"

# Find top n hub proteins.
n <- 3
HubProteins <- melt(kMEdat, id = c("Protein", "Module"))
HubProteins <- subset(HubProteins, HubProteins$Module == HubProteins$variable)
HubProteins$variable <- NULL
colnames(HubProteins) <- c("Protein", "Module", "kME")
HubProteins <- HubProteins %>%
  group_by(Module) %>%
  top_n(n, kME)

# Generate boxplots to examine expresion of hub proteins.
plot_list <- ggplotProteinBoxes(cleanDat,
                                interesting.proteins = as.character(HubProteins$Protein),
                                dataType = "Relative Abundance", traits = traits_temp,
                                order = c(2, 7, 4, 6, 5, 8, 1, 9, 3, 10)
)

# Example plot.
plot_list[[1]]

# Build a df with statistical results.
# colnames should match x axis of plot.
stats <- lapply(results, function(x)
  as.data.frame(cbind(Uniprot = x$Uniprot, FDR = x$FDR)))
names(stats) <- names(results)
df <- stats %>% reduce(left_join, by = "Uniprot")
colnames(df)[c(2:ncol(df))] <- names(stats)

# Annotate rows as gene|uniprot
Uniprot <- df$Uniprot
Gene <- mapIds(
  org.Mm.eg.db,
  keys = Uniprot,
  column = "SYMBOL",
  keytype = "UNIPROT",
  multiVals = "first"
)
rownames(df) <- paste(Gene, Uniprot, sep = "|")
df$Uniprot <- NULL
stats <- df

# Define proteins of interest.
prots <- HubProteins$Protein
length(prots)

# Generate boxplots to examine expresion of hub proteins.
plot_list <- ggplotProteinBoxes(cleanDat,
                                interesting.proteins = as.character(prots),
                                dataType = "Relative Abundance", traits = traits_temp,
                                order = c(2, 7, 4, 6, 5, 8, 1, 9, 3, 10)
)
# Example plot.
plot <- plot_list[[1]]
plot
annotate_stars(plot, stats)

# Loop to add stars.
plot_list <- lapply(plot_list, function(x) annotate_stars(x, stats))

# Loop to add custom colors.
colors <- rep(c("gray", "yellow", "blue", "green", "purple"), each = 2)
plot_list <- lapply(plot_list, function(x) x + scale_fill_manual(values = colors))

# Example plot:
plot_list[[1]]

# Loop to add module assignment annotation.
#results_modules <- data.frame(
#  geneNames = rownames(cleanDat),
#  dynamicColors = net$colors
#)
#for (i in 1:length(plot_list)) {
#  plot <- plot_list[[i]]
#  namen <- results_modules$dynamicColors[match(plot$labels$title, results_modules$geneNames)]
#  plot <- plot + ggtitle(paste(plot$labels$title, namen, sep = "|"))
#  plot_list[[i]] <- plot
#}

# Store as list.
mhplots <- split(plot_list, HubProteins$Module)

# Print to pdf (all plots).
#if (saveplots == TRUE) {
#  file <- paste0(outputfigsdir, "/", tissue, "_WGCNA_Analysis_HubProtein_Expression.pdf")
#  ggsavePDF(plot_list, file)
#}

## Write hub data to excel.
#file <- paste0(outputtabsdir, "/", tissue, "_WGCNA_Analysis_Module_Hubs.xlsx")
#write.excel(as.data.frame(HubProteins), file)

# Function to loop and save top three hubs of desired color.
#my_func <- function(color) {
#  for (i in 1:3) {
#    plot <- mhplots[[color]][[i]]
#    file <- paste0(outputfigsdir, "/", outputMatName, color, "_", i, "_Module_Hub.tiff")
#    ggsave(file, plot, width = 3, height = 2.5, units = "in")
#  }
#}

# Save plots:
if (saveplots == TRUE) {
for (i in 1:length(mhplots)) {
  color <- names(mhplots)[i]
  for (x in 1:length(mhplots[[i]])) {
    plot <- mhplots[[i]][[x]]
    file <- paste0(outputfigsdir, "/", outputMatName, group, "_", color, x, "_Module_Hub.tiff")
    ggsave(file, plot, width = 3, height = 3, units = "in")
    }
}
}

#-------------------------------------------------------------------------------
#' Send DEP Community with Co-expression modules to Cytoscape.
#-------------------------------------------------------------------------------

######################### INSURE YOU HAVE SELECTED THE CORRECT GENOTYPE ########
# Which subgraph?
N <- 2
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")
geno <- groups[N]
file <- paste0(Rdatadir, "/", "DEP_KNN_Communities.Rds")
KNN_communities <- readRDS(file)
prots <- sapply(KNN_communities, "[", 1)
names(prots) <- groups
prots <- prots[[geno]]
length(prots)
################################################################################

# Load TAMPOR cleanDat from file: #3022 rows.
datafile <- paste(Rdatadir, tissue, "TAMPOR_data_outliersRemoved.Rds", sep = "/")
cleanDat <- readRDS(datafile)
cleanDat <- log2(cleanDat)
cleanDat[1:5, 1:5]
dim(cleanDat) # 1 outlier removed.

# Subset cleanDat based on prots of interest.
subg_name <- geno
subDat <- subset(cleanDat, rownames(cleanDat) %in% prots)
dim(subDat)

# Load the data.
dir <- "D:/Documents/R/Synaptopathy-Proteomics/Tables/Network"
file <- paste(dir, "SIF.xlsx", sep = "/")
sif <- read_excel(file, sheet = 1)

# Create a data frame with all node attributes.
nodes <- data.frame(
  Entrez = unlist(meta$entrez[!is.na(meta$entrez)]),
  Uniprot = unlist(meta$uniprot[!is.na(meta$entrez)]),
  Symbol = unlist(meta$gene[!is.na(meta$entrez)])
)

# Gather DEP stats.
dat <- stats_df[, grep(geno, colnames(stats_df))]
rownames(dat) <- stats_df$Entrez
dat$sigProt <- apply(dat, 1, function(x) any(x < 0.05))
idx <- match(nodes$Entrez, rownames(dat))
nodes$sigProt <- dat$sigProt[idx]
sum(nodes$sigProt == TRUE)

# Make igraph object.
g <- graph_from_data_frame(d = sif, vertices = nodes, directed = FALSE)

# Coerce to simple graph--remove duplicate edges and self-loops.
g <- simplify(g)
is.simple(g)

# Number of nodes and edges.
length(V(g))
length(E(g))

# Map prots to entrez.
v <- unlist(lapply(prots, function(x) unlist(protein2entrez[x])))

# Subset the graph.
subg <- induced_subgraph(g, v)

# Send graph to cytoscape?
send_to_cytoscape <- TRUE

# Build df of node attributes.
idx <- match(unlist(entrez2protein[names(V(subg))]), rownames(subDat))
hex <- lapply(as.list(net$colors[idx]), function(x) col2hex(x))
df <- data.frame(
  Protein = unlist(entrez2protein[names(V(subg))]),
  Entrez = names(V(subg)),
  Color = net$colors[idx],
  HexColor = unlist(hex)
)
rownames(df) <- df$Entrez

# Send to cytoscape with RCy3!
if (send_to_cytoscape == TRUE) {
  cytoscapePing()
  quiet(RCy3::createNetworkFromIgraph(subg, subg_name))

  # Load node attribute table in cytoscape.
  loadTableData(df)

  # Create custom syle to customize appearance.
  style.name <- subg_name
  colvec <- df$HexColor

  defaults <- list(
    NODE_FILL = col2hex("grey"),
    NODE_SHAPE = "Ellipse",
    NODE_SIZE = 55,
    EDGE_WIDTH = 2.0,
    EDGE_TRANSPARENCY = 120
  )

  nodeLabels <- mapVisualProperty("node label", "Symbol", "p")
  nodeFills <- mapVisualProperty("node fill color", "HexColor", "p", colvec)

  # edgeWidth <- mapVisualProperty('edge width','weight','p')
  createVisualStyle(style.name, defaults, list(nodeLabels, nodeFills))
  lockNodeDimensions(TRUE, style.name)
  setVisualStyle(style.name)

  # Apply perfuse force directed layout.
  layoutNetwork(layout.name = "force-directed")
  setNodeColorDefault(col2hex("grey"), style.name)
}

#-------------------------------------------------------------------------------
#' ## GO Enrichment analysis of Modules (colors).
#-------------------------------------------------------------------------------
#+ eval = FALSE

## Prepare a matrix of class labels (colors) to pass to enrichmentAnalysis().
geneNames <- rownames(subDat)
dynamicColors <- net$colors
results_modules <- as.data.frame(cbind(geneNames, dynamicColors))

# Reshape the module colors data.
colors <- as.data.frame(results_modules)
mytable <- table(colors)
colors <- as.data.frame.matrix(mytable)
head(colors)

# Convert 0 and 1 to column names.
logic <- colors == 1 # 1 will become TRUE, and 0 will become FALSE.
# Loop through each column to replace 1 with column header (color).
for (i in 1:ncol(logic)) {
  col_header <- colnames(colors)[i]
  colors[logic[, i], i] <- col_header
}

## Map Uniprot IDs to Entrez.
# Get Uniprot IDs and gene symbols from rownames
# Uniprot_IDs <- as.character(rownames(colors))
Uniprot_IDs <- as.character(colsplit(rownames(colors), "\\|", c("Symbol", "UniprotID"))[, 2])
# Map Uniprot IDs to Entrez IDs:
entrez <- mapIds(org.Mm.eg.db,
  keys = Uniprot_IDs,
  column = "ENTREZID",
  keytype = "UNIPROT",
  multiVals = "first"
)

# Insure that colors is a matrix.
colors <- as.matrix(colors)
head(colors)

# look at the number of genes assigned to each cluster.
table(colors)

# The colors matrix and vector of cooresponding entrez IDs
# will be passed to enrichmentAnalysis().

# Build a GO annotation collection:
if (!exists(deparse(substitute(musGOcollection)))) {
  musGOcollection <- buildGOcollection(organism = "mouse")
}

# Perform GO analysis for each module using hypergeometric (Fisher.test) test.
# As implmented by the WGCNA function enrichmentAnalysis().
# FDR is the BH adjusted p-value.
# Insure that the correct background (used as reference for enrichment)
# has been selected!
# useBackgroud = "given" will use all given genes as reference background.

GOenrichment <- enrichmentAnalysis(
  classLabels = colors,
  identifiers = entrez,
  refCollection = musGOcollection,
  useBackground = "given", # options are: given, reference (all), intersection, and provided.
  threshold = 0.05,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "0"
)

# Save GOenrichment as RDS
# file <- paste(Rdatadir,tissue,"GOenrichment_data.RDS",sep="/")
# saveRDS(GOenrichment,file)

# Create some space by clearing some memory.
collectGarbage()

# Collect the results.
results_GOenrichment <- list()
for (i in 1:length(GOenrichment$setResults)) {
  results_GOenrichment[[i]] <- GOenrichment$setResults[[i]]$enrichmentTable
}
length(results_GOenrichment)

# Combined result
GO_result <- do.call(rbind, results_GOenrichment)
idx <- GO_result$Bonferroni < 0.05
nsig <- length(idx[idx == TRUE])
print(paste("There are", nsig, "GO terms that exhibit significant enrichment among all modules."))

# Top Term for each module.
topGO <- subset(GO_result, rank == 1)

# Number of modules with sig. GO enrichment.
table(topGO$FDR < 0.1)[2]
table(topGO$Bonferroni < 0.1)[2]
sigGO <- subset(topGO, FDR < 0.1)

# These are the unique sig. GO terms.
unique(sigGO$shortDataSetName)

# Add names to list of results.
names(results_GOenrichment) <- colnames(colors)
names(results_GOenrichment)

# Add topGo to list.
results_GOenrichment <- c(list(topGO = topGO), results_GOenrichment)
names(results_GOenrichment)

# Write results to file.
# file <- paste0(outputtabsdir,"/",tissue,"_WGCNA_Analysis_Module_GOenrichment_Results.xlsx")
# write.excel(results_GOenrichment,file)

#-------------------------------------------------------------------------------
#' ## Visualize GO enrichemnt. 
#-------------------------------------------------------------------------------

saveplots <- FALSE

# Single plot.
#ggplotGOscatter(results_GOenrichment, color = "brown", topN = 5)

# All plots.
colors <- as.list(unique(net$colors))
plot_list <- lapply(colors, function(x) ggplotGOscatter(results_GOenrichment, x, topN = 5))
names(plot_list) <- colors
goplots <- plot_list

# Save top modules as tiff.
if (saveplots==TRUE) {
for (i in 1:length(plot_list)) {
  color <- names(plot_list)[i]
  file <- paste0(outputfigsdir, "/", outputMatName, "GO_Scatter_", group, "_", color, ".tiff")
  ggsave(file, plot_list[[color]])
}
}

#-------------------------------------------------------------------------------
#' ## Build complete GO semantic similarity graph.
#-------------------------------------------------------------------------------
# Build a GO semantic similarity graph or load it from file.

# Choose a GO ontology.
type <- 1 # MF, BP, CC
ontology <- c("MF","BP","CC")[type]

# Evaluate GO similarity for all genes.
build_GO_similarity_network <- FALSE # otherwise it will be loaded from file. 

if (build_GO_similarity_network == TRUE){
  # Build GO database.
  if (!exists("msGOMF")){ msGOMF <- godata('org.Mm.eg.db', ont= "MF")}
  if (!exists("msGOBP")){ msGOBP <- godata('org.Mm.eg.db', ont= "BP")}
  if (!exists("msGOCC")){ msGOCC <- godata('org.Mm.eg.db', ont= "CC")}
  msGO <- list(msGOMF,msGOBP,msGOCC)[[type]]
  # Build GO similarity matrix. 
  goSim <- mgeneSim(genes = unlist(meta$entrez),
                    semData = msGO, measure="Wang",verbose=TRUE)
  file <- paste0(Rdatadir,"/","New_GO_similarity_network","_",ontology,".Rds")
  saveRDS(goSim,file)
}else{
  # Read from file. 
  print(paste("Loaded GO",ontology, "network from file!"))
  file <- paste0(Rdatadir,"/","New_GO_similarity_network","_",ontology,".Rds")
  goSim <- readRDS(file)
}

#-------------------------------------------------------------------------------
#' ## Test preservation of module's biological coherence with NetRep.
#-------------------------------------------------------------------------------

# NetRep input...
# GO Similarity matrix.
adjm <- goSim
dim(adjm)
# Remap adjm names to symbol|uniprot.
idx <- match(rownames(adjm),unlist(meta$entrez))
rownames(adjm) <- meta$protein[idx]
colnames(adjm) <- meta$protein[idx]

# Insure adjm and data have matching dimensions.
idx <- rownames(subDat) %in% rownames(adjm)
newDat <- subDat[idx,]
colors <- net$colors[idx]
idx <- rownames(adjm) %in% rownames(subDat)
adjm <- adjm[idx,idx]

# Check.
dim(newDat)[1] == dim(adjm)[1] & length(colors) == dim(adjm)[1]

# Sort.
idx <- match(rownames(adjm), rownames(newDat))
newDat <- newDat[idx,]
colors <- colors[idx]
names(colors) <- rownames(newDat)

# Check.
all(rownames(newDat)==rownames(adjm))
adjm[1:5,1:5] # What about values == 1?
dim(adjm)[1]

# Correlation matrix.
r <- bicor(t(newDat)) # Data has already be log-transformed.

# NetRep input:
data_list <- list(data = t(newDat))  # The protein expression data. 
correlation_list <- list(data = r)   # The bicor correlation matrix. 
network_list <- list(data = adjm)    # The go similarity matrix.   
module_labels <- colors              # Module labels. 
names(module_labels) <- rownames(newDat)

# Calculate module statistics. 
preservation <- NetRep::modulePreservation(
  network = network_list, 
  data = data_list, 
  correlation = correlation_list, 
  moduleAssignments = list(data = module_labels),
  modules = NULL, 
  backgroundLabel = "grey", 
  discovery = "data", 
  test = "data",
  selfPreservation = TRUE, 
  #nThreads-1, 
  #nPerm = 0, 
  null = "overlap", 
  alternative = "greater", 
  simplify = TRUE,
  verbose = TRUE)

# Collect stats. 
preservation <- preservation[c("observed","p.values")]
df <- data.frame(p.value = preservation$p.values[,1]) # just average edge weight.
df$p.adjust <- p.adjust(df$p.value, method = "bonferroni")

# Are modules biologically cohesive?
df # Shank3, yes. Syngap1 no. Ube3a no... 
