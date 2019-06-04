#' ---
#' title: PPI Network Analysis. 
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

# Sometimes, if you have not cleared the workspace of all loaded packages,
# you man incounter problems.
# To remove all packages, you can call the following:
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
CodeVersion <- "Network_Analysis"

# Define tisue type: cortex = 1; striatum = 2.
type <- 3
tissue <- c("Cortex", "Striatum", "Combined")[type]

# Set the working directory.
rootdir <- "D:/Documents/R/Synaptopathy-Proteomics"
setwd(rootdir)

# Set any other directories.
functiondir <- paste(rootdir, "Functions", sep = "/")
datadir <- paste(rootdir, "Input", sep = "/")
Rdatadir <- paste(rootdir,"RData", sep = "/")

# Create code-version specific figure and tables folders if they do not already exist.
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
functiondir <- paste(rootdir, "Functions", sep = "/")
my_functions <- paste(functiondir, "TMT_Preprocess_Functions.R", sep = "/")
source(my_functions)

# Define prefix for output figures and tables.
outputMatName <- paste(tissue, "_Network_Analysis_", sep = "")

# Globally set ggplots theme.
ggplot2::theme_set(theme_gray())

#-------------------------------------------------------------------------------
#' ## Load the WGNCA and TMT data.
#-------------------------------------------------------------------------------

# Data is...
# Load TAMPOR cleanDat from file:
datafile <- paste(Rdatadir,tissue,"TAMPOR_data_outliersRemoved.Rds",sep="/")
cleanDat <- readRDS(datafile)
cleanDat <- log2(cleanDat)
cleanDat[1:5,1:5]
dim(cleanDat)

# Load PPI data
dir <- paste(datadir,"PPI Network", sep="/")
file <- paste(dir,"Cortex_SIF_031019.txt", sep="/")
sif <- read.table(file,header=TRUE, sep=",")

# Load WGCNA network and meta Modules.
file <- paste(Rdatadir,"Network_and_metaModules.Rds",sep="/")
data <- readRDS(file)
net <- data$net
meta <- data$meta

# Add gene names.
meta$gene <- mapIds(org.Mm.eg.db, keys=meta$entrez, column="SYMBOL", 
                    keytype="ENTREZID", multiVals="first")

# Add uniprot.
meta$uniprot <- sapply(strsplit(meta$protein,"\\|"),"[",2)

# Load TAMPOR statistical results.
file <- paste(outputtabs,"Final_TAMPOR",
              "Combined_TMT_Analysis_TAMPOR_GLM_Results.xlsx", sep = "/")
results <- lapply(as.list(c(1:8)),function(x) read_excel(file,x))
names(results) <- excel_sheets(file)

# Load convergent modules.
file <- paste(Rdatadir,"module_overlap.Rds", sep="/")
module_overlap <- readRDS(file)

#-------------------------------------------------------------------------------
#' ## Load compiled PPI network
#-------------------------------------------------------------------------------

# Load the data.
dir <- "D:/Documents/R/Synaptopathy-Proteomics/Tables/Network"
file <- paste(dir,"SIF.xlsx",sep="/")
sif <- read_excel(file,sheet = 1)

# Node attributes. 
nodes <- data.frame(Entrez = unlist(meta$entrez[!is.na(meta$entrez)]),
                    Symbol = unlist(meta$gene[!is.na(meta$entrez)]),
                    Module = unlist(meta$module[!is.na(meta$entrez)]),
                    MetaModule = unlist(meta$metaModule[!is.na(meta$entrez)]))
nodes$MetaModule <- paste0("MM",nodes$MetaModule)

# Add hex colors for modules.
colors <- unlist(lapply(as.list(meta$module),function(x) col2hex(x)))
colors <- colors[!is.na(meta$entrez)]
nodes$ModulColor <- colors

# Make igraph object. 
g <- graph_from_data_frame(d=sif, vertices=nodes, directed=FALSE)

# Coerce to simple graph--remove duplicate edges and self-loops.
g <- simplify(g)
is.simple(g)

# Number of nodes and edges. 
length(V(g)) # All but three unmapped genes. 
length(E(g))

# Number of connected components.
connected_components <- components(g, mode = c("weak", "strong"))
connected_components$csize

#-------------------------------------------------------------------------------
#' Significantly dysregulated proteins...
#-------------------------------------------------------------------------------

# Build a df with statistical results.
stats <- lapply(results,function(x) 
  as.data.frame(cbind(Uniprot=x$Uniprot,FDR=x$FDR)))
names(stats) <- names(results)
df <- stats %>% reduce(left_join, by = "Uniprot")
colnames(df)[c(2:ncol(df))] <- names(stats)
idx <- match(df$Uniprot,meta$uniprot)
df <- add_column(df,Entrez = meta$entrez[idx], .after = 1)

# Proteins with any significant change.
df$sigProt <- apply(df,1,function(x) any(x[c(3:10)]<0.05))

sum(df$sigProt) # Total
round(100*sum(df$sigProt)/nrow(df),3) # Percent

# Create table for cytoscape.
file <- paste(outputtabsdir,"SigProts.csv",sep = "/")
write.csv(df,file)

#-------------------------------------------------------------------------------
#' ## Generate module subgraphs.
#-------------------------------------------------------------------------------

# Send graphs to cytoscape?
send_to_cytoscape <- FALSE

# Generate subgraphs for all modules.
module_membership <- split(meta$entrez,meta$module)

# Empty lists for output of loop
module_subraphs <- list()
cc <- list()

# Loop:
n <- length(module_membership)

# Open Cytoscape application.
#system('"C:/Program Files/Cytoscape_v3.6.1/Cytoscape.exe"', wait = FALSE)
if (send_to_cytoscape == TRUE){ cytoscapePing() }

for (i in 1:n){
  
  # Subset of nodes.
  # Insure that all are in graph.
  v <- module_membership[[i]]
  v <- v[v %in% names(V(g))]
  
  # Make subgraph
  subg <- induced_subgraph(g,v)
  module_subraphs[[i]] <- subg
  
  # Calculate subg's mean clustering coefficent or transitivity.
  # This is a measure of a graphs interconnectedness. 
  cc[[i]] <- mean(transitivity(subg, type = "local", isolates = "zero"))
  
  # Send to cytoscape.
  if (send_to_cytoscape == TRUE){
    print(paste("Working on subgraph", i,"..."))
    # Change vertex names to gene symbol. 
    idx <- match(vertex_attr(subg)$name,meta$entrez)
    subg <- set.vertex.attribute(subg,"name",value = meta$gene[idx])
    color <- col2hex(names(module_membership)[i])
    setNodeShapeDefault('Ellipse')
    lockNodeDimensions(TRUE)
    quiet(RCy3::createNetworkFromIgraph(subg,names(module_membership)[i]))
  }
}

# Gather clustering coefficients. 
result <- data.frame(Module = names(module_membership),
                     MeanClusteringCoefficient = unlist(cc))

result$EdgeDensity <- unlist(
  lapply(module_subraphs,function(x) length(E(x))/length(V(x))))

#-------------------------------------------------------------------------------
#' ## Evaluate GO semantic similarity (~biological cohesiveness) for all modules 
#-------------------------------------------------------------------------------
# Utilize the GOSemSim package to create GO similarity matrixes for all module 
# subgraphs.

# Build GO database.
if (!exists("msGOMF")){ msGOMF <- godata('org.Mm.eg.db', ont= "MF")}
if (!exists("msGOBP")){ msGOBP <- godata('org.Mm.eg.db', ont= "BP")}
if (!exists("msGOCC")){ msGOCC <- godata('org.Mm.eg.db', ont= "CC")}

msGO <- list(msGOMF,msGOBP,msGOCC)
names(msGO) <- c("MF","BP","CC")

# Empty list for output of loop.
out <- list()
n <- length(module_membership)

# Loop:
for (i in 1:n){
  print(paste0("Working on subgraph: ", i, "..."))
  goSim <- lapply(msGO, function(x)
    mgeneSim(genes = module_membership[[i]], 
             semData = x, 
             measure="Wang",
             verbose=FALSE))
  
  # Calculate average edge weight.
  my_func <- function(x){
    diag(x) <- NA
    aew <- mean(x,na.rm=TRUE)
    return(aew)
  }
  foo <- lapply(goSim, function(x) my_func(x))  
  out[[i]] <- foo
  
  # Evaluate centrality as mean clustering coefficient in the GO similarity graph. 
  #subg <- lapply(goSim, function(x) 
  #  graph_from_adjacency_matrix(x, mode = "undirected", weighted = TRUE, diag = FALSE))
  #cc <- lapply(subg, function(x)
  #  transitivity(x, type = "average", isolates = "zero"))
  #out[[i]] <- cc
  }

# Extract the results...
GO_Similarity <- data.frame(
  MF = unlist(sapply(out,"[", 1)),
  BP = unlist(sapply(out,"[", 2)),
  CC = unlist(sapply(out,"[", 3)))

rownames(GO_Similarity) <- names(module_membership)

GO_Similarity$Module <- names(module_membership)
GO_Similarity$ModuleCC <- result$MeanClusteringCoefficient

# Examine the relationship between module size and function coherence.
module_sizes <- as.data.frame(table(net$colors))
idx <- match(GO_Similarity$Module,module_sizes$Var1)

df <- data.frame(Module = GO_Similarity$Module,
                 GOCoherence = GO_Similarity$MF,
                 Size = module_sizes$Freq[idx])
df <- df[!df$Module=="grey",]

plot(df$GOCoherence,df$Size)
cor(df$GOCoherence,df$Size,method = "spearman")

#-------------------------------------------------------------------------------
#' ## Build complete GO semantic similarity graph.
#-------------------------------------------------------------------------------
# Build a GO semaintic similarity graph or load it from file.

# Build GO database.
if (!exists("msGOMF")){ msGOMF <- godata('org.Mm.eg.db', ont= "MF")}
if (!exists("msGOBP")){ msGOBP <- godata('org.Mm.eg.db', ont= "BP")}
if (!exists("msGOCC")){ msGOCC <- godata('org.Mm.eg.db', ont= "CC")}

# Choose a GO ontology.
type <- 1 # MF, BP, CC
ontology <- c("MF","BP","CC")[type]
msGO <- list(msGOMF,msGOBP,msGOCC)[[type]]

# Evaluate GO similarity for all genes.
build_GO_similarity_network <- FALSE

if (build_GO_similarity_network == TRUE){
  goSim <- mgeneSim(genes = unlist(module_membership),
                    semData = msGO, measure="Wang",verbose=TRUE)
  file <- paste0(Rdatadir,"/","GO_similarity_network","_",ontology,".Rds")
  saveRDS(goSim,file)
}else{
  # Read from file. 
  print(paste("Loaded GO",ontology, "network from file!"))
  file <- paste0(Rdatadir,"/","GO_similarity_network","_",ontology,".Rds")
  goSim <- readRDS(file)
}

#-------------------------------------------------------------------------------
#' ## Modularity of the GO similarity graph.
#-------------------------------------------------------------------------------

# The GO semantic similarity matrix.
goSim[1:5,1:5]

# Remove grey. 
out <- subset(meta$entrez,meta$module=="grey")
idx <- rownames(goSim) %in% out
dm <- goSim[!idx,!idx]

gs <- graph_from_adjacency_matrix(dm, mode = "undirected", weighted = TRUE)

# Vector v of node module membership.
nodes <- get.vertex.attribute(gs,"name")
idx <- match(nodes,meta$entrez)
v <- as.numeric(as.factor(meta$module[idx]))

# Modularity...
modularity(gs, membership = v, weights = edge_attr(gs, "weight"))

#-------------------------------------------------------------------------------
# Test preservation of module biological coherence.
#-------------------------------------------------------------------------------
#  Consider a permutation test with 99 permutations, 
#  5 reaching the statistically significant threshold. 
#  Assume a two group experiment with 6 in each group.
#  Input total.nperm=462

# NetRep input:
# GO Similarity matrix.
adjm <- goSim
idx <- match(rownames(adjm),meta$entrez)
rownames(adjm) <- meta$protein[idx]
colnames(adjm) <- meta$protein[idx]

# Insure adjm and data have matching dimensions.
idx <- rownames(cleanDat) %in% rownames(adjm)
subDat <- cleanDat[idx,]

# Sort subDat so that it matches adjm.
idx <- match(rownames(adjm),rownames(subDat))
subDat <- subDat[idx,]
all(rownames(subDat)==rownames(adjm))

# Correlation matrix.
r <- bicor(t(subDat)) # Data has already be log-transformed.

# Module labels.
idx <- match(rownames(subDat),meta$protein)
colors <- meta$module[idx]

data_list <- list(data = t(subDat))  # The protein expression data. 
correlation_list <- list(data = r)   # The bicor correlation matrix. 
network_list <- list(data = adjm)    # The go similarity matrix.   
module_labels <- colors              # Module labels. 
names(module_labels) <- rownames(subDat)

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
df <- data.frame(p.value = preservation$p.values[,1])
df$p.adjust <- p.adjust(df$p.value, method = "bonferroni")

#-------------------------------------------------------------------------------
#' ## Find first degree neighbors with 2+ degree to seed nodes.
#-------------------------------------------------------------------------------

# Create a Dictionary-like object of genes and entrez ids.
entrez <- as.list(meta$gene)
genes <- as.list(meta$entrez)
names(entrez) <- meta$entrez
names(genes) <- meta$gene
head(entrez,2); head(genes,2)

seeds <- list()
network_size <- list()
degree_to_stay <- 2
send_to_cytoscape = FALSE

# Loop to create networks.
# Cytoscape should be open before proceeding!
for (i in 2:length(module_overlap)){
  
  # Get subset of nodes (v) in modules overlap.
  # We will use these to seed a network.
  v <- meta$entrez[meta$module %in% module_overlap[[i]]]
  length(v)
  # Insure that all nodes are in the network.
  #table(v %in% vertex_attr(g, "name"))
  v <- v[v %in% vertex_attr(g, "name")]
  seeds[[i]] <- v
  names(seeds)[[i]] <- names(module_overlap)[[i]]
  
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
  out <- as.character(genes[grep("Ywha*",names(genes))])
  dist[rownames(dist) %in% out,] <- 0
  
  # Calculate degree to seed nodes (sum).
  dist$SeedDegree <- apply(dist[,-ncol(dist)],1,function(x) sum(x))
  
  # We will keep nodes that have at least 2 connections with seed nodes.
  keep <- dist$SeedDegree>=degree_to_stay
  dist <- dist[keep,]
  keepers <- unique(c(v,rownames(dist)))
  subg <- induced_subgraph(g,keepers)
  network_size[[i]] <- length(V(subg))
  
  # Send to cytoscape.
  if (send_to_cytoscape == TRUE){
    cytoscapePing()
    print(paste("Working on subgraph", i,"..."))
    quiet(RCy3::createNetworkFromIgraph(subg,names(module_overlap)[i]))
  }
  # Status message. 
  if (i==length(module_overlap)){print("Complete!")}
}

# Number of seed nodes for each graph.
result <- data.frame(Seed_Number = unlist(lapply(seeds,function(x) length(x))),
                     Network_Size = unlist(network_size))

#-------------------------------------------------------------------------------
#' ## PPI graphs of DEP's from each genotype:tissue.
#-------------------------------------------------------------------------------

# Generate PPIs graphs using DEPs from each genotype as seed nodes.
# Add nodes with 2+ connections to these seed nodes for context.
# Do not consider connections to 1433 proteins.

# Create a Dictionary-like object mapping uniprot IDs to Entrez.
uniprot <- as.list(meta$entrez)
names(uniprot) <- meta$uniprot
genes <- as.list(meta$entrez)
names(genes) <- meta$gene

# Defaults for analysis. 
degree_to_stay <- 2
send_to_cytoscape = FALSE
combine_tissues = FALSE

# Empty list for results.
results <- list()

# Define significantly DEPs.
sigProts <- list()
for (i in 1:length(stats)){
  df <- stats[[i]]
  sigProts[[i]] <- df$Uniprot[df$FDR < 0.05]
}
names(sigProts) <- names(stats)

# Combine tissues for each genotype.
if (combine_tissues == TRUE){
  genos <- c("Shank2","Shank3","Syngap1","Ube3a")
  u <- list()
  for (geno in genos){
    idx <- grep(geno,names(sigProts))
    u[[geno]] <- union(sigProts[[idx[1]]],sigProts[[idx[2]]])
    }
  sigProts <- u
}

# Map sigProts to Entrez.
sigEntrez <- lapply(sigProts,function(x) unlist(uniprot[x]))

# Loop to create networks.
# Cytoscape should be open before proceeding!
for (i in 1:length(sigEntrez)){
  
  print(paste("Working on subgraph", i,"..."))
  
  # Get subset of nodes (v) in modules overlap.
  # We will use these to seed a network.
  v <- sigEntrez[[i]]
  length(v)
  
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
  out <- as.character(genes[grep("Ywha*",names(genes))])
  dist[rownames(dist) %in% out,] <- 0
  
  # Calculate degree to seed nodes (sum).
  dist$SeedDegree <- apply(dist[,-ncol(dist)],1,function(x) sum(x))
  
  # We will keep nodes that have at least 2 connections with seed nodes.
  # degree_to_stay = 2
  keep <- dist$SeedDegree>=degree_to_stay
  dist <- dist[keep,]
  keepers <- unique(c(v,rownames(dist)))
  subg <- induced_subgraph(g,keepers)
  
  # Build df of node attributes. 
  df <- data.frame(Node = names(V(subg)),
                   sigProt = names(V(subg)) %in% sigEntrez[[i]])
  rownames(df) <- names(genes)[match(df$Node,genes)]
  
  # Change vertex names to gene symbol. 
  idx <- match(names(V(subg)),meta$entrez)
  subg <- set.vertex.attribute(subg,"name",value = meta$gene[idx])
  
  # How many nodes. 
  print(paste0("Total Nodes: ", length(V(subg))))
  
  # result
  results[[i]] <- list(seeds= seeds, 
                       subg = subg, 
                       nodes = names(V(subg)))

  # Send to cytoscape.
  if (send_to_cytoscape == TRUE){
    cytoscapePing()
    quiet(RCy3::createNetworkFromIgraph(subg,names(sigEntrez)[i]))
    # Load node attribute table in cytoscape.
    loadTableData(df)
    setNodeShapeDefault('Ellipse')
    lockNodeDimensions(TRUE)
    
  }
}

# Examine overlap in ppi networks...
names(results) <- names(sigProts)
nodes <- sapply(results,"[", 3)
names(nodes)

# Build a matrix showing overlap.
col_names <- names(nodes)
row_names <- names(nodes)

# All possible combinations.
contrasts <- expand.grid(col_names,row_names)
colnames(contrasts) <- c("ConditionA","ConditionB")
contrasts$ConditionA <- as.vector(contrasts$ConditionA)
contrasts$ConditionB <- as.vector(contrasts$ConditionB)

# Loop to calculate intersection for all contrasts. 
int <- list()
for (i in 1:dim(contrasts)[1]){
  a <- unlist(nodes[contrasts$ConditionA[i]])
  b <- unlist(nodes[contrasts$ConditionB[i]])
  int[[i]] <- intersect(a,b)
}

contrasts$Name <- paste(contrasts$ConditionA,
                        contrasts$ConditionB,sep="_U_")
names(int) <- contrasts$Name

# Add intersection to contrasts.
contrasts$Intersection <- lapply(int,function(x) length(x))

# Calculate percent intersection.
int <- list()
for (i in 1:dim(contrasts)[1]){
  a <- unlist(nodes[contrasts$ConditionA[i]])
  b <- unlist(nodes[contrasts$ConditionB[i]])
  int[[i]] <- length(intersect(a,b))/length(unique(c(a,b)))
}

# Add to dm
contrasts$Percent <- unlist(int)

# Make overlap matrix.
dm <- matrix(contrasts$Intersection,nrow=8,ncol=8)
rownames(dm) <- colnames(dm) <- row_names

# heirarchical clustering of dissimilarity matrix calculated as 1-percent_overlap.
dm <- matrix(contrasts$Percent,nrow=8,ncol=8)
rownames(dm) <- colnames(dm) <- row_names
diss <- 1 - dm
hc <- hclust(as.dist(diss), method = "average")
p1 <- ggdendrogram(hc, rotate=FALSE)
p1

# Remove upper tri and melt.
dm[lower.tri(dm)] <- NA

# Calculate percent overlap.
dm2 <- sweep(dm,1,apply(dm, 1, function(x) max(x, na.rm=TRUE)), FUN = "/")

# Melt
df <- melt(dm,na.rm = TRUE)
df$percent <- round(melt(dm2, na.rm = TRUE)$value,2)

# Generate plot.
plot <- ggplot(df, aes(Var2, Var1, fill = percent)) +
  geom_tile(color = "white") + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 2.5) +
  scale_fill_gradient2(name="Percent Overlap") + 
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
    legend.direction = "horizontal") + 
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  coord_fixed()

plot

#-------------------------------------------------------------------------------
#' ## PPI graph of WPCNA modules...
#-------------------------------------------------------------------------------

# Build data frame.
df <- data.frame(Protein = meta$protein,
                 Color = meta$metaModule)
# Map uniprot IDs to Entrez.
df$Uniprot <- sapply(strsplit(df$Protein,"\\|"),"[",2)
df$Gene <- sapply(strsplit(df$Protein,"\\|"),"[",1)
df$Entrez <- mapIds(org.Mm.eg.db, keys=df$Uniprot, column="ENTREZID", 
                    keytype="UNIPROT", multiVals="first")
df$Entrez2 <- mapIds(org.Mm.eg.db, keys=df$Gene, column="ENTREZID", 
                     keytype="SYMBOL", multiVals="first")
head(df)

# Function to get non-NA value as Entrez ID. 
my_func <- function(x){
  if (!is.na(x[5])){
    y = x[5]
  }else{
    y = x[6]
  }
}
df$id <- apply(df,1,function(x) my_func(x))
df$moduleColor <- net$colors[match(df$Protein,rownames(cleanDat))]


# Subset sif based on proteins in df.
sif$logic <- sif$EntrezA %in% df$id & sif$EntrezB %in% df$id
sif <- sif[sif$logic,]

# Build PPI graph.
ppi_df <- sif[,c(3,4,1,2)]
g <- graph_from_data_frame(ppi_df,directed = FALSE)

# Subset of graph based on module colors.
modules <- split(df$id,df$Color)
v <- modules$darkgrey

# Insure all of v are in g
v <- v[v %in% vertex_attr(g,"name")]
subg <- induced_subgraph(g,v)
# Change vertex names to gene symbol. 
idx <- match(vertex_attr(subg)$name,df$id)
subg <- set.vertex.attribute(subg,"name",value = df$Gene[idx])
plot(subg)

# Send module subgraphs to cytoscape. 
cytoscapePing()

for (i in 1:length(modules)){
  print(paste("Working on subgraph",i,"..."))
  v <- modules[[i]]
  color <- names(modules)[i]
  # Insure all of v are in g
  v <- v[v %in% vertex_attr(g,"name")]
  subg <- induced_subgraph(g,v)
  # Change vertex names to gene symbol. 
  idx <- match(vertex_attr(subg)$name,df$id)
  subg <- set.vertex.attribute(subg,"name",value = df$Gene[idx])
  #setNodeColorDefault(col2hex(color)) # Must be hex.
  setNodeShapeDefault('Ellipse')
  lockNodeDimensions(TRUE)
  quiet(RCy3::createNetworkFromIgraph(subg,color))
}




#-------------------------------------------------------------------------------

# Load Gen Louvain communities.
file <- paste(dir,"Cortex_GLouvain_Communities.csv",sep="/")
pcomm <- read.csv(file)
colnames(pcomm) <- c("Protein","Community")

# Annotate with Uniprot and Entrez.
Uniprot <- mapIds(
  org.Mm.eg.db,
  keys=pcomm$Protein,
  column="UNIPROT", 
  keytype="SYMBOL", 
  multiVals="first")

Entrez <- mapIds(
  org.Mm.eg.db,
  keys=pcomm$Protein,
  column="ENTREZID", 
  keytype="SYMBOL", 
  multiVals="first")

pcomm <- add_column(pcomm,Uniprot,.after=1)
pcomm <- add_column(pcomm,Entrez,.after=2)

#-------------------------------------------------------------------------------
# Visualize WPCNA communities. 
#-------------------------------------------------------------------------------
#sif

# Create igraph object.
df <- sif[,c(1,2)]
g <- graph_from_data_frame(df,directed = FALSE)

# Get all proteins in a module.
color <- "steelblue"
prots <- subset(results_modules,dynamicColors==color)
v <- as.character(sapply(strsplit(prots$geneNames,"\\|"),"[",1))
# Keep only those in the graph.
v <- v[v %in% V(g)$name]
print(paste(length(v),"of the", dim(prots)[1], 
            color,"proteins identified in the PPI graph."))

# Subset graph
subg <- induced.subgraph(g,v)
plot(subg)

# Convert to Cytoscape graph object. 
# Insure that a connection to Cytoscape has been established. 
# Cytoscape must be open.
# Multiple graphs can be created.
cytoscapePing()
#setVisualStyle('default')
setNodeColorDefault(col2hex(color)) # Must be hex.
setNodeShapeDefault('Ellipse')
lockNodeDimensions(TRUE)
quiet(RCy3::createNetworkFromIgraph(subg,color))

# Get all proteins in a module.
color <- "brown"
prots <- subset(results_modules,dynamicColors==color)
v <- as.character(prots$geneNames)
length(v)

# Keep only those in the graph.
v <- v[v %in% V(PPIgraph)$name]
length(v)

# Subset PPIgraph.
subg <- induced.subgraph(PPIgraph,v)
plot(subg)

# Convert to Cytoscape graph object. 
# Insure that a connection to Cytoscape has been established. 
# Cytoscape must be open.
# Multiple graphs can be created.
cytoscapePing()
#setVisualStyle('default')
setNodeColorDefault(col2hex(color)) # Must be hex.
setNodeShapeDefault('Ellipse')
lockNodeDimensions(TRUE)
quiet(RCy3::createNetworkFromIgraph(subg,color))
exportImage("Fooman",'PDF')

#-------------------------------------------------------------------------------
# Is any community disproportionately affected?
# Asses enrichment with the hypergeometric test.

# Annotated with sig or NS.
pcomm$sig <- as.numeric(pcomm$Uniprot %in% sig_genes[[1]])

stats <- subset(pcomm) %>% group_by(Community) %>% 
  dplyr::summarise(Sig.Count = sum(sig),
                   Size = length(sig))
# Add N, total protein number and m, total number of sig prots.
stats$N <- sum(stats$Size)
stats$m <- sum(stats$Sig.Count)

# Use apply to calc hypergeometric pvalue for enrichment.
stats$pval <- apply(stats,1, 
                    function(x) phyper(x[2], x[5], x[4]-x[5], x[3], lower.tail = FALSE))
stats$p.adj <- nrow(stats)*stats$pval


# Add FE column.
stats <- add_column(stats,stats$Sig.Count/(((stats$m/stats$N)*stats$Size)),.after=5)
colnames(stats)[6] <- "FE"

####
# Load Protein complexes.
file <- paste(dir,"Cortex_ClusterOne_Complexes.csv",sep="/")
pcomp <- read.csv(file)
colnames(pcomp) <- c("Protein","Complex")

# Annotate with Uniprot and Entrez.
Uniprot <- mapIds(
  org.Mm.eg.db,
  keys=pcomp$Protein,
  column="UNIPROT", 
  keytype="SYMBOL", 
  multiVals="first")

Entrez <- mapIds(
  org.Mm.eg.db,
  keys=pcomp$Protein,
  column="ENTREZID", 
  keytype="SYMBOL", 
  multiVals="first")

pcomp <- add_column(pcomp,Uniprot,.after=1)
pcomp <- add_column(pcomp,Entrez,.after=2)

##Filter protein complexes.
pcomp <- pcomp[!pcomp$Complex==0,]
# Remove small complexes (size<3).
size <- subset(pcomp) %>% group_by(Complex) %>% dplyr::summarise(Size=length(Complex))
pcomp$Size <- size[match(pcomp$Complex,size$Complex),2]
pcomp <- pcomp[!pcomp$Size<3,]

# Hypergeometric test.
pcomp$sig <- as.numeric(pcomp$Uniprot %in% sig_genes[[1]])

stats <- subset(pcomp) %>% group_by(Complex) %>% 
  dplyr::summarise(Sig.Count = sum(sig),
                   Size = length(sig))
# Add N, total protein number and m, total number of sig prots.
stats$N <- sum(stats$Size)
stats$m <- sum(stats$Sig.Count)

# Apply to calc pvalue.
stats$pval <- apply(stats,1, 
                    function(x) phyper(x[2], x[5], x[4]-x[5], x[3], lower.tail = FALSE))
stats$p.adj <- nrow(stats)*stats$pval
stats <- as.data.frame(stats[order(stats$pval),])
pcomp <- as.data.frame(pcomp)

#fixme: size.size column is wack.
pcomp[pcomp$Complex==318,]


#-------------------------------------------------------------------------------
#' ## Which modules are enriched for DEPs?
#-------------------------------------------------------------------------------

# Melt pvalues.
x <- melt(t(stats))
colnames(x) <- c("Condition","Protein","FDR")
x$Module <- results_modules$dynamicColors[match(x$Protein,results_modules$geneNames)]

p <- x$FDR
g <- x$Module

df <- as.data.frame(cbind(g,p))

df.stats <- subset(df) %>% group_by(g) %>% 
  dplyr::summarise(
    Module.Size = length(p),
    Sig.Count = sum(p<0.05))
df.stats$Total.sig <- sum(df.stats$Sig.Count)
df.stats$N <- sum(df.stats$Module.Size)

df.stats$p.val <- phyper(df.stats$Sig.Count,
                         df.stats$Total.sig,
                         df.stats$N-df.stats$Sig.Count,
                         df.stats$Module.Size, lower.tail = FALSE)
df.stats$FDR <- p.adjust(df.stats$p.val,method = "bonferroni")
df.stats$FE.pval <- -log(df.stats$p.val)

# Rank by pvalue.
df.stats <- df.stats[order(df.stats$p.val),]
df.stats$g <- factor(df.stats$g,levels=unique(df.stats$g))

# Generate plot.
plot <- ggplot(df.stats, aes(x=g,y=FE.pval,fill=g)) + geom_col() +
  scale_fill_manual(values=levels(df.stats$g)) + xlab(NULL) + 
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
    axis.title.x = element_text(color = "black", size = 10),
    axis.title.y = element_text(color = "black", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none")

plot + annotate("text", x = df.stats$g, y = 1, 
                label = df.stats$Sig.Count, size = 3)
#fixme: what is the NA?

# Save to pdf.
file <- paste0(outputfigsdir,"/",outputMatName,"Module_DEP_BarPlot.pdf")
ggsavePDF(plot,file)

#-------------------------------------------------------------------------------
#' ## Cytoscape networks with RCy3.
#-------------------------------------------------------------------------------

# The adjacency network.
r <- bicor(t(cleanDat))
net_bicor <- ((1+r)/2)^power
net_bicor[1:5,1:5]

# Create igraph object of protein co-expression (PCEgraph). 
PCEgraph <- graph_from_adjacency_matrix(
  adjmatrix = net_bicor, 
  mode = c("directed"), 
  weighted = TRUE, 
  diag = FALSE)

# Get all proteins in a module.
color <- "brown"
prots <- subset(results_modules,dynamicColors==color)
v <- as.character(prots$geneNames)
length(v)

# Get hubs.
prots <- subset(HubProteins$Protein,HubProteins$Module==color)
v <- as.character(prots)

# Subset GCEgraph.
subg <- induced.subgraph(GCEgraph,v)
plot(subg)

# Convert to Cytoscape graph object. 
# Insure that a connection to Cytoscape has been established. 
# Cytoscape must be open.
# Multiple graphs can be created.
cytoscapePing()
#setVisualStyle('default')
setNodeColorDefault(col2hex(color)) # Must be hex.
setNodeShapeDefault('Ellipse')
lockNodeDimensions(TRUE)
quiet(RCy3::createNetworkFromIgraph(subg,color))

# Get all proteins in a module.
color <- "brown"
prots <- subset(results_modules,dynamicColors==color)
v <- as.character(prots$geneNames)
length(v)

# Keep only those in the graph.
v <- v[v %in% V(PPIgraph)$name]
length(v)

# Subset PPIgraph.
subg <- induced.subgraph(PPIgraph,v)
plot(subg)

# Convert to Cytoscape graph object. 
# Insure that a connection to Cytoscape has been established. 
# Cytoscape must be open.
# Multiple graphs can be created.
cytoscapePing()
#setVisualStyle('default')
setNodeColorDefault(col2hex(color)) # Must be hex.
setNodeShapeDefault('Ellipse')
lockNodeDimensions(TRUE)
quiet(RCy3::createNetworkFromIgraph(subg,color))
exportImage("Fooman",'PDF')

#-------------------------------------------------------------------------------
#' ## Which modules are enriched for PPIs?
#-------------------------------------------------------------------------------

df <- as.data.frame(cbind(Protein = rownames(cleanDat), 
                          Module = net$colors))
# Keep only those in the graph.
df <- df[df$Protein %in% V(PPIgraph)$name,]

# Subset PPI graph.
v2 <- V(PPIgraph)$name
v2 <- v2[v2 %in% df$Protein]
PPIgraph <- induced.subgraph(PPIgraph,v2)

# Total degree.
sum(degree(PPIgraph))

# Split into lists.
v <- split(df,df$Module)

# Subset PPIgraph by module color.
subg_list <- lapply(v,function(x) induced.subgraph(PPIgraph,x$Protein))

# Calculate the number of PPIs (degrees) in a subgraph.
deg_list <- lapply(subg_list,function(x) sum(degree(x)))

do.call(rbind,deg_list)
