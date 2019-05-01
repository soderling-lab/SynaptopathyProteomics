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
})

# Define version of the code.
CodeVersion <- "params"

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

# Globally set ggplots theme.
ggplot2::theme_set(theme_gray())

#-------------------------------------------------------------------------------
#' ## Load data.
#-------------------------------------------------------------------------------

# Data is...
# Load TAMPOR cleanDat from file:
datafile <- paste(Rdatadir,tissue,"TAMPOR_data_outliersRemoved.Rds",sep="/")
cleanDat <- readRDS(datafile)
cleanDat <- log2(cleanDat)
cleanDat[1:5,1:5]
dim(cleanDat)

# Directory for PPI data.
dir <- paste(datadir,"PPI Network", sep="/")

# Load SIF.
file <- paste(dir,"Cortex_SIF_031019.txt", sep="/")
sif <- read.table(file,header=TRUE, sep=",")

#-------------------------------------------------------------------------------
#' ## Evaluate modularity of PPI graph based on WGCNA partitions.
#-------------------------------------------------------------------------------

# Load colors.
# Load Parameters and network stats.
files <- list.files(Rdatadir, pattern = "Stats")
file <- paste0(Rdatadir,"/",files)
file
out <- readRDS(file)
length(out)

# Network stats
result <- do.call(rbind,sapply(out,"[",1))
result$iter <- c(1:nrow(result))
rownames(result) <- paste0("params_",1:nrow(result))

# Module colors.
modcolors <- sapply(out,"[",2)
names(modcolors) <- rownames(result)

# Preservation stats.
modstats <- sapply(out,"[",3)
names(modstats) <- rownames(result)

# Build data frame.
df <- data.frame(Protein = rownames(cleanDat),
                 Color = NA)
# Map uniprot IDs to Entrez.
df$Uniprot <- sapply(strsplit(df$Protein,"\\|"),"[",2)
df$Gene <- sapply(strsplit(df$Protein,"\\|"),"[",1)
df$Entrez <- mapIds(org.Mm.eg.db, keys=df$Uniprot, column="ENTREZID", 
                    keytype="UNIPROT", multiVals="first")
df$Entrez2 <- mapIds(org.Mm.eg.db, keys=df$Gene, column="ENTREZID", 
                     keytype="SYMBOL", multiVals="first")

# Function to get non-NA value as Entrez ID. 
my_func <- function(x){
  if (!is.na(x[5])){
    y = x[5]
  }else{
    y = x[6]
  }
}
df$id <- apply(df,1,function(x) my_func(x))

out <- list()
# Loop to calculate modularity.
for (i in 1:length(modcolors)){
  print(i)
  # Add colors.
  df$Color <- modcolors[[i]]

  # Remove sif rows that are not mapped to genes in df.
  sif$logic <- sif$EntrezA %in% df$id & sif$EntrezB %in% df$id
  sif <- sif[sif$logic,]

  # Make igraph.
  sif_df <- sif[,c(3,4)]
  g <- graph_from_data_frame(sif_df, directed = FALSE)

  nodes <- vertex_attr(g, "name")
  colors <- df$Color[match(nodes,df$id)]

  # Calculate modularity, q.
  membership <- as.numeric(as.factor(colors))
  q1 <- modularity(g, membership)

  # Without "grey" nodes.
  v <- nodes[!colors=="grey"]
  subg <- induced_subgraph(g,v)
  membership <- as.numeric(as.factor(colors))
  membership <- membership[!colors=="grey"]
  q2 <- modularity(subg, membership)
  dm <- cbind(q1,q2)
  out[[i]] <- dm
}

result <- as.data.frame(do.call(rbind,out))
result$params <- c(1:nrow(result))
sub <- result[result[,2]==max(result[,2]),]

#-------------------------------------------------------------------------------
#' ## PPI graph of WGCNA modules.
#-------------------------------------------------------------------------------

# Build data frame.
df <- data.frame(Protein = rownames(cleanDat),
                 Color = net$colors)
# Map uniprot IDs to Entrez.
df$Uniprot <- sapply(strsplit(df$Protein,"\\|"),"[",2)
df$Gene <- sapply(strsplit(df$Protein,"\\|"),"[",1)
df$Entrez <- mapIds(org.Mm.eg.db, keys=df$Uniprot, column="ENTREZID", 
                    keytype="UNIPROT", multiVals="first")
df$Entrez2 <- mapIds(org.Mm.eg.db, keys=df$Gene, column="ENTREZID", 
                     keytype="SYMBOL", multiVals="first")

# Function to get non-NA value as Entrez ID. 
my_func <- function(x){
  if (!is.na(x[5])){
    y = x[5]
  }else{
    y = x[6]
  }
}
df$id <- apply(df,1,function(x) my_func(x))

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

# Sig modules.
sigMods <- module_summary$Module[module_summary$p.adj<0.1]
modules <- modules[names(modules)[names(modules) %in% sigMods]]

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
  setNodeColorDefault(col2hex(color)) # Must be hex.
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
# V#isualize WPCNA communities. 
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
