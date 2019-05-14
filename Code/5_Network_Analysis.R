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
#' ## Load data.
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

# Add gene name.
meta$gene <- mapIds(org.Mm.eg.db, keys=meta$entrez, column="SYMBOL", 
                    keytype="ENTREZID", multiVals="first")

# Load TAMPOR statistical results.
file <- paste(outputtabs,"Final_TAMPOR",
              "Combined_TMT_Analysis_TAMPOR_GLM_Results.xlsx", sep = "/")
results <- lapply(as.list(c(1:8)),function(x) read_excel(file,x))
names(results) <- excel_sheets(file)

# Load convergent modules.
file <- paste(Rdatadir,"module_overlap.Rds", sep="/")
module_overlap <- readRDS(file)

#-------------------------------------------------------------------------------
#' ## Build SIF and NOA files.
#-------------------------------------------------------------------------------

# Add Tissue.Genotype column names to results.
cols <- do.call(rbind,strsplit(names(results),"\\."))[,-2]
col_names <- apply(cols, 1, function(x) paste(x,collapse="."))

# Loop to rename columns of dataframes in list. 
for (i in 1:length(results)){
  df <- results[[i]]
  y <- col_names[i]
  colnames(df)[c(4:ncol(df))] <- paste(y,colnames(df)[c(4:ncol(df))])
  results[[i]] <- df
}

# Bind dataframes in list together.
df_NOA <- results %>% reduce(left_join, by = c("Uniprot","Entrez","Gene"))

# Remove candiate columns.
df_NOA <- df_NOA[,-grep("candidate",colnames(df_NOA))]

# Percentage of genes mapped to sif.
table(df_NOA$Entrez %in% sif$EntrezA | df_NOA$Entrez %in% sif$EntrezB)[2]/
  length(unique(c(sif$EntrezA,sif$EntrezB)))

# Annotate as Cortex.Sig, Genotype.Sig, and Tissue.Genotype.Sig
idx <- grep("FDR",colnames(df_NOA))
logic <- df_NOA[,idx] < 0.05
ids <- paste(sapply(strsplit(colnames(logic)," "),"[",1),"sig",sep = ".")

# Loop to convert TRUE to Tissue.Genotype.Sig
out <- list()

# Loop through each column to replace 1 with column header (color).
for (i in 1:ncol(logic)){
  col_header <- ids[i]
  temp <- logic[,i]
  temp[temp] <- col_header
  out[[i]] <- temp
}

# Bind togehter, replace FALSE, add column names. 
out <- do.call(cbind,out)
out[out==FALSE] <- "NA"
colnames(out) <- paste0("cat",c(1:ncol(out)))

# Add to NOA table. 
df_NOA <- cbind(df_NOA,out)

# Export NOA and sif
file <- paste(outputtabsdir,"NOA.csv",sep="/")
write.csv(df_NOA,file)

file <- paste(outputtabsdir,"SIF.csv",sep="/")
write.csv(sif,file)

## Add WGCNA modules to NOA.
gene <- sapply(strsplit(rownames(cleanDat),"\\|"),"[",1)
uniprot <- sapply(strsplit(rownames(cleanDat),"\\|"),"[",2)

# Map Uniprot to entrez.
entrez <- mapIds(org.Mm.eg.db, keys=uniprot, column="ENTREZID", 
                 keytype="UNIPROT", multiVals="first")
table(is.na(entrez))

# Map un-mapped genes to entrez.
entrez[is.na(entrez)] <- mapIds(org.Mm.eg.db, keys=gene[is.na(entrez)], 
                                column="ENTREZID", keytype="SYMBOL", multiVals="first")
table(is.na(entrez)) # Good only 3 remaining ids are not mapped.

# Check, meta and cleanDat in matching order.
all(meta$protein == rownames(cleanDat))

# Add meta modules.
df_NOA2 <- data.frame(Uniprot = uniprot,
                 Entrez = entrez,
                 Gene = gene,
                 Module = net$colors,
                 MetaModule = meta$metaModule)

# Change meta Module ids to MM#
df_NOA2$MetaModule <- paste0("MM",df_NOA2$MetaModule)

# Add hex colors for modules.
colors <- unlist(lapply(as.list(df_NOA2$Module),function(x) col2hex(x)))
df_NOA2$Module_color <- colors

# Save to csv.
file <- paste(outputtabsdir,"NOA2.csv",sep="/")
write.csv(df_NOA2,file)

# Add column for interacton type to sif.
sif <- add_column(sif,type = "ppi",.after=ncol(sif))

# Export this as csv.
file <- paste(outputtabsdir,"SIF.csv",sep="/")
write.csv(sif,file)

sif2 <- data.frame(NodeA = df_NOA2$Gene,
                   NodeB = df_NOA2$Module,
                   EntrezA = df_NOA2$Entrez,
                   EntrezB = df_NOA2$Module,
                   type = "module-gene")

# Bind and export as csv.
file <- paste(outputtabsdir,"SIF2.csv",sep="/")
write.csv(sif2,file)

# Build a WPCE network.
data <- cleanDat
all(rownames(data)==meta$protein)
rownames(data) <- meta$entrez
# Omit un-mapped (NA) rows.
data <- data[!is.na(rownames(data)),] 
r <- bicor(t(data))
adjm <- ((1+r)/2)^12 # Signed network.

# Create igraph object.
graph <- graph_from_adjacency_matrix(
  adjmatrix = adjm, 
  mode = c("undirected"), 
  weighted = TRUE, 
  diag = FALSE)

#-------------------------------------------------------------------------------
#' ## Find first degree neighbors with 2+ degree to seed nodes.
#-------------------------------------------------------------------------------

# Create a Dictionary-like object of genes and entrez ids.
entrez <- as.list(meta$gene)
genes <- as.list(meta$entrez)
names(entrez) <- meta$entrez
names(genes) <- meta$gene

# Create igraph object.
g <- graph_from_data_frame(sif[,c(3,4,1,2)], directed = FALSE)

seeds <- list()
network_size <- list()
degree_to_stay <- 2

# Loop to create networks.
for (i in 1:length(module_overlap)){
  
  # Get subset of nodes (v) in modules overlap.
  # We will use these to seed a network.
  v <- meta$entrez[meta$module %in% module_overlap[[i]]]
  
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
  # Number of nodes in each network.
  #unlist(lapply(subg,function(x) length(V(x))))
  #sum(unlist(lapply(subg,function(x) length(V(x)))))

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
  
  # Calculate closeness centrality.
  cc <- closeness(uniong, 
                  vids = V(uniong), 
                  mode = "all",
                  weights = NULL, 
                  normalized = TRUE)
  # Add to dist.
  dist$closeness <- cc[match(rownames(dist),names(cc))]
  
  # Calculate degree to seed nodes (sum).
  dist$SeedDegree <- apply(dist[,-ncol(dist)],1,function(x) sum(x))
  
  # We will keep nodes that have at least 2 degrees with seed nodes.
  keep <- dist$SeedDegree>=degree_to_stay
  dist <- dist[keep,]
  
  # Rank by closeness centrality. 
  dist$ccRank <- rank(dist$closeness)/nrow(dist)
  
  # Keep top 25%
  dist <- dist[order(dist$ccRank,decreasing = TRUE),]
  n <- round(.25 * nrow(dist))
  sub <- dist[c(1:n),]
  keepers <- unique(c(rownames(sub),v))
  subg <- induced_subgraph(g,keepers)
  network_size[[i]] <- length(V(subg))
  
  # Send to cytoscape.
  cytoscapePing()
  print(paste("Working on subgraph", i,"..."))
  quiet(RCy3::createNetworkFromIgraph(subg,names(module_overlap)[i]))
  if (i==length(module_overlap)){print("Complete!")}
}

# Number of seed nodes for each graph.
result <- data.frame(Seed_Number = unlist(lapply(seeds,function(x) length(x))),
                     Network_Size = unlist(network_size))

#-------------------------------------------------------------------------------
#' ## PPI graph of WGCNA meta Module communities.
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
