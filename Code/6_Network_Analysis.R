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

# Map 3 un-mapped entrez IDs by hand.
meta[is.na(meta$entrez),]
not_mapped = list("Ndufb1|P0DN34" = 102631912,
                  "F8a1|Q00558"   = 14070,
                  "Pc|Q05920"     = 18563)
for (i in 1:length(not_mapped)){
  idx <- meta$protein == names(not_mapped)[i]
  meta$entrez[idx] <- not_mapped[[i]]
}
# Check. 
sum(is.na(meta$entrez))

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

# keep only nodes in sif
#nodes <- subset(nodes, nodes$Entrez %in% unique(c(sif$musEntrezA,sif$musEntrezB)))
#dim(nodes)

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

# Calculate single metric for node color/size.
# Score = log2( PC * -log10(pvalue))
fc <- as.matrix(df[,seq(3,ncol(df),by=3)])
pval <- as.matrix(df[seq(4,ncol(df),by=3)])
score <- fc*(1-pval)
colnames(score) <- gsub("logFC","score",colnames(score))
df <- cbind(df,score)

# Combine with node attributes. 
# use uniprot to combine as this is the MOST stable identifier for all proteins.
nodes <- merge(nodes, df, by = c("Uniprot"))

# Clean-up the df.
nodes$Entrez.y <- NULL
names(nodes)[2] <- "Entrez"
nodes <- nodes[,c(2,1,3:ncol(nodes))] # reorder so that entrez is first col
dim(nodes)

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
connected_components$csize[1] # The largest connected component. 

#-------------------------------------------------------------------------------
#' ## Build subraphs for DEPs and their communities. 
#-------------------------------------------------------------------------------

# Generate PPIs graphs using DEPs from each tissue:genotype contrast as seed 
# nodes. Add nodes with 2+ connections to these seed nodes for biological 
# context. Do not consider connections to 1433* chaperone proteins.

# Create a Dictionary-like object mapping uniprot IDs to Entrez.
uniprot <- as.list(meta$entrez)
names(uniprot) <- meta$uniprot
genes <- as.list(meta$entrez)
names(genes) <- meta$gene
uniprot2gene <- as.list(meta$gene)
names(uniprot2gene) <- meta$uniprot

# Defaults for analysis. 
degree_to_stay = 2
send_to_cytoscape = FALSE # This will only be done for 1 randomly seeded subg. 
combine_tissues   = FALSE
generate_random   = FALSE # 1,000 randomlly seeded subgs will be generated for each DEP community.

# Custom colors:
colors <- as.list(c("#FFF200", "#00A2E8", "#22B14C", "#A349A4"))
names(colors) <- c("Shank2","Shank3","Syngap1","Ube3a")

# Empty list for results.
community_results <- list()
rand_results <- list()

# Define significantly DEPs.
sigProts <- list()
for (i in 1:length(stats)){
  df <- stats[[i]]
  sigProts[[i]] <- df$Uniprot[df$FDR < 0.05]
}
names(sigProts) <- names(stats)

# Add Combined.KO.Shank2
sigProts$Combined.KO.Shank2 <- unique(c(sigProts$Cortex.KO.Shank2,
                                        sigProts$Striatum.KO.Shank2))

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

# Reorder
#sigEntrez <- sigEntrez[c(2,6,3,7,1,5,4,8)]
#names(sigEntrez)

# Loop to create networks.
# Cytoscape should be open before proceeding!
for (i in 1:length(sigEntrez)){
  
  print(paste("Working on subgraph", i,"..."))
  
  # Get subset of nodes (v) in modules overlap.
  # We will use these to seed a network.
  v <- sigEntrez[[i]] # number of significant proteins for this exp.
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
  print(paste("Additional nodes:", length(keepers) - length(seeds)))
  subg <- induced_subgraph(g,keepers)
  
  # Build df of node attributes. 
  idx <- nodes$Entrez %in% names(V(subg))
  df <- nodes[idx,]
  rownames(df) <- df$Entrez
  idy <- grepl(names(sigProts)[i],colnames(df)) & grepl("FDR",colnames(df))
  df$sigProt <- df[,idy] < 0.05
  
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
    geno <- strsplit(names(sigEntrez)[i], "\\.")[[1]][3]
    style.name <- names(sigEntrez)[i]
    colvec <- c(col2hex("gray"), unlist(colors[geno]))
    
    defaults <- list(NODE_FILL = col2hex("grey"),
                     NODE_SHAPE="Ellipse",
                     NODE_SIZE=55,
                     EDGE_WIDTH = 2.0,
                     EDGE_TRANSPARENCY=120)
     
    nodeLabels <- mapVisualProperty('node label','Symbol','p')
    nodeFills <- mapVisualProperty('node fill color','sigProt','d',c(FALSE,TRUE), colvec)
    
    #edgeWidth <- mapVisualProperty('edge width','weight','p')
    createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills))
    lockNodeDimensions(TRUE, style.name)
    setVisualStyle(style.name)
    
    # Apply perfuse force directed layout. 
    layoutNetwork(layout.name = "force-directed")
    setNodeColorDefault(col2hex("grey"), style.name)
    
    # Apply node color mapping to score. (this will overwrite other mapping)                    
    column <- paste("score",names(sigProts)[i])
    control.points <- c(min(df[column]), 0.0, max(df[column]))
    cols <- c(col2hex("blue"), col2hex("white"), col2hex("red"))
    setNodeColorMapping(column, control.points, cols, style.name = style.name)
  }
}

# Name the results.  
names(community_results) <- names(sigEntrez)

# Save network images...
#full.path=paste(getwd(),'vignette_image',sep='/')
#exportImage(full.path, 'PDF') #.pdf

# Convert pdf to tiff...
#pdf <- "D:/Documents/R/Synaptopathy-Proteomics/Tables/Combined/Network_Analysis/Cortex.HET.Syngap1.pdf"
#library(pdftools)
#tiff <- pdf_render_page(pdf, dpi = 600, numeric = TRUE)

# save to bitmap formats
#library(tiff)
#writeTIFF(tiff, "test.tiff", compression = "none")
#-------------------------------------------------------------------------------

# Try breaking down modules by functional groups...


#-------------------------------------------------------------------------------
#' ## Examine overlap between DEP communities and randomly seeded graphs.
#-------------------------------------------------------------------------------

# Examine overlap between subgraphs and random subgraphs. 
names(results) <- names(sigEntrez)
if (length(rand_results>0)){ names(rand_results) <- names(sigEntrez)}

# collect nodes of graphs. 
nodes <- sapply(results,"[",3)
rand_nodes <- sapply(rand_results,"[",3)

# Function to calculate length, union, and intersection.
func <- function(x,y){
  len_x <- length(x)
  len_y <- length(y)
  u <- sum(x %in% y)
  i <- length(unique(c(x,y)))
  return(list(length_x = len_x, length_y = len_y, 
              union = u, intersection = i))
}
res <- mapply(func,nodes,rand_nodes)

# Collect every four... elements of list, and combine as dm. 
b <- split(res, rep(seq(from=1, to = 8, by = 1),each = 4))
dm <- do.call(rbind,lapply(b,function(x) unlist(x)))
rownames(dm) <- names(sigEntrez)
colnames(dm) <- c("len(V(subg))","len(V(rand))", "intersection","union")
df <- as.data.frame(dm)
df$percent <- round(100*(df$intersection/df$union),2)
df$seeds <- unlist(lapply(results,function(x) length(x[[1]])))

#-------------------------------------------------------------------------------
#' ## Examine pairwise overlap in DEP communities.
#-------------------------------------------------------------------------------

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
A <- list()
B <- list()
for (i in 1:dim(contrasts)[1]){
  a <- unlist(nodes[contrasts$ConditionA[i]])
  b <- unlist(nodes[contrasts$ConditionB[i]])
  A[[i]] <- length(a)
  B[[i]] <- length(b)
  int[[i]] <- intersect(a,b)
}

contrasts$Name <- paste(contrasts$ConditionA,
                        contrasts$ConditionB,sep="_U_")
names(int) <- contrasts$Name

# Add intersection to contrasts.
contrasts$Intersection <- lapply(int,function(x) length(x))
contrasts$A <- unlist(A)
contrasts$B <- unlist(B)
contrasts$C <- contrasts$A + contrasts$B

# Calculate percent intersection.
int <- list()
for (i in 1:dim(contrasts)[1]){
  a <- unlist(nodes[contrasts$ConditionA[i]])
  b <- unlist(nodes[contrasts$ConditionB[i]])
  int[[i]] <- length(intersect(a,b))/length(unique(c(a,b)))
}

# Add to contrasts dm
contrasts$Percent <- unlist(int)

# Make overlap matrix.
dm <- matrix(contrasts$Intersection,nrow=8,ncol=8)
rownames(dm) <- colnames(dm) <- row_names

# heirarchical clustering of dissimilarity matrix calculated as 1-percent_overlap.
dm <- matrix(contrasts$Percent,nrow=8,ncol=8)
rownames(dm) <- colnames(dm) <- row_names
diss <- 1 - dm
hc <- hclust(as.dist(diss), method = "average")
dendro <- ggdendrogram(hc, rotate=TRUE, labels = FALSE)
# Strip labels. 
dendro <- dendro +  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
dendro

# Remove upper tri and melt.
dm[lower.tri(dm)] <- NA
df <- melt(dm,na.rm = TRUE)

# Add intersection.
idx <- match(paste(df$Var1,df$Var2), paste(contrasts$ConditionA,contrasts$ConditionB))
df$intersection <- unlist(contrasts$Intersection[idx])

# Generate plot.
# Fix colors. 
# Fix rounding of percent overlap.
# Make theme consistent with other plot. 

# Order df based on dendrogram. 
levels(df$Var1) <- hc$labels[hc$order]
levels(df$Var2) <- hc$labels[hc$order]

plot <- ggplot(df, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "black") + 
  geom_text(aes(Var2, Var1, label = round(value,1)), color = "black", size = 1.0) +
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

plot <- plot + theme(legend.position = "none")
plot

# Save heatmap and dendrogram. 
file <- paste0(outputfigsdir,"/",outputMatName,"DEP_Community_Overlap_matirx.tiff")
ggsave(file,plot, width = 3, height = 3, units = "in", dpi = 300)

file <- paste0(outputfigsdir,"/",outputMatName,"DEP_Community_Overlap_dendro.tiff")
ggsave(file,dendro, width = 3, height = 3, units = "in", dpi = 300)

#-------------------------------------------------------------------------------
#' ## Generate custom gene lists.
#-------------------------------------------------------------------------------

file <- paste(datadir,"PPI Network" ,"MGIBatchReport_060619.xlsx", sep = "/")
mgi_go <- if(!exists("mgi_go")) read_excel(file,sheet=1)

# Keep Entrez genes.
mgi_go <- subset(mgi_go, mgi_go$`Input Type`=="Entrez Gene")

# Search for "ribosome" terms.
idx <- unlist(lapply(mgi_go$Term,function(x) grepl("ribosome|ribosomal",x)))
ribosome <- mgi_go[idx,]
length(unique(ribosome$Input))

# Search for proteasome.
idx <- unlist(lapply(mgi_go$Term,function(x) grepl("proteasome|proteolysis",x)))
proteasome <- mgi_go[idx,]
length(unique(proteasome$Input))

# Search for 
sum(community_results$Cortex.KO.Ube3a$nodes.entrez %in% proteasome$Input)

# Load node attribute table in cytoscape.
foo <- unique(proteasome$Input)
man <- as.data.frame(foo)
row.names(man) <- man$foo
man$category <- rep("proteasome",nrow(man))

loadTableData(man)

#-------------------------------------------------------------------------------
#' ## Build df of significantly dysregulated proteins.
#-------------------------------------------------------------------------------

# Build a df with statistical results.
#stats <- lapply(results,function(x) data.frame(Uniprot=x$Uniprot, FDR=x$FDR))
#names(stats) <- names(results)
#df <- stats %>% reduce(left_join, by = "Uniprot")
#colnames(df)[c(2:ncol(df))] <- names(stats)
#idx <- match(df$Uniprot,meta$uniprot)
#df <- add_column(df,Entrez = meta$entrez[idx], .after = 1)

stats <- lapply(results, function(x)
  data.frame(Uniprot = x$Uniprot, FDR = x$FDR))
names(stats) <- names(results)
df <- stats %>% reduce(left_join, by = "Uniprot")
colnames(df)[c(2:ncol(df))] <- names(stats)
head(df)

# Proteins with any significant change.
df$sigProt <- apply(df,1,function(x) any(as.numeric(x[c(2:9)])<0.05))

sum(df$sigProt) # Total
round(100*sum(df$sigProt)/nrow(df),3) # Percent

# Create table for cytoscape.
#file <- paste(outputtabsdir,"SigProts.csv",sep = "/")
#write.csv(df,file)

#-------------------------------------------------------------------------------
#' ## Generate subgraphs of WPCNA modules. 
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
#' ## Evaluate GO semantic similarity (~biological cohesiveness) of all modules. 
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
type <- 2 # MF, BP, CC
ontology <- c("MF","BP","CC")[type]
msGO <- list(msGOMF,msGOBP,msGOCC)[[type]]

# Evaluate GO similarity for all genes.
build_GO_similarity_network <- FALSE # otherwise, load from file. 

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
#' ## Break down DEP communities into functional subclusters....
#-------------------------------------------------------------------------------

# Create go similarity graph. 
gs <- graph_from_adjacency_matrix(goSim, mode = "undirected", weighted = TRUE)


subg <- community_results$Cortex.HET.Syngap1$subg


nodes <- community_results$Cortex.HET.Syngap1$nodes.entrez
n <- V(gs)[names(V(gs)) %in% nodes]

subg <- induced_subgraph(gs,n)


gl <- cluster_louvain(subg, E(gs))
gl$membership


man <- as.matrix(as_adjacency_matrix(gs, attr = "weight"))
foo <- as.matrix(as_adjacency_matrix(subg, attr = "weight"))

connectivity <- colSums(goSim)
# Examine fit. 
fit <- ggplotScaleFreePlot(connectivity)
fit

dim(foo)
foo[1:5,1:5]
max(foo)


#-------------------------------------------------------------------------------
#' ## Evaluate modularity of the GO similarity graph.
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
#' ## Test preservation of module's biological coherence with NetRep.
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
#' ## Build WGCNA module community sugraphs.
#-------------------------------------------------------------------------------

# Find nodes with 2+ connections to seed nodes.

# Create a Dictionary-like object of genes and entrez ids for mapping protein 
# identifiers.
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
#' ## Write DEP communities to file.
#-------------------------------------------------------------------------------

x <- sapply(results,"[",4)
names(x) <- gsub(".nodes.entrez","", names(x))
file <- "DEP_Communities.xlsx"
write.excel(x,file)

writeClipboard(as.matrix(meta$entrez))

#-------------------------------------------------------------------------------
#' ## Evaluate GO semantic similarity for DEP communities. 
#-------------------------------------------------------------------------------

# gene to entrez dictionary.
gene2entrez <- as.list(meta$entrez)
names(gene2entrez) <- meta$gene

# goSim is the GO semantic similarity matrix for all genes. 

# Empty list for output of loop.
out <- list()

# List of genes to pass to loop. 
#gene_list <- rand_community_nodes
gene_list <- sapply(results,"[",3)
class(gene_list)
length(gene_list)
#class(gene_list[[1]])
n <- length(gene_list)

# Loop through gene list, calculate GOSemSim:
for (i in 1:n){
  print(paste0("Working on subgraph: ", i, "..."))
  genes <- gene_list[[i]]
  #entrez <- genes
  entrez <- unlist(gene2entrez[genes])
  
  # Subset adjacency matrix. 
  idx <- colnames(goSim) %in% entrez
  subdm <- goSim[idx,idx]
  diag(subdm) <- NA
  
  # heirarchical clustering
  #diss <- 1 - subdm
  #hc <- hclust(as.dist(diss), method = "average")
  #dendro <- ggdendrogram(hc, rotate=FALSE, labels = FALSE)
  #dendro
  
  #mean(subdm,na.rm=TRUE)
  # Calculate average edge weight. 
  out[[i]] <- mean(subdm, na.rm = TRUE)
}
  
res <- do.call(rbind,out)
hist(res)


# Evaluate centrality as mean clustering coefficient in the GO similarity graph. 
  #subg <- lapply(goSim, function(x) 
  #  graph_from_adjacency_matrix(x, mode = "undirected", weighted = TRUE, diag = FALSE))
  #cc <- lapply(subg, function(x)
  #  transitivity(x, type = "average", isolates = "zero"))
  #out[[i]] <- cc


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
#' ## Examine GO enrichment for each DEP community.
#-------------------------------------------------------------------------------

## Prepare a matrix of class labels (colors) to pass to enrichmentAnalysis().
communityNodes <- lapply(results,function(x) x$nodes)

logic <- list()
for (i in 1:length(communityNodes)){
  logic[[i]] <- meta$gene %in% communityNodes[[i]]
}

labels <- do.call(cbind, logic)
rownames(labels) <- meta$entrez
colnames(labels) <- names(communityNodes)

# Convert TRUE to column names. 
logic <- labels == TRUE # 1 will become TRUE, and 0 will become FALSE.
# Loop through each column to replace 1 with column header.
for (i in 1:ncol(logic)){
  col_header <- colnames(labels)[i]
  labels[logic[,i],i] <- col_header
}

# Insure that labels is a matrix.
labels <- as.matrix(labels)
head(labels)

# look at the number of genes assigned to each cluster. 
table(labels)

# The labels matrix and vector of cooresponding entrez IDs 
# will be passed to enrichmentAnalysis().

# Build a GO annotation collection:
if (!exists(deparse(substitute(musGOcollection)))){
  musGOcollection <- buildGOcollection(organism = "mouse")}

# Creates some space by clearing some memory.
collectGarbage()

# Perform GO analysis for each module using hypergeometric (Fisher.test) test.
# As implmented by the WGCNA function enrichmentAnalysis().
# FDR is the BH adjusted p-value. 
# Insure that the correct background (used as reference for enrichment)
# has been selected!
# useBackgroud = "given" will use all given genes as reference background.

GOenrichment <- enrichmentAnalysis(
  classLabels = labels,
  identifiers = rownames(labels), # entrez ids
  refCollection = musGOcollection,
  useBackground = "given", # options are: given, reference (all), intersection, and provided. 
  threshold = 0.05,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "FALSE")

# Create some space by clearing some memory.
collectGarbage()

# Collect the results. 
results_GOenrichment <- list()
for (i in 1:length(GOenrichment$setResults)){
  results_GOenrichment[[i]] <- GOenrichment$setResults[[i]]$enrichmentTable
}
length(results_GOenrichment)
names(results_GOenrichment) <- colnames(labels)

# Write results to file. 
file <- paste0(outputtabsdir,"/",tissue,"_WGCNA_Analysis_DEP_Community_GOenrichment_Results.xlsx")
write.excel(results_GOenrichment,file)

## Visualize GO terms with GOscatter plots. 
# Module name should be a color for the plot.
colorvec <- c("green1","yellow1","blue1","purple1",
              "green2","yellow2","blue2","purple2")
godata <- results_GOenrichment
names(godata) <- colorvec

# Specify the top percent of terms to print with topN.
plots <- list()
for (i in 1:length(godata)){
  color <- colorvec[i]
  plots[[i]] <- ggplotGOscatter(godata, color = color, topN = 0.05)
}
names(plots) <- names(results_GOenrichment)

#fixme: add specific colors, change dot alpha may improve appearance. 
# ugg shank2...
plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]
plots[[5]]
plots[[6]]
plots[[7]]
plots[[8]]

# Loop to save plots.
for (i in 1:length(plots)){
  plot <- plots[[i]] + theme(legend.position = "none")
  namen <- gsub("\\.","_",names(plots)[i])
  file <- paste0(outputfigsdir,"/",outputMatName,namen,"DEP_Community_GO_scatter.tiff")
  ggsave(file, plot, width = 3, height = 3, units = "in", dpi = 300)
}

#-------------------------------------------------------------------------------
#' ## Evaluate go term similarity of DEP communities.
#-------------------------------------------------------------------------------

# Entrez genes for each DEP community (cluster).
clusters <- sapply(results,"[",4)
names(clusters) <- gsub(".nodes.entrez","",names(clusters))

# Evaluate GO semanitic similarity between DEP communities (clusters).
gosim <- mclusterSim(clusters, semData=msGOBP, measure="Wang", combine="BMA")

# hclustering
diag(gosim) <- NA
diss <- 1 - gosim
hc <- hclust(as.dist(diss), method = "average")
dendro <- ggdendrogram(hc)
dendro


#-------------------------------------------------------------------------------
#' ## Cluster GO terms for each DEP community by semantic similarity. 
#-------------------------------------------------------------------------------

# Build GO database.
if (!exists("msGOMF")){ msGOMF <- godata('org.Mm.eg.db', ont= "MF")}
if (!exists("msGOBP")){ msGOBP <- godata('org.Mm.eg.db', ont= "BP")}
if (!exists("msGOCC")){ msGOCC <- godata('org.Mm.eg.db', ont= "CC")}

msGO <- list(msGOMF,msGOBP,msGOCC)
names(msGO) <- c("MF","BP","CC")

# Collect genes associated with each DEP community.
#genes <- unlist(sapply(results,"[",3))
#clusters <- list()

#for (i in 1:length(results)){
#  genes <- unlist(results[[i]][3])
#  entrez <- unlist(gene2entrez[genes])
#  clusters[[i]] <- entrez
#}
#names(clusters) <- names(results)

# Build collection of GO terms. 
#x <- org.Mm.egGO
#mapped_genes <- mappedkeys(x)
#musGO <- as.list(x[mapped_genes])

# Loop through clusters and get GO terms associated with this group of proteins.
#out <- list()
#for (i in 1:length(clusters)){
#  idx <- clusters[[i]]
#  out[[i]] <- unlist(lapply(musGO[idx],function(x) names(x)))
#}

goIDs <- lapply(results_GOenrichment,function(x) x$dataSetID)
namen <- lapply(results_GOenrichment,function(x) x$inGroups)

for (i in 1:length(goIDs)){
  govec <- goIDs[[i]]
  names(govec) <- gsub("GO|GO.","",namen[[i]])
  goIDs[[i]] <- govec
}

# For each group of go terms, evaluate similarity.
out <- list()
for (i in c(1,2,3,4,5,7,8)){
  go1 <- goIDs[[i]]
  go1 <- subset(go1, names(go1) == "BP")
  gosim <- mgoSim(go1, go1, semData=msGO$BP, measure="Wang", combine=NULL)
  
  # hclustering.
  diss <- 1 - gosim
  hc <- hclust(as.dist(diss), method = "average")
  dendro <- ggdendrogram(hc, rotate=FALSE, labels = TRUE)
  out[[i]] <- list(dendrogram = dendro, hclust = hc, goSimilarity = gosim)
}

names(out) <- names(goIDs)
lapply(out,function(x) x$dendrogram)

out$Cortex.HET.Syngap1$dendrogram
hc <- out$Cortex.HET.Syngap1$hclust
v <- cutree(hc,k=4)
length(unique(v))
results_GOenrichment$Cortex.HET.Syngap1$dataSetID %in% names(v)
out$Cortex.HET.Syngap1$dendrogram


# Get GO terms associated with DEP community.
# Evaluate semantic similarity amongst these terms.
# Cluster and display on PPI graph. 
#mgoSim(go1, go2, semData=hsGO, measure="Wang", combine=NULL)

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
``
