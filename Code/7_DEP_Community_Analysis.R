#' ---
#' title: WPCNA analysis of DEP communities. 
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
CodeVersion <- "DEP_Community_Analysis"

# Define tisue type: cortex = 1; striatum = 2.
type <- 3
tissue <- c("Cortex", "Striatum", "Combined")[type]

# Set the working directory.
rootdir <- "D:/Documents/R/Synaptopathy-Proteomics"
setwd(rootdir)

# Set any other directories.
functiondir <- paste(rootdir, "Code", sep = "/")
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

# Load DEP communities. 
file <- paste0(Rdatadir,"/","DEP_Communities.Rds")
DEP_communities <- readRDS(file)

#-------------------------------------------------------------------------------
#' ## Adding KNN to DEP communities.
#-------------------------------------------------------------------------------

# Create a dictionary like object for mapping entrez to gene|uniprot.
entrez2protein <- as.list(meta$protein)
names(entrez2protein) <- meta$entrez

# Calculate protein co-expression (correlation) matrix. 
r <- bicor(t(cleanDat))
diag(r) <- NA
r[lower.tri(r)] <- NA

# Loop to get K(3)NN for each seed node. 
out <- list()
for (i in 1:length(DEP_communities)){
  
  # Collect the data
  community <- names(DEP_communities)[i]
  community_data <- DEP_communities[[community]]

  s <- community_data$seeds
  m <- community_data$nodes.entrez
  nodes <- unlist(entrez2protein[s])
  c_nodes <- unlist(entrez2protein[m])
  subg <- community_data$subg

  # Define a function to get KNN neighbors.
  fx <- function(x, n = 3){ 
    x.sorted <- sort(x, decreasing = TRUE)
    x.knn <- x.sorted[c(1:n)]
    return(as.list(x.knn))
    }
  
  # Get KNN for all seed proteins.  
  data_knn <- apply(r,2,fx)
  data_knn <- data_knn[nodes]
  df <- melt(do.call(rbind, lapply(data_knn,function(x) unlist(x))))
  colnames(df) <- c("NodeA","NodeB","Bicor")
  
  # Define proteins of interest as seeds + community nodes + knn.
  prots <- unique(c(df$NodeA,df$NodeB,nodes,c_nodes))
  
  df.summary <- data.frame(
    seeds = length(s),
    community = length(m),
    combined = length(prots))
  out[[i]] <- list(proteins = prots, summary = df.summary)
  }

names(out) <- names(DEP_communities)

# Check numbers.
x <- sapply(out,"[",2)
do.call(rbind,x)

# Write to file.
file <- paste0(Rdatadir,"/","DEP_KNN_Communities.Rds")
saveRDS(out,file)

#-------------------------------------------------------------------------------
#' ## Start WGCNA. Choosing a soft thresholding power, Beta.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Estimate powers?
estimatePower <- TRUE

# Load TAMPOR cleanDat from file: #3022 rows.
datafile <- paste(Rdatadir,tissue,"TAMPOR_data_outliersRemoved.Rds",sep="/")
cleanDat <- readRDS(datafile)
cleanDat <- log2(cleanDat)
cleanDat[1:5,1:5]
dim(cleanDat) # 1 outlier removed.

# Subset cleanDat based on prots of interest. 
prots <- out$Cortex.HET.Syngap1$proteins
subg_name <- "Cortex.HET.Syngap1"
subDat <- subset(cleanDat, rownames(cleanDat) %in% prots)
dim(subDat)

# Strip module, metamodule and modulColor attributes from subg.
attr <- names(vertex_attr(subg, index = V(subg)))[c(4,5,6)]
for (i in 1:length(attr)){
  subg <- delete_vertex_attr(subg, attr[i])
}

# Check.
length(names(vertex_attr(subg, index = V(subg)))) # 35

# Load combined sample info.
traitsfile <- paste(Rdatadir,tissue,"Combined_Cortex_Striatum_traits.Rds",sep="/")
sample_info <- readRDS(traitsfile)
sample_info[1:5,1:5]
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
if (estimatePower==TRUE){
  sft <- pickSoftThreshold(t(subDat),
                           powerVector = powers, corFnc = "bicor",
                           blockSize = 15000, verbose = 3, networkType = "signed")
  
  # Create table. 
  mytable <- round(sft$fitIndices,2)
  mytable$truncated.R.sq <- NULL
  table <- tableGrob(mytable, rows = NULL)
  grid.arrange(table)
  
  # Save table as tiff.
  #file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology_Table.tiff")
  #ggsave(file,table, width = 3, height = 2.5, units = "in")
  
  # Figure. ggplotScaleFreeFit() generates three plots.
  plots <- ggplotScaleFreeFit(sft)
  plots$Grid
  
  # Save as tiff.
  #file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology_1.tiff")
  #ggsave(file,plots$ScaleFreeFit, width = 3, height = 2.5, units = "in")
  
  # Save as tiff.
  #file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology_2.tiff")
  #ggsave(file,plots$MeanConnectivity, width = 3, height = 2.5, units = "in")
  
  # Save plots and table as PDF.
  #plot_list <- list(table,plots$ScaleFreeFit, plots$MeanConnectivity)
  #file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology.pdf")
  #ggsavePDF(plot_list,file)
}

# Choose a scale free power.
power <- 7 # Why does the plot look weird? What does this mean? 
  
# Calculate node connectivity in the Weighted co-expression network.  
connectivity <- softConnectivity(
  datExpr=t(subDat), 
  corFnc = "bicor",
  weights = NULL,
  type = "signed",
  power = power, 
  blockSize = 15000, 
  minNSamples = NULL, 
  verbose = 0, 
  indent = 0)

# Examine fit. 
fit <- ggplotScaleFreePlot(connectivity)
fit$stats
plot <- fit$ggplot
plot

# Save fig. 
#file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeFit",".tiff")
#ggsave(file,plot, width = 3, height = 2.5, units = "in")

#-------------------------------------------------------------------------------
#' ## Build the WPCNA Network.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Load blockwiseModules Parameters.
files <- list.files(Rdatadir,pattern = "Stats")
file <- paste0(Rdatadir,"/",files[1])
sampled_params <- do.call(rbind,sapply(readRDS(file),"[",1))
rownames(sampled_params) <- paste0("params_",(1:nrow(sampled_params)))
sampled_params$iter <- c(1:nrow(sampled_params))
dim(sampled_params)[1]

# Choose network building parameters.
params_iter <- 738 #573 #630 #481 #223 #216 #981 #981 # 10 #691 10, 507 530
params <- sampled_params[params_iter,]
params[,c(1:7)]

# Network Parameters table.
mytable <- params[,c(1:7)]
table <- tableGrob(mytable, rows = NULL)
#file <- paste0(outputfigsdir,"/",outputMatName,"Network_Parameters.tiff")
#ggsave(file,table,width = 9, height = 1)

# Network parameters:
power <- power
networkType <- "signed"
corType <- "bicor"
enforceMMS <- TRUE         # Should minimal modules size be inforced? Done after network building.

# Other key parameters for optimization:
minModSize <- params$minModSize          # Minimum module size. seq(1,50,by=1) 
deepSplit <- params$deepSplit            # Sensitivity for module splitting [0-4]. 4 is most sensitive. Increasing results in more modules.  
mergeCutHeight <- params$mergeCutHeight  # Cut height for module detection. Was 0.07. Increasing results in more modules.  
reassignThresh <- params$reassignThresh  # pvalue threshold for reassigning nodes to modules. 
minKMEtoStay <- params$minKMEtoStay      # minimum module connectivity score for assigning to a module. 
minCoreKMESize <- params$minCoreKMESize  # minimim number of genes in a modules with minKMEtoStay.
pamStage <- params$pamStage

# Insure that minCoreKME is not less than minModSize.
if (minCoreKMESize < minModSize){
  minCoreKMESize <- minModSize
  print(paste0("Warning: minCoreMKESize set to minModsize"," (",minModSize,")"))
}

# Defaults:
maxBlockSize <- 12000                    # maximum block size for module detection. 
detectCutHeight = 0.995                  # dendrogram cut height for module detection. 

# Call blockwiseModules to build WGCNA network. 
# Setting saveTOM = FALSE will really slow thing down.
net <- blockwiseModules(t(subDat),
                        power = power, 
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
                        maxBlockSize = maxBlockSize)

# Check the number of modules. 
nModules_original <- length(unique(net$colors))
nModules_original 

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
adjm <- ((1+r)/2)^power # Signed network. 

data_list <- list(data = t(subDat))  # The protein expression data. 
correlation_list <- list(data = r)     # The bicor correlation matrix. 
network_list <- list(data = adjm)      # The weighted, signed co-expresion network.   
module_labels <- net$colors            # Module labels. 
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
  nThreads-1, 
  #nPerm = 100000, # Increase nPerm to 100,000 in order to stabilize the result with large number of modules. 
  null = "overlap", 
  alternative = "greater", 
  simplify = TRUE,
  verbose = TRUE)

# Collect stats. 
preservation <- preservation[c("observed","p.values")]

# Get the maximum permutation test p-value.
maxPval <- apply(preservation$p.values, 1, function(x) max(x,na.rm=TRUE))

# Modules removed if adjusted pvalue is greater than alpha = 0.05.
alpha = 0.05
modules_out <- names(maxPval)[maxPval>alpha/nModules_original]
nModules_out <- length(modules_out)
names(net$MEs)[grepl(paste(modules_out,collapse="|"),names(net$MEs))]<-"MEgrey"

# Set non-significant modules to grey.
net$colors[net$colors %in% modules_out] <- "grey"

# Total number of modules, excluding grey:
nModules <- length(unique(net$colors))-1
params$nModules <- nModules
nModules

# Check percent grey.
percent_grey <- round(100*sum(net$colors=="grey")/length(net$colors),2)
print(paste("Percent grey nodes =", percent_grey))
params$PercentGrayNodes <- percent_grey

# Table of key network stats.
mytable <- data.frame(
  nNodes = sum(net$colors!="grey"),
  PercentGrey = round(params$PercentGrayNodes,2),
  nModules = nModules,
  MedianCoherence = round(params$medianModCoherence,3))
#NetworkModularity = round(params$q2,3))
table <- tableGrob(mytable, rows = NULL)
grid.arrange(table)
#file <- paste0(outputfigsdir,"/",outputMatName,"Key_Network_Stats.tiff")
#ggsave(file,table,width = 6, height = 1)

#-------------------------------------------------------------------------------
#' ## Enforce min module size. Recalculate MEs.
#-------------------------------------------------------------------------------

# Fraction of un-assigned proteins. 
print(paste("Percent grey = ",round(100*table(net$colors)["grey"]/length(net$goodGenes),2),
            " (n=",table(net$colors)["grey"],")",sep=""))

# Module Summary (excluding grey)
nModules <- length(table(net$colors)) - 1
print(paste("Total number of modules =", nModules,"(without grey)"))
modules <- cbind(colnames(as.matrix(table(net$colors))), table(net$colors))
orderedModules <- cbind(Mnum = paste("M", seq(1:nModules), sep = ""), Color = labels2colors(c(1:nModules)))
modules <- modules[match(as.character(orderedModules[, 2]), rownames(modules)), ]
moduleColors <- as.data.frame(cbind(orderedModules, Size = modules))

# If necessary, enforce minModSize. Set these to grey.
# Loop to enforce minModSize:
if (enforceMMS) {
  removedModules <- orderedModules[which(modules < minModSize), "Color"]
  print(paste("There are", length(removedModules), "modules that contain less than",
              minModSize,"proteins."))
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
                        softPower = power,
                        scale = TRUE,
                        verbose = 0, indent = 0)
net$MEs <- MEs$eigengenes

# module assignments.
moduleColors <- data.frame(protein = rownames(subDat),
                           module  = net$colors)

# Check number of modules.
nModules_original
nModules
nModules_out

#-------------------------------------------------------------------------------
#' Send DEP Community with Co-expression modules to Cytoscape.
#-------------------------------------------------------------------------------

# Send graph to cytoscape?
send_to_cytoscape = TRUE

# Build df of node attributes. 
# FIX THIS!!! 
idx <- match(unlist(entrez2protein[names(V(subg))]), rownames(subDat))
hex <- lapply(as.list(net$colors[idx]), function(x) col2hex(x))
df <- data.frame(Protein = unlist(entrez2protein[names(V(subg))]),
                 Entrez = (V(subg)),
                 Color = net$colors[idx],
                 HexColor = unlist(hex))
rownames(df) <- df$Entrez

# Send to cytoscape with RCy3!
if (send_to_cytoscape == TRUE){
  cytoscapePing()
  quiet(RCy3::createNetworkFromIgraph(subg,subg_name))
  
  # Load node attribute table in cytoscape.
  loadTableData(df)
  
  # Create custom syle to customize appearance. 
  geno <- strsplit(subg_name, "\\.")[[1]][3]
  style.name <- subg_name
  colvec <- df$HexColor
  
  defaults <- list(NODE_FILL = col2hex("grey"),
                   NODE_SHAPE="Ellipse",
                   NODE_SIZE=55,
                   EDGE_WIDTH = 2.0,
                   EDGE_TRANSPARENCY=120)
  
  nodeLabels <- mapVisualProperty('node label','Symbol','p')
  nodeFills <- mapVisualProperty('node fill color','HexColor','p', colvec)
  
  #edgeWidth <- mapVisualProperty('edge width','weight','p')
  createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills))
  lockNodeDimensions(TRUE, style.name)
  setVisualStyle(style.name)
  
  # Apply perfuse force directed layout. 
  layoutNetwork(layout.name = "force-directed")
  setNodeColorDefault(col2hex("grey"), style.name)
  
  # Apply node color mapping to score. (this will overwrite other mapping)                    
  #column <- paste("score",names(sigProts)[i])
  #control.points <- c(min(df[column]), 0.0, max(df[column]))
  #cols <- c(col2hex("blue"), col2hex("white"), col2hex("red"))
  #setNodeColorMapping(column, control.points, cols, style.name = style.name)
}

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
                                         numericMeta$Sample.Model, sep=".")

# Pool WT within a tissue.
idx <- grepl("Cortex.WT",numericMeta$Tissue.Sample.Model)
numericMeta$Tissue.Sample.Model[idx] <- "Cortex.WT"
idy <- grepl("Striatum.WT",numericMeta$Tissue.Sample.Model)
numericMeta$Tissue.Sample.Model[idy] <- "Striatum.WT"
unique(numericMeta$Tissue.Sample.Model)

# Pull out as many numeric traits for correlation later as we can.
numericMeta$Group <- tempVec <- as.vector(as.data.frame(do.call(rbind, 
                                                                strsplit(colnames(cleanDat), "\\.")))[, 1])
numericMeta$Group <- as.numeric(as.factor(numericMeta$Group))

# Tissue
numericMeta$Tissue <- as.numeric(as.factor(numericMeta$Tissue))-1

# Genotype groupings (genetic background). 
numericMeta$Model <- gsub(" ","",numericMeta$Model)
numericMeta$Syngap1 <- as.numeric(numericMeta$Model=="Syngap1")
numericMeta$Ube3a <- as.numeric(numericMeta$Model=="Ube3a")
numericMeta$Shank2 <- as.numeric(numericMeta$Model=="Shank2")
numericMeta$Shank3 <- as.numeric(numericMeta$Model=="Shank3")

# Sex
numericMeta$Sex <- as.numeric(as.factor(numericMeta$Sex))-1

# Make Syngap1 WT v KO column.
temp_df <- as.data.frame(do.call(rbind,strsplit(numericMeta$Sample.Model, "\\.")))
numericMeta$Syngap1_KO <- NA
numericMeta$Syngap1_KO[which(temp_df$V1=="WT")] <- 0
numericMeta$Syngap1_KO[which(temp_df$V1=="HET" & temp_df$V2=="Syngap1")] <- 1

# Make Ube3a WT v HET column.
temp_df <- as.data.frame(do.call(rbind,strsplit(numericMeta$Sample.Model, "\\.")))
numericMeta$Ube3a_KO <- NA
numericMeta$Ube3a_KO[which(temp_df$V1=="WT")] <- 0
numericMeta$Ube3a_KO[which(temp_df$V1=="KO" & temp_df$V2=="Ube3a")] <- 1

# Make Shank2 WT v KO column.
temp_df <- as.data.frame(do.call(rbind,strsplit(numericMeta$Sample.Model, "\\.")))
numericMeta$Shank2_KO <- NA
numericMeta$Shank2_KO[which(temp_df$V1=="WT")] <- 0
numericMeta$Shank2_KO[which(temp_df$V1=="KO" & temp_df$V2=="Shank2")] <- 1

# Make Shank3 WT v KO column.
temp_df <- as.data.frame(do.call(rbind,strsplit(numericMeta$Sample.Model, "\\.")))
numericMeta$Shank3_KO <- NA
numericMeta$Shank3_KO[which(temp_df$V1=="WT")] <- 0
numericMeta$Shank3_KO[which(temp_df$V1=="KO" & temp_df$V2=="Shank3")] <- 1

# Control vs. disease.
numericMeta$SampleType <- gsub(" ","", numericMeta$SampleType)
numericMeta$ASD <- NA
numericMeta$ASD[numericMeta$SampleType=="WT"] <- 0
numericMeta$ASD[numericMeta$SampleType=="KO"] <- 1
numericMeta$ASD[numericMeta$SampleType=="HET"] <- 1

# Tissue:Genotype groupings (genetic background). 
f1 <- as.factor(numericMeta$TissueType)
f2 <- as.factor(numericMeta$Model)
mod <- model.matrix(~0+f1:f2)
colnames(mod) <- apply(expand.grid(levels(f1), levels(f2)), 1, paste, collapse=".")
numericMeta <- cbind(numericMeta,mod)

# Tissue:Sex
f1 <- as.factor(numericMeta$TissueType)
f2 <- as.factor(numericMeta$SexType)
mod <- model.matrix(~0+f1:f2)
colnames(mod) <- apply(expand.grid(levels(f1), levels(f2)), 1, paste, collapse=".")
numericMeta <- cbind(numericMeta,mod)

# Generate Tissue Specific Sample.Model contrasts
g <- paste(numericMeta$TissueType,numericMeta$Sample.Model,sep=".")
allContrasts <- combn(unique(g),2)
# Keep if tissue and model are equal.
keepers <- function(x){
  y <- do.call(rbind,strsplit(x,"\\."))
  logic <- y[1,1] == y[2,1] & y[1,3] == y[2,3]
  return(logic)}
keep <- apply(allContrasts,2,function(x) keepers(x))
contrasts <- allContrasts[,keep]
# Pool WT within a tissue.
contrasts[grepl("Cortex.WT",contrasts)] <- "Cortex.WT"
contrasts[grepl("Striatum.WT",contrasts)] <- "Striatum.WT"
# Coerce to list.
contrasts_list <- apply(contrasts,2,as.list)
# lapply through contrasts list and generate model vector.
mod_list <- lapply(contrasts_list,function(x) match(numericMeta$Tissue.Sample.Model,x)-1)
mod <- do.call(cbind,mod_list)
# Add names
colnames(mod) <- sapply(contrasts_list,"[",1)
# Fix Cortex.KO.Shank3 column name.
colnames(mod)[grep("Cortex.WT",colnames(mod))] <- "Cortex.KO.Shank3"

# Add to numericMeta
numericMeta <- cbind(numericMeta,mod)

# Tissue specific disease status (Control v ASD).
g <- numericMeta$SampleType
g[grepl("KO|HET",g)] <- "ASD"
g[grepl("WT",g)] <- "Control"
f1 <- paste(numericMeta$TissueType,g,sep=".")
allContrasts <- list(c("Cortex.ASD","Cortex.Control"),
                     c("Striatum.ASD","Striatum.Control"))
mod_list <- lapply(allContrasts,function(x) match(f1,x)-1)
mod <- do.call(cbind,mod_list)
colnames(mod) <- colnames(mod) <- sapply(allContrasts,"[",1)
numericMeta <- cbind(numericMeta,mod)

# Sex specific disease status (control v ASD).
g <- numericMeta$SampleType
g[grepl("KO|HET",g)] <- "ASD"
g[grepl("WT",g)] <- "Control"
numericMeta$M.ASD <- as.numeric(numericMeta$SexType == "M" & g == "ASD")
numericMeta$F.ASD <- as.numeric(numericMeta$SexType == "F" & g == "ASD")

# Determine numerical indices. The columns of numericMeta with numerical data.
# Warnings OK; This determines which traits are numeric and if forced to numeric values, 
# non-NA values do not sum to 0.
numericIndices <- unique(c(which(!is.na(apply(numericMeta, 2, function(x) sum(as.numeric(x))))), 
                           which(!(apply(numericMeta, 2, function(x) sum(as.numeric(x), na.rm = T))) == 0)))

#-------------------------------------------------------------------------------
#' ## VerboseBoxplots - vizualize module expression across traits.
#-------------------------------------------------------------------------------

# Should plots be saved?
saveplots = FALSE
savegroups = FALSE

# Calculate Module EigenProteins.
MEList <- moduleEigengenes(t(subDat), colors = net$colors)
MEs <- orderMEs(MEList$eigengenes)
colnames(MEs) <- gsub("ME", "", colnames(MEs)) 
rownames(MEs) <- rownames(numericMeta)

# Insure traits are in matching order.
traits <- sample_info
idx <- match(rownames(MEs),rownames(traits))
traits <- traits[idx,]
all(rownames(traits)==rownames(MEs))

# Define groups, the biological groups of interest. 
groups <- paste(numericMeta$TissueType,numericMeta$Sample.Model,sep=".")
groups[grepl("Cortex.WT",groups)] <- "WT.Cortex"
groups[grepl("Striatum.WT",groups)] <- "WT.Striatum"
unique(groups)

# Calculate Kruskal-Wallis pvalues for all modules (columns of MEs df).
KWtest <- apply(MEs,2,function(x) kruskal.test(x,as.factor(groups)))

# Extract pvalues from the list of KW tests. 
# The pvalue is the 3rd element of each list.
KW_results <- as.data.frame(do.call(rbind, sapply(KWtest,"[",3)))

# Clean up the result. Drop grey!
colnames(KW_results) <- "p.value"
rownames(KW_results) <- gsub(".p.value","",rownames(KW_results))
KW_results <- add_column(KW_results,Module = rownames(KW_results),.before = 1)
rownames(KW_results) <- NULL
KW_results <- KW_results[!KW_results$Module=="grey",]

## Which comparisons are significant? 
# Bonferroni correction for nModules.
KW_results$p.adj <- p.adjust(as.numeric(KW_results$p.value),method = "bonferroni")
# Benjamini-Hochberg FDR correction. 
KW_results$FDR <- p.adjust(as.numeric(KW_results$p.value), method = "BH")
KW_results <- KW_results[order(KW_results$p.value),]

# nsig modules. 
KWsigModules <- KW_results$Module[KW_results$FDR<0.05]
print(paste("nModules with KW FDR < 0.05 =", length(KWsigModules)))

# Split Module EigenProtein (MEs) dm into a list of column vectors for lapply. 
ME_list <- split(as.matrix(MEs), rep(1:ncol(MEs), each = nrow(MEs)))
names(ME_list) <- colnames(MEs)

# Add vector of groups
ME_list <- lapply(ME_list,function(x) data.frame(x = x, groups = groups))

# Define levels for order of bars in plot.
#levels <- c("M.WT","M.ASD","F.WT","F.ASD")
levels <- c("WT.Cortex","WT.Striatum",
            "Cortex.KO.Shank2","Striatum.KO.Shank2",
            "Cortex.KO.Shank3","Striatum.KO.Shank3",
            "Cortex.HET.Syngap1","Striatum.HET.Syngap1",
            "Cortex.KO.Ube3a","Striatum.KO.Ube3a")

# Generate contrasts matrix for comparisons of interest.
contrasts <- makePairwiseContrasts(list("M.WT","F.WT"),list("M.ASD","F.ASD"))
g1 <- list("WT.Cortex","WT.Striatum")
g2 <- list(
  c("Cortex.KO.Shank2","Cortex.KO.Shank3","Cortex.HET.Syngap1","Cortex.KO.Ube3a"),
  c("Striatum.KO.Shank2","Striatum.KO.Shank3","Striatum.HET.Syngap1","Striatum.KO.Ube3a"))
contrasts <- makePairwiseContrasts(g1,g2)

# Loop through ME_list and generate verboseBoxPlot.
# lapply wont work here because the name is not preserved when you call lapply()...
# method is the p.adj method for the Dunn's test p-value.
plot_data <- list()
for (i in 1:dim(MEs)[2]){
  x <- ME_list[[i]]$x
  g <- ME_list[[i]]$groups
  color <- names(ME_list)[[i]]
  plot <- ggplotVerboseBoxplot(x,g,
                               levels,
                               contrasts,
                               color,
                               stats=TRUE,
                               method="dunn", 
                               correction_factor = 1)
  plot_data[[i]] <- plot
  names(plot_data)[[i]] <- color
}

# Extract plots. 
plot_list <- sapply(plot_data,"[",1)

# Store as VBplots
vbplots <- plot_list

# Extract post-hoc test stats.
Dtest_stats <- sapply(plot_data,"[",3)
names(Dtest_stats) <- sapply(strsplit(names(Dtest_stats),"\\."),"[",1)

# Loop to add module column.
for (i in 1:length(Dtest_stats)){
  df <- Dtest_stats[[i]]
  namen <- names(Dtest_stats)[i]
  df <- add_column(df,Module = namen,.before=1)
  Dtest_stats[[i]] <- df
}

# Number of significant post-hoc tests.
nsig <- do.call(rbind,lapply(Dtest_stats,function(x) sum(x$P.adj<0.05)))
idx <- match(KW_results$Module,rownames(nsig))
KW_results$DunnettTestNsig <- nsig[idx,]

# Convert Dtest stats into data frame. By casting long data into wide.
# dcast should work for multiple value.var, but it isnt...
df <- do.call(rbind,
              lapply(Dtest_stats, 
                     function(x) reshape2::dcast(x, Module ~ Comparison, value.var=c("Z"))))
df2 <- do.call(rbind,
               lapply(Dtest_stats, 
                      function(x) reshape2::dcast(x, Module ~ Comparison, value.var=c("P.unadj"))))

# Fix column names. 
colnames(df)[2:ncol(df)] <- paste(colnames(df)[2:ncol(df)],"Z")
colnames(df2)[2:ncol(df2)] <- paste(colnames(df2)[2:ncol(df2)],"P.value")

# Bind as single df. Drop grey.
data <- cbind(df,df2)
data <- data[!data$Module=="grey",]
dim(data)

# Add KW P.value and FDR.
idx <- match(data$Module,KW_results$Module)
data <- add_column(data,KW.P.value = KW_results$p.value[idx],.after = 1)
data <- add_column(data,KW.FDR = KW_results$FDR[idx],.after = 2)

# Subset only sig... 
data <- subset(data,data$KW.P.value<0.05)
data <- subset(data,data$KW.FDR<0.05)

# Count instances of significant change.
idx <- grep("P.value",colnames(data))[-1] # Remove kw p.value
logic <- data[,idx] < 0.05
ids <- sapply(strsplit(colnames(logic),"\\ "),"[",3)
out <- list()
for (i in 1:ncol(logic)){
  col <- logic[,i]
  col[col] <- ids[i]
  out[[i]] <- col
}
dm <- do.call(cbind,out)
dm[dm=="FALSE"] <- ""
colnames(dm) <- ids
data$SigTests <- apply(dm,1,function(x) paste(unique(x),collapse=" "))
KWdf <- data

# Subset kwdata
KW_sub <- subset(KW_results,KW_results$FDR<0.05 & KW_results$DunnettTestNsig>0)

# Summarize modules with convergance...
df <- count(KWdf$SigTests)
df <- df[-1,]
foo <- split(KWdf,KWdf$SigTests)
man <- lapply(foo,function(x) unique(x$Module))
module_overlap <- man

# Save this to file.
#file <- paste(Rdatadir,"module_overlap.Rds", sep="/")
#saveRDS(module_overlap,file)

# Write groups of tiffs to file.
# Fix names of vbplots.
names(vbplots) <- sapply(strsplit(names(vbplots),"\\."),"[",1)

# Make a directory for storing all vbplots.
#new_dir <- paste(outputfigsdir,"VBplot_Groups",sep="/")
#dir.create(new_dir)

# Loop over groups of plots, and save them in common directory.
#if (savegroups == TRUE){
#  print("Saving plots!")
#  for (i in 1:length(module_overlap)){
#    namen <- sub(" ","",names(module_overlap)[i])
#    namen <- sub("\\ ","_",gsub("\\.","",namen))
#    group <- module_overlap[[i]]
#    output_directory <- paste(new_dir,namen,sep="/")
#    dir.create(output_directory)
#    plot_names <- module_overlap[[i]]
#    # Loop to save plots as individual tiffs.
#    for (k in 1:length(plot_names)){
#      plot_name <- plot_names[k]
#      plot <- vbplots[[plot_name]]
#      file <- paste0(output_directory,"/",plot_name,".tiff")
#      ggsave(file,plot,width = 8.17, height = 5.17)
#    }
#  }
#}

# Summarize number of post-hoc changes in significant modules...
cols <- c(1,grep(" P.value",colnames(data)))
df <- melt(data[,cols], id = "Module")
colnames(df) <- c("Module","Genotype","P.Value")
sig_counts <- df %>% group_by(Genotype) %>% 
  dplyr::summarise(Count = sum(P.Value<0.05))
sig_counts <- as.data.frame(sig_counts)
sig_counts$Genotype <- as.character(sig_counts$Genotype)
sig_counts$Genotype <- sapply(strsplit(sig_counts$Genotype," - "),"[",2)
sig_counts$Genotype <- gsub(" P.value","", sig_counts$Genotype)
knitr::kable(sig_counts)

# Save top sig plots...
#for (i in 1:dim(KW_sub)[1]){
#  file <- paste0(outputfigsdir,"/",outputMatName, "_", KW_sub$Module[i],"_verboseBoxplot",".tiff")
#  plot <- vbplots[[KW_sub$Module[i]]]
#  ggsave(file, plot, width = 3.25, height = 2.5, units = "in")
#}

# Edit plots such that x-axis labels are red if post-hoc test is significant.
for (i in 1:dim(KW_sub)[1]){
  plot <- vbplots[[KW_sub$Module[i]]]
  stats <- Dtest_stats[[KW_sub$Module[i]]]
  b <- ggplot_build(plot)
  foo <- b$data[[3]]
  sig <- foo$label[order(foo$x)]
  a <- rep("black",8)
  a[grepl("\\*",sig)] <- "red"
  a <- c(rep("black",2),a)
  plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
  file <- paste0(outputfigsdir,"/",outputMatName, "_", KW_sub$Module[i],"_verboseBoxplot",".tiff")
  ggsave(file, plot, width = 3.25, height = 2.5, units = "in")
}

