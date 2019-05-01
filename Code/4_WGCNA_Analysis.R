#' ---
#' title: WGCNA Analysis of TMT data. 
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
CodeVersion <- "Semi_Final"

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
outputMatName <- paste(tissue, "_WGCNA_Analysis_", sep = "")

# Globally set ggplots theme.
ggplot2::theme_set(theme_gray())

#-------------------------------------------------------------------------------
#' ## Load the EdgeR statistical results.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Fixme:
# InterBatch statistical comparisons with EdgeR GLM:
#file <- paste0(rootdir,"/","Tables/Combined/v14_TAMPOR/Combined_TMT_Analysis_TAMPOR_GLM_Results.xlsx")
#results <- lapply(as.list(c(1:8)),function(x) read_excel(file,x))
#names(results) <- excel_sheets(file)

#-------------------------------------------------------------------------------
#' ## Start WGCNA. Choosing a soft thresholding power, Beta.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Estimate powers?
estimatePower <- FALSE

# Data is...
# Load TAMPOR cleanDat from file: #2918 or 3022 rows.
datafile <- paste(Rdatadir,tissue,"TAMPOR_data_outliersRemoved.Rds",sep="/")
cleanDat <- readRDS(datafile)
cleanDat <- log2(cleanDat)
cleanDat[1:5,1:5]
dim(cleanDat)

# Load combined sample info.
traitsfile <- paste(Rdatadir,tissue,"Combined_Cortex_Striatum_traits.Rds",sep="/")
sample_info <- readRDS(traitsfile)
sample_info[1:5,1:5]
dim(sample_info)

# Allow parallel WGCNA calculations:
allowWGCNAThreads()
parallelThreads <- 9
clusterLocal <- makeCluster(c(rep("localhost", parallelThreads)), type = "SOCK")
registerDoParallel(clusterLocal)

## Determine soft power, beta.
# Vector of powers to test:
powers <- seq(4, 20, by = 1.0)

# Soft Power selection
if (estimatePower==TRUE){
  sft <- pickSoftThreshold(t(cleanDat),
                           powerVector = powers, corFnc = "bicor",
                           blockSize = 15000, verbose = 3, networkType = "signed")
  
  # Create table. 
  mytable <- round(sft$fitIndices,2)
  mytable$truncated.R.sq <- NULL
  table <- tableGrob(mytable, rows = NULL)
  grid.arrange(table)
  
  # Save table as tiff.
  file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology_Table.tiff")
  ggsave(file,table)
  
  # Figure. ggplotScaleFreeFit() generates three plots.
  plots <- ggplotScaleFreeFit(sft)
  plots$Grid
  
  # Save as tiff.
  file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology.tiff")
  ggsave(file,plots$Grid)
  
  # Save plots and table as PDF.
  #plot_list <- list(table,plots$ScaleFreeFit, plots$MeanConnectivity)
  #file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology.pdf")
  #ggsavePDF(plot_list,file)
}

#-------------------------------------------------------------------------------
#' ## Build the WG(P)CNA Network.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Load blockwiseModules Parameters.
files <- list.files(Rdatadir,pattern = "Stats")
files
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
file <- paste0(outputfigsdir,"/",outputMatName,"Network_Parameters.tiff")
ggsave(file,table,width = 9, height = 1)

# Network parameters:
power <- 12 #9
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
net <- blockwiseModules(t(cleanDat),
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

#-------------------------------------------------------------------------------
#' ## Enforce module preservation.
#-------------------------------------------------------------------------------
#' Remove modules that are not preserved (i.e. have insignificant module 
#' preservation statistics.)
#' 
#+ eval = FALSE

# Input for NetRep:
r <- bicor(t(cleanDat))
adjm <- ((1+r)/2)^power # Signed network. 

data_list <- list(data = t(cleanDat))
correlation_list <- list(data = r)       
network_list <- list(data = adjm)       
module_labels <- net$colors     
names(module_labels) <- rownames(cleanDat)

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
  #nThreads, 
  #nPerm = 1000, # nPerm will be determined by the function. 
  null = "overlap", 
  alternative = "greater", 
  simplify = TRUE,
  verbose = TRUE)
# Collect stats. 
preservation <- preservation[c("observed","p.values")]

# Get the maximum permutation test p-value.
maxPval <- apply(preservation$p.values, 1, function(x) max(x,na.rm=TRUE))
nModules.original <- length(unique(net$colors))-1
nModules.original

# Modules removed if adjusted pvalue is greater than alpha.
alpha = 0.05
modules_out <- names(maxPval)[maxPval>alpha/nModules.original]
length(modules_out)
names(net$MEs)[grepl(paste(modules_out,collapse="|"),names(net$MEs))]<-"MEgrey"

# Drop nsig Modules (Set color to grey)
net$colors[net$colors %in% modules_out] <- "grey"
nModules <- length(unique(net$colors))-1
params$nModules <- nModules

# Check percent grey.
percent_grey <- round(100*sum(net$colors=="grey")/length(net$colors),2)
print(paste("Percent grey nodes =", percent_grey))
params$PercentGrayNodes <- percent_grey

# Table of key network stats.
mytable <- data.frame(
  nNodes = sum(net$colors!="grey"),
  PercentGrey = params$PercentGrayNodes,
  nModules = params$nSigModules,
  MedianCoherence = round(params$medianModCoherence,3),
  NetworkModularity = round(params$q2,3))
table <- tableGrob(mytable, rows = NULL)
grid.arrange(table)
file <- paste0(outputfigsdir,"/",outputMatName,"Key_Network_Stats.tiff")
ggsave(file,table,width = 6, height = 1)

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
MEs <- moduleEigengenes(t(cleanDat), 
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
moduleColors <- data.frame(protein = rownames(cleanDat),
                           module  = net$colors)

#-------------------------------------------------------------------------------
#' ## Visualize the WGCNA adjaceny network.
#-------------------------------------------------------------------------------

# Generate dissimilarity matrix (1-TOM(adj)).
# Exclude grey proteins.
idx <- net$colors != "grey"
diss <- 1 - TOMsimilarityFromExpr(
  datExpr = t(cleanDat[idx,]), 
  corType = "bicor",
  networkType = "signed",
  power = power,
  TOMType = "signed",
  TOMDenom = "mean",
  verbose = 0)

# Perform hierarchical clustering.   
GeneNames <- rownames(cleanDat)
colnames(diss) <- rownames(diss) <- GeneNames[idx]
hier <- hclust(as.dist(diss), method = "average")

# Generate dendrogram with module color labels.
file <- paste0(outputfigsdir,"/",outputMatName,"Network_Dendrogram.tiff")
tiff(file)
plotDendroAndColors(hier, 
                    net$colors[idx], 
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

dev.off()

#Visualize the TOM plot.
#diag(diss) <- NA
#diss[lower.tri(diss)] <- NA
#sizeGrWindow(7,7)
#TOMplot(diss, hier, as.character(net$colors[idx]))

# Visualize clustering by MDS.
#cmdMDS <- cmdscale(as.dist(diss))
#df <- data.frame(x     = cmdMDS[,1],
#                 y     = cmdMDS[,2])
#df$color = net$colors[match(rownames(df),rownames(cleanDat))]

# Examine MDS WGNA dissimilarity matrix.  
#plot <- ggplot(df,aes(x=x,y=y)) + geom_point(aes(colour = color),size = 1.5) + 
#  scale_color_manual(values = df$color) + theme(legend.position = "none") + 
#  ggtitle("Module MDS Plot") + xlab("Leading LogFC dim 1") + ylab("Leading LogFC dim 2") +
#  theme(
#    plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
#   axis.title.x = element_text(color = "black", size = 11, face = "bold"),
#    axis.title.y = element_text(color = "black", size = 11, face = "bold")
#  )

#plot

#-------------------------------------------------------------------------------
#' ## Evaluate Network quality (modularity) of WGCNA partition.
#-------------------------------------------------------------------------------

# Calculate network modularity.
# Calculate module cohesiveness (percent variance explained). 

# Calculate the adjacency network.
r <- bicor(t(cleanDat))
adjm <- ((1+r)/2)^power #signed.
#adjm <- abs(r)^power     #un-signed.

# Create igraph object. 
graph <- graph_from_adjacency_matrix(
  adjmatrix = adjm, 
  mode = c("undirected"), 
  weighted = TRUE, 
  diag = FALSE)

# Calculate modularity, q.
membership <- as.numeric(as.factor(net$colors))
q1 <- modularity(graph, membership, weights = edge_attr(graph, "weight"))
q1

# Without "grey" nodes.
v <- rownames(cleanDat)[!net$colors=="grey"]
subg <- induced_subgraph(graph,v)
membership <- as.numeric(as.factor(net$colors))
membership <- membership[!net$colors=="grey"]
q2 <- modularity(subg, membership, weights = edge_attr(subg, "weight"))
q2

# Module Coherence (proportion of variance explained; PVE)
# Exclude grey.
pve <- propVarExplained(datExpr=t(cleanDat), colors=net$colors, 
                        MEs=net$MEs, corFnc = "cor", corOptions = "use = 'p'")
median(pve[!names(pve)=="PVEgrey"])
mean(pve[!names(pve)=="PVEgrey"])

# Examine distribution of Module sizes.
df <- data.frame(Size = as.matrix(table(net$colors)))
df$Color <- rownames(df)
df <- df[order(df$Size,decreasing = TRUE),]
df <- df[!df$Color=="grey",]
df$Color <- factor(df$Color,levels = unique(df$Color))

# Plot
plot <- ggplot(df, aes(x = Color, y = Size, fill = Color)) + geom_col() +
  ylab("Module Size") + scale_fill_manual(values=levels(df$Color)) + 
  xlab("Module (Color)") + scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
        axis.title.x = element_text(color = "black", size = 11, face = "bold"),
        axis.text.x = element_blank(), axis.ticks.x=element_blank(),
        axis.title.y = element_text(color = "black", size = 11, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5))
plot

# Build stats annotation.
stats <- as.matrix(summary(df$Size))
df <- add_column(as.data.frame(stats),rownames(stats),.before=1)
tt <- ttheme_default(base_size = 11, core=list(bg_params = list(fill = "white")))
tab <- tableGrob(df, cols = NULL, rows=NULL, theme=tt)
g <- gtable_add_grob(tab, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                     t = 1, b = nrow(tab), l = 1, r = ncol(tab))

# Calculate ranges for positioning the annotation layer at the top right corner. 
pos <- ggranges(plot)$TopRight
xmin <- pos$xmin
xmax <- pos$xmax
ymin <- pos$ymin
ymax <- pos$ymax

# Add annotation.
# fixme: Should remove grey bar!
plot <- plot + annotation_custom(g, xmin = xmin, xmax = xmax, 
                                 ymin = ymin, ymax = ymax)
plot

# Save as tiff.
file <- paste0(outputfigsdir,"/",outputMatName,"Module_Sizes.tiff")
ggsave(file,plot)

#-------------------------------------------------------------------------------
#' ## Show that the expression of interacting proteins are highly correlated. 
#-------------------------------------------------------------------------------

# Illustrate the correaltion of two proteins.
saveplots = FALSE

# Calculate bicor adjacency matrix.
R <- cor(t(cleanDat))
diag(R) <- NA
R[lower.tri(R)] <- NA
R <- na.omit(melt(R))
colnames(R) <- c("Protein1","Protein2","Bicor")
R <- R[order(R$Bicor,decreasing=TRUE),]
prot1 <- "Wdr7|Q920I9"
prot2 <- "Rogdi|Q3TDK6"

# Generate a protein scatter plot.
#ggplotProteinScatterPlot(cleanDat,prot1,prot2)
#file <- "Protein_ScatterPlot.tiff"
#ggsave(file)


# Load simple interaction file (SIF). An edge list of known PPIs among all 
# identified proteins (Cortex + striatum). 
files <- c(
  paste(datadir,"PPI Network","Cortex_SIF_031019.txt",sep="/"),
  paste(datadir,"PPI Network","Striatum_SIF_031219.txt",sep="/")
)
sif_list <- lapply(as.list(files),function(x) read.table(x,header = TRUE, sep=","))

# Merge SIFS.
sif <- do.call(rbind,sif_list)
head(sif)

# Map rownames of cleanDat (UniprotIDs) to Entrez.
Uniprot <- sapply(strsplit(rownames(cleanDat),"\\|"),"[",2)
Entrez <- mapIds(org.Mm.eg.db, keys=Uniprot, column="ENTREZID", 
                 keytype="UNIPROT", multiVals="first")
length(Entrez)

# Subset SIF Keep only those proteins mapped to cleanDat:Entrez.
idx <- sif$EntrezA %in% Entrez & sif$EntrezB %in% Entrez
sif <- sif[idx,]

# Subset Entrez IDs. Keep only those mapped to sif:Entrez.
idx <- Entrez %in% sif$EntrezA | Entrez  %in% sif$EntrezB
Entrez <- Entrez[idx]
length(Entrez)

# Create iGraph.
edgeList <- cbind(sif$EntrezA,sif$EntrezB)
PPIgraph <- graph_from_edgelist(edgeList, directed = FALSE)

# Calculate bicor adjacency matrix.
dat <- cleanDat
Uniprot <- sapply(strsplit(rownames(dat),"\\|"),"[",2)
rownames(dat) <- mapIds(org.Mm.eg.db, keys=Uniprot, column="ENTREZID", 
                        keytype="UNIPROT", multiVals="first")
  
R <- bicor(t(dat))
#idx <- melt(upper.tri(R,diag=FALSE))
#idx <- idx$value
diag(R) <- NA
R <- na.omit(melt(R))
colnames(R) <- c("ProteinA","ProteinB","Bicor")
R$ID <- paste(R$ProteinA,R$ProteinB,sep="_")
head(R)

# Add TRUE/FALSE if known interacting pair.
sif$IntA <- paste(sif$EntrezA,sif$EntrezB,sep="_")
sif$IntB <- paste(sif$EntrezB,sif$EntrezA,sep="_")
R$InteractionPair <- as.numeric(R$ID %in% sif$IntA | R$ID %in% sif$IntB)
table(R$InteractionPair)

# The data is really large, sample N data points from 
# the distributions of interacting and non-interacting proteins.
# set.seed to insure reproducible random samples.
set.seed(7)
x <- cbind(sample(R$Bicor[R$InteractionPair==0],2000),FALSE)
set.seed(7)
y <- cbind(sample(R$Bicor[R$InteractionPair==1],2000),TRUE)
df <- as.data.frame(rbind(x,y))
colnames(df) <- c("bicor","group")
df$group <- as.factor(df$group)

# Generate a plot.
plot <- ggplot(df, aes(x=group, y=bicor, fill=group)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=20, outlier.size=1)

# Annotate with Wilcoxon Rank sum p-value.
wrs <- wilcox.test(y[,1], x[,1],alternative = "greater")
pval <- paste("P =", formatC(wrs$p.value,format="e",digits=2))
plot <- plot + annotate("text",x = 1.5, y = 1, label = "*", size = 12, color="black") + 
  ggtitle(paste("Wilcox Rank Sum", pval)) + 
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold")
  )
plot  

# Save to file.
if (saveplots==TRUE){
  file <- paste0(outputfigsdir,"/",outputMatName,"InteractingProteinBicorBoxplot.pdf")
  ggsavePDF(plot,file)
}

# Save as tiff.
file <- paste0(outputfigsdir,"/",outputMatName,"Interacting_Protein_Bicor.tiff")
ggsave(file,plot)

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
numericMeta$Tissue.Sample.Model <- paste(numericMeta$TissueType,numericMeta$Sample.Model,sep=".")

# Pool WT within a tissue.
numericMeta$Tissue.Sample.Model[grepl("Cortex.WT",numericMeta$Tissue.Sample.Model)] <- "Cortex.WT"
numericMeta$Tissue.Sample.Model[grepl("Striatum.WT",numericMeta$Tissue.Sample.Model)] <- "Striatum.WT"
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
#' ## Calculate Module EigenProteins.
#-------------------------------------------------------------------------------
#' 
#' Eigenproteins (also, eigengenes) are the first principle component of a 
#' module. They are a summary of a module. An eigenprotein is not an actual 
#' protein.
#' 
#+ eval = FALSE

# Calculate Module EigenProteins (MEs)
MEs <- tmpMEs <- data.frame()
MEList <- moduleEigengenes(t(cleanDat), colors = net$colors)
MEs <- orderMEs(MEList$eigengenes)

# Remove prefix.
colnames(MEs) <- gsub("ME", "", colnames(MEs)) 
rownames(MEs) <- rownames(numericMeta)
MEs[1:5,1:5]

# Write to excel.
file <- paste0(outputtabsdir,"/",tissue,"_WGCNA_Analysis_ModuleEigenprotein_expression.xlsx")
write.excel(MEs,file,rowNames = TRUE)

#-------------------------------------------------------------------------------
#' ## Determine module membership (kME). 
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Module membership (kME) is the strength of association (correlation) between a 
# proteins expression vector and that of module EigenProteins.
tmpMEs <- net$MEs
colnames(tmpMEs) <- paste("ME", colnames(MEs), sep = "")

# The protein-wise correlation with module Eigenproteins. 
kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc = "bicor")
colnames(kMEdat) <- gsub("kME","",colnames(kMEdat))
kMEdat[1:5,1:5]

# Clean up a little.
idx <- match(rownames(kMEdat),rownames(cleanDat))
kMEdat <- add_column(kMEdat,Protein = rownames(kMEdat),.before = 1)
kMEdat <- add_column(kMEdat,Module = net$colors[idx],.after = 1)
rownames(kMEdat) <- NULL

# Write to excel with custom function write.excel().
file <- paste0(outputtabsdir,"/",tissue,"_WGCNA_Analysis_ModuleMembership_kME.xlsx")
write.excel(kMEdat,file)

#-------------------------------------------------------------------------------
#' ## Identify module Hub proteins.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Should plots be saved?
saveplots = FALSE

# Modify traits$Sample.Model for grouping all WT's together. 
traits <- sample_info
traits_temp <- traits
# Column groupings are determined by traits$Sample.Model. 
traits_temp$Sample.Model <- paste(traits_temp$Tissue,
                                  traits_temp$Sample.Model,sep=".")
traits_temp$Sample.Model[grepl("Cortex.WT",traits_temp$Sample.Model)] <- "Cortex.WT"
traits_temp$Sample.Model[grepl("Striatum.WT",traits_temp$Sample.Model)] <- "Striatum.WT"

# Find top n hub proteins.
n <- 3
HubProteins <- melt(kMEdat,id=c("Protein","Module"))
HubProteins <- subset(HubProteins,HubProteins$Module==HubProteins$variable)
HubProteins$variable <- NULL
colnames(HubProteins) <- c("Protein","Module","kME")
HubProteins <-  HubProteins %>% group_by(Module) %>% top_n(n, kME)

# Generate boxplots to examine expresion of hub proteins.
plot_list <- ggplotProteinBoxes(cleanDat,
                                interesting.proteins=as.character(HubProteins$Protein),
                                dataType="Relative Abundance",traits=traits_temp,
                                order = c(2,7,4,6,5,8,1,9,3,10))

# Example plot.
plot_list[[1]]

# Build a df with statistical results.
# colnames should match x axis of plot.
stats <- lapply(results,function(x) 
  as.data.frame(cbind(Uniprot=x$Uniprot,FDR=x$FDR)))
names(stats) <- names(results)
df <- stats %>% reduce(left_join, by = "Uniprot")
colnames(df)[c(2:ncol(df))] <- names(stats)

# Annotate rows as gene|uniprot
Uniprot <- df$Uniprot
Gene <- mapIds(
  org.Mm.eg.db,
  keys=Uniprot,
  column="SYMBOL", 
  keytype="UNIPROT", 
  multiVals="first")
rownames(df) <- paste(Gene,Uniprot,sep="|")
df$Uniprot <- NULL
stats <- df

# Define proteins of interest.
prots <- HubProteins$Protein
length(prots)

# Generate boxplots to examine expresion of hub proteins.
plot_list <- ggplotProteinBoxes(cleanDat,
                                interesting.proteins=as.character(prots),
                                dataType="Relative Abundance",traits=traits_temp,
                                order = c(2,7,4,6,5,8,1,9,3,10))
# Example plot.
plot <- plot_list[[1]]
plot
annotate_stars(plot,stats)

# Loop to add stars.
plot_list <- lapply(plot_list,function(x) annotate_stars(x,stats))

# Loop to add custom colors.
colors <- rep(c("gray","yellow","blue","green","purple"),each=2)
plot_list <- lapply(plot_list,function(x) x + scale_fill_manual(values=colors))

# Example plot:
plot_list[[1]]

# Loop to add module assignment annotation.
results_modules <- data.frame(geneNames     = rownames(cleanDat),
                              dynamicColors = net$colors)
for (i in 1:length(plot_list)){
  plot <- plot_list[[i]]
  namen <- results_modules$dynamicColors[match(plot$labels$title,results_modules$geneNames)]
  plot <- plot + ggtitle(paste(plot$labels$title,namen,sep="|"))
  plot_list[[i]] <- plot
}

# Example plot.
plot_list[[1]]

# Store as list.
mhplots <- split(plot_list,HubProteins$Module)

# Print to pdf.
if (saveplots==TRUE){
  file <- paste0(outputfigsdir,"/",tissue,"_WGCNA_Analysis_HubProtein_Expression.pdf")
  ggsavePDF(plot_list,file)
}

## Write hub data to excel. 
file <- paste0(outputtabsdir,"/",tissue,"_WGCNA_Analysis_Module_Hubs.xlsx")
write.excel(as.data.frame(HubProteins),file)

# Function to loop and save top three hubs of desired color.
my_func <- function(color){
  for (i in 1:3){
    plot <- mhplots[[color]][[i]]
    file <- paste0(outputfigsdir,"/",outputMatName,color,"_",i,"_Module_Hub.tiff")
    ggsave(file,plot)
  }
}

# Use lapply to generate plots for desired colors.

#-------------------------------------------------------------------------------
#' ## Determine protein-wise correlation with traits (Gene Significance). 
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Calculate protein-wise correlation (Pearson) with traits. 
# In the WGCNA parlance this is also known as Gene significance (GS).
geneSignificance <- cor(sapply(numericMeta[, numericIndices], as.numeric), 
                        t(cleanDat), use = "pairwise.complete.obs")
# Add rownames. 
rownames(geneSignificance) <- colnames(numericMeta)[numericIndices]
geneSignificance[1:5,1:5]

# Convert correlation coefficients to colors. 
geneSigColors <- t(numbers2colors(t(geneSignificance), signed = TRUE, 
                                  lim = c(-1, 1), naColor = "black"))
# Add rownames (protein name).
rownames(geneSigColors) <- colnames(numericMeta)[numericIndices]
geneSigColors[1:5,1:5]

# What are the top correlates (biomarkers?) of each trait?
n <- 3
GS <- abs(as.data.frame(t(geneSignificance)))
GS <- melt(add_column(GS,Protein=rownames(GS),.before=1))
colnames(GS)[c(2,3)] <- c("Trait","cor")     
biomarkers <- GS %>% group_by(Trait) %>% top_n(n, cor)

# Write to excel.
data <- list(abs(t(geneSignificance)),as.matrix(biomarkers))
names(data) <- c("geneSignificance","Biomarkers")
file <- paste0(outputtabsdir,"/",tissue,"_WGCNA_Analysis_geneSignificance.xlsx")
write.excel(data,file,rowNames = TRUE)

#-------------------------------------------------------------------------------
#' ## Examine module significance (MS).
#-------------------------------------------------------------------------------

# Module significance is the average GS of all proteins in a module. 
# Which module is the most important for a given trait (e.g. Ube3a_KO)?

# Should plots be saved?
saveplots = FALSE

# Module Assignment summary.
geneNames <- rownames(cleanDat)
dynamicColors <- net$colors
results_modules <- as.data.frame(cbind(geneNames,dynamicColors))

# Add Module Assignment to GS data.
idx <- match(GS$Protein,results_modules$geneNames)
GS$Module <- results_modules$dynamicColors[idx]
# Remove grey.
GS <- GS[!GS$Module=="grey",]
head(GS)

# Divide GS into list by Traits.
GS_list <- split(GS, GS$Trait)

# Single plot.
colnames(numericMeta)[numericIndices]
trait <- "Cortex.ASD"
x <- GS_list[[trait]]$cor
g <- GS_list[[trait]]$Module
ggplotModuleSignificanceBoxplot(x,g,trait,stats=FALSE)

# Plots for main interesting traits.
trait_list <- as.list(names(GS_list)[21:length(GS_list)])
plots <- lapply(trait_list, function(x) 
  ggplotModuleSignificanceBoxplot(GS_list[[x]]$cor,GS_list[[x]]$Module,x,stats=FALSE))
names(plots) <- unlist(trait_list)
msplots <- plots

# Save to pdf.
if (saveplots==TRUE){
  file <- paste0(outputfigsdir,"/",tissue,"_WGCNA_Anaysis_MS_Boxplots.pdf")
  ggsavePDF(plots=MSplot,file)
}

# Calculate MS as the mean GS for a module.
MS <- subset(GS) %>% group_by(Module, Trait) %>% dplyr::summarise(MS = mean(cor))

# Seperate into elements of a list by trait.
MS_list <- split(MS,MS$Trait)
names(MS_list) <- unique(MS$Trait)

# Write to excel.
file <- paste0(outputtabsdir,"/",tissue,"_WGCNA_Analysis_ModuleSignificance.xlsx")
write.excel(MS_list,file)

#-------------------------------------------------------------------------------
#' ## Module significance - Which modules are enriched for DEPs?
#-------------------------------------------------------------------------------

## Insure statistical results are loaded.
# InterBatch statistical comparisons with EdgeR GLM:
file <- paste0(rootdir,"/","Tables/Combined/v14_TAMPOR/Combined_TMT_Analysis_TAMPOR_GLM_Results.xlsx")
results <- lapply(as.list(c(1:8)),function(x) read_excel(file,x))
names(results) <- excel_sheets(file)

# Build a df with statistical results.
# colnames should match x axis of plot.
stats <- lapply(results,function(x) 
  as.data.frame(cbind(Uniprot=x$Uniprot,FDR=x$FDR)))
names(stats) <- names(results)
df <- stats %>% reduce(left_join, by = "Uniprot")
colnames(df)[c(2:ncol(df))] <- names(stats)

# Annotate rows as gene|uniprot
Uniprot <- df$Uniprot
Gene <- mapIds(
  org.Mm.eg.db,
  keys=Uniprot,
  column="SYMBOL", 
  keytype="UNIPROT", 
  multiVals="first")
rownames(df) <- paste(Gene,Uniprot,sep="|")
df$Uniprot <- NULL
stats <- df

# Tidy up pvalues.
x <- melt(t(stats))
colnames(x) <- c("Condition","Protein","FDR")
x$Module <- results_modules$dynamicColors[match(x$Protein,results_modules$geneNames)]
p <- x$FDR
g <- x$Module
df <- as.data.frame(cbind(g,p))

# Calculate summary statistics.
df.stats <- subset(df) %>% group_by(g) %>% 
  dplyr::summarise(
    Module.Size = length(p),
    Sig.Count = sum(p<0.05))
df.stats$Total.sig <- sum(df.stats$Sig.Count)
df.stats$N <- sum(df.stats$Module.Size)
df.stats <- df.stats[!df.stats$g=="grey",] # Remove grey.
df.stats <- na.omit(df.stats)

# Perform hypergeometric test for enrichment. 
# FDR is Bonferroni!
df.stats$p.val <- phyper(df.stats$Sig.Count,
                         df.stats$Total.sig,
                         df.stats$N-df.stats$Sig.Count,
                         df.stats$Module.Size, lower.tail = FALSE)
df.stats$FDR <- p.adjust(df.stats$p.val,method = "bonferroni")
df.stats$FE.pval <- -log(df.stats$p.val)
df.stats <- df.stats[order(df.stats$p.val),]
df.stats$g <- factor(df.stats$g, levels = df.stats$g)

# Plot
# Red line is sig at BH p.adjust. 
plot <- ggplot(df.stats, aes(x = g, y = FE.pval, fill = g)) + 
  geom_col() + 
  ylab("FE P-value") + scale_fill_manual(values=levels(df.stats$g)) + 
  xlab("Module (Color)") + scale_y_continuous(expand = c(0, 0)) + 
  geom_hline(yintercept = -log(0.05/dim(df.stats)[1]), 
             color = "red", linetype = "dashed") + 
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1))
plot

# Save as tiff.
file <- paste0(outputfigsdir,"/",outputMatName,"Module_Enrichment_DEPs.tiff")
ggsave(file,plot)

#-------------------------------------------------------------------------------
#' ## Determine module EigenProtein cor with traits (EigenProtein significance). 
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Each ME v trait correlation network could be vizualized...
# Calculate Eigengene-trait correlations:
# Warnings are okay.
MEcorTraits <- bicorAndPvalue(MEs, numericMeta[, numericIndices])
names(MEcorTraits) #nObs is Number of obs. 

# Clean up the data.
# MEcorTraits is a list. Melt each df in the list.
data_temp <- lapply(MEcorTraits,function(x) melt(x))

# Loop to add column names.
for (i in 1:length(data_temp)){
  colnames(data_temp[[i]]) <- c("Module","Trait",names(data_temp)[i])
}

# Merge the elements of the list.
df <- data_temp %>% reduce(left_join, by = c("Module","Trait"))

# Which module is the top correlates of every traits?
n <- 3
topModules <- df %>% group_by(Trait) %>% top_n(n, -log(p))

# Write to excel.
file <- paste0(outputtabsdir,"/",tissue,"_WGCNA_Analysis_topModules_corTraits.xlsx")
write.excel(as.data.frame(topModules),file)

# Convergence... a module that is highly correlated with several traits.
# Also, the degree of a module node in the module-trait network.
MSdm <- MEcorTraits$bicor
colnames(MSdm)

# Get comparisons of interest. 
cols <- c(
  "Cortex.KO.Shank2","Striatum.KO.Shank2",
  "Cortex.KO.Shank3","Striatum.KO.Shank3",
  "Cortex.HET.Syngap1","Striatum.HET.Syngap1",
  "Cortex.KO.Ube3a","Striatum.KO.Ube3a")
idx <- match(cols,colnames(MSdm))
dm <- MSdm[,idx]

ESdf <- as.data.frame(dm)
rownames(ESdf) <- gsub("ME","",rownames(ESdf))
ESdf$Degree <- rowSums(abs(ESdf))
sum(ESdf$Degree)
max(ESdf$Degree)

# Plot.
# Fixme: need to clean this plot up. 
df <- data.frame(Color  = rownames(ESdf),
                 Degree = ESdf$Degree)
x <- df$Degree
g <- df$Color
ggplot(df,aes(x=Color,y=Degree,fill=Color)) + geom_col()

# Calculate module coherence (Percent variation explained by ME).
pve <- as.data.frame(
  propVarExplained(datExpr=t(cleanDat), colors=net$colors,
                   MEs=net$MEs, corFnc = "cor", corOptions = "use = 'p'"))
pve <- add_column(pve, Module = gsub("PVE","",rownames(pve)), .before=1)
rownames(pve) <- NULL
colnames(pve)[2] <- "PVE"

#-------------------------------------------------------------------------------
#' ## VerboseBoxplots - vizualize module expression across traits.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Should plots be saved?
saveplots = FALSE

# Calculate Module EigenProteins.
MEs <- tmpMEs <- data.frame()
MEList <- moduleEigengenes(t(cleanDat), colors = net$colors)
MEs <- orderMEs(MEList$eigengenes)
# let's be consistent in case prefix was added, remove it.
colnames(MEs) <- gsub("ME", "", colnames(MEs)) 
rownames(MEs) <- rownames(numericMeta)

# Insure traits are in matching order.
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
KWsigModules <- KW_results$Module[KW_results$p.adj<0.05]
print(paste("nModules with p.adj < 0.05 =", length(KWsigModules)))
print(paste("nModules with p.adj < 0.10 =", length(KW_results$Module[KW_results$p.adj<0.1])))
print(paste("nModules with FDR < 0.05 =", length(KW_results$Module[KW_results$FDR<0.05])))

# Split Module EigenProtein (MEs) dm into a list of column vectors for lapply. 
ME_list <- split(as.matrix(MEs), rep(1:ncol(MEs), each = nrow(MEs)))
names(ME_list) <- colnames(MEs)

# Add vector of groups
ME_list <- lapply(ME_list,function(x) data.frame(x = x, groups = groups))

# Example plot. Specify stats=TRUE to return the KW and Dunn test results.
# If overall p-value (unadjusted KW) is significant , then title is red.
# Signicance of Dunn's post-hoc test comparisons to pooled WT group 
# are indicated with astricks.
x <- ME_list[[1]]$x
g <- ME_list[[1]]$groups
color <- names(ME_list)[1]
color
# Define levels for order of bars in plot.
#levels <- c("WT.Cortex","ASD.Cortex","WT.Striatum","ASD.Striatum")
levels <- c("WT.Cortex","WT.Striatum",
            "Cortex.KO.Shank2","Striatum.KO.Shank2",
            "Cortex.KO.Shank3","Striatum.KO.Shank3",
            "Cortex.HET.Syngap1","Striatum.HET.Syngap1",
            "Cortex.KO.Ube3a","Striatum.KO.Ube3a")

# Generate contrasts matrix for comparisons of interest.
#contrasts <- makePairwiseContrasts(list("WT.Cortex","WT.Striatum"),list("ASD.Cortex","ASD.Striatum"))
#contrasts
g1 <- list("WT.Cortex","WT.Striatum")
g2 <- list(
  c("Cortex.KO.Shank2","Cortex.KO.Shank3","Cortex.HET.Syngap1","Cortex.KO.Ube3a"),
  c("Striatum.KO.Shank2","Striatum.KO.Shank3","Striatum.HET.Syngap1","Striatum.KO.Ube3a"))
contrasts <- makePairwiseContrasts(g1,g2)
contrasts

# Generate plot.
plot_data <- ggplotVerboseBoxplot(x,g,levels,contrasts,color,stats=TRUE, 
                                  method = "dunn")
plot_data$plot
#plot_data$dunnett
plot_data$dunn

# Loop through ME_list and generate verboseBoxPlot.
# lapply wont work here because the name is not preserved when you call lapply()...
# method is the p.adj method for the Dunn's test p-value.
plot_data <- list()
for (i in 1:dim(MEs)[2]){
  x <- ME_list[[i]]$x
  g <- ME_list[[i]]$groups
  color <- names(ME_list)[[i]]
  plot <- ggplotVerboseBoxplot(x,g,levels,contrasts,color,stats=TRUE,method="dunn")
  plot_data[[i]] <- plot
  names(plot_data)[[i]] <- color
}

# Extract plots. 
plot_list <- sapply(plot_data,"[",1)

# Store as VBplots
vbplots <- plot_list

# Save all plots.
if (saveplots==TRUE){
  file <- paste0(outputfigsdir,"/",tissue,"_WGCNA_Analysis_VerboseBoxPlots.pdf")
  ggsavePDF(plot_list,file)
}

# Extract Dunnett test stats.
Dtest_stats <- sapply(plot_data,"[",3)
names(Dtest_stats) <- sapply(strsplit(names(Dtest_stats),"\\."),"[",1)

# Loop to add module column.
for (i in 1:length(Dtest_stats)){
  df <- Dtest_stats[[i]]
  namen <- names(Dtest_stats)[i]
  df <- add_column(df,Module = namen,.before=1)
  Dtest_stats[[i]] <- df
}

# Number of significant tests.
nsig <- do.call(rbind,lapply(Dtest_stats,function(x) sum(x$P.adj<0.05)))
idx <- match(KW_results$Module,rownames(nsig))
KW_results$DunnettTestNsig <- nsig[idx,]

# Save tiffs of top plots.
idx <- subset(KW_results,p.adj<0.1)$Module
for (i in 1:length(idx)){
  color <- idx[i]
  file <- paste0(outputfigsdir,"/",outputMatName,"Verbose_boxplot_",color,".tiff")
  ggsave(file,vbplots[[paste(color,"plot",sep=".")]])
}

#-------------------------------------------------------------------------------
#' ## VerboseScatter plots.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Should plots be saved?
saveplots = FALSE

# Calculate MEs
MEList <- moduleEigengenes(t(cleanDat), colors = net$colors)
MEs <- orderMEs(MEList$eigengenes)

# X data = kME data. 
MMdata <- signedKME(t(cleanDat), MEs, corFnc = "bicor")
colnames(MMdata) <- gsub("kME","",colnames(MMdata))
MMdata[1:5,1:5]
class(MMdata)
MMdata <- MMdata[,!colnames(MMdata) == "grey"] # Remove grey.

# Y data = Gene significnce. The absolute value of protein-wise cor to traits.
geneSignificance <- cor(sapply(numericMeta[, numericIndices], as.numeric), 
                        t(cleanDat), use = "pairwise.complete.obs")
rownames(geneSignificance) <- colnames(numericMeta)[numericIndices]
GSdata <- as.data.frame(abs(t(geneSignificance)))
GSdata[1:5,1:5]

# moduleGenes
geneNames <- rownames(cleanDat)
color <- net$colors
moduleGenes <- as.data.frame(cbind(geneNames,color))
head(moduleGenes)

# Plot an example with ggplotVerboseScatterPlot:
ggplotVerboseScatterPlot(MMdata,GSdata,moduleGenes,"blue","Ube3a_KO")

## Loop to generate all plots for a given group (trait).
modules <- colnames(MMdata)
plots <- list()
trait <- "Ube3a_KO"

# Loop
for (color in seq_along(modules)){
  plots[[color]] <- ggplotVerboseScatterPlot(MMdata,
                                             GSdata,
                                             moduleGenes,
                                             modules[color],
                                             trait,
                                             stats=TRUE)
}
names(plots) <- modules

# Gather stats. The function ggplotVerboseScatterPlot returns linear fit statistics.
stats <- as.data.frame(do.call(rbind,sapply(plots,"[",2)))
stats$MMnum <- c(1:nrow(stats))
stats <- stats[order(stats$p.value),]
head(stats)

# Significant correlations.
sub <- stats[stats$p.value<0.05,]
sub[order(sub$R2, decreasing=TRUE),]

# Exmaple of a significant correlation.
x <- strsplit(rownames(sub)[1],"\\|")[[1]][1]
plots[[x]]

## Create verbose scatter plots for multipe traits.

# Make the loop above into a function:
make_ggplotVerboseScatterPlots <- function(MMdata,GSdata,moduleGenes,modules,trait){
  for (color in seq_along(modules)){
    plots[[color]] <- ggplotVerboseScatterPlot(MMdata,GSdata,moduleGenes,modules[color],trait)
  }
  names(plots) <- modules
  return(plots)
}

# Define the traits of interest.
colnames(GSdata)
#all_traits <- c("ASD","Shank2_KO","Shank3_KO","Syngap1_KO","Ube3a_KO")
all_traits <- colnames(GSdata)[c(13:ncol(GSdata))]
modules <- colnames(MMdata)
plots <- list()

# Use lapply to generate plots.
plots <- lapply(as.list(all_traits),
                function(x) 
                  make_ggplotVerboseScatterPlots(MMdata,GSdata,moduleGenes,modules,x))
names(plots) <- all_traits

# plots is a list of all verbose scatter plots for a trait.
length(plots)
names(plots)

# Gather plots and stats.
plot_list <- lapply(plots,function(x) sapply(x,"[",1))
stats_list <- lapply(plots,function(x) do.call(rbind,sapply(x,"[",2)))
names(plot_list) <- all_traits

# Sort stats by R2.
stats_list <- lapply(stats_list,
                     function(x) x[order(as.data.frame(x)$R2, decreasing = TRUE),])

# Store as list of plots.
vsplots <- plot_list

# Write plots to file.
if (saveplots==TRUE){
  for (i in 1:length(plot_list)){
    plots <- plot_list[[i]]
    file <- paste0(outputfigsdir,"/",tissue,"_WGCNA_Analysis_",names(plot_list)[i],"_verboseScatterplots.pdf")
    ggsavePDF(plots,file)
  }
}

# Save stats to excel.
stats_list <- c(list(do.call(rbind,stats_list)),stats_list)
names(stats_list)[1] <- "All_Traits"
file <- paste0(outputtabsdir,"/",tissue,"_WGCNA_Analysis_VerboseScatter_Stats.xlsx")
write.excel(stats_list,file, rowNames = TRUE)

# top striatum.asd?
#foo <- stats[grepl("M.ASD",rownames(stats)),]

# Significant correlations?
stats <- as.data.frame(stats_list$All_Traits)
stats$p.adj <- p.adjust(stats$p.value,method = "bonferroni")
sub <- subset(stats,p.adj<0.05 & slope >0)
sub

# Turquoise plot. 
plot <- vsplots$Cortex.ASD$turquoise.plot
plot 

# Save as tiff.
file <- paste0(outputfigsdir,"/",outputMatName,"Turquose_Verbose_Scatter.tiff")
ggsave(file,plot)

#-------------------------------------------------------------------------------
#' ## Build module (Eigengene) network.
#-------------------------------------------------------------------------------

# WGCNA Function:
#file <- paste0(outputfigsdir,"/",outputMatName,"Module_EigenGene_Network.pdf")
#CairoPDF(file = file, width = 16, height = 12)
plotEigengeneNetworks(MEs[,!colnames(MEs)=="grey"], "Eigengene Network",
                      excludeGrey = TRUE, greyLabel = "grey",
                      marHeatmap = c(3, 4, 2, 2), 
                      marDendro = c(0, 4, 2, 0), 
                      plotDendrograms = TRUE, 
                      plotAdjacency = FALSE,
                      printAdjacency = FALSE,
                      xLabelsAngle = 90, 
                      heatmapColors = blueWhiteRed(50))
#dev.off()  


# Generate distance matrix. 
r <- bicor(MEs[!colnames(MEs)=="MEgrey"])
diss <- 1 - r

# ME network 
MEnet <- r
diag(MEnet) <- NA
MEnet_list <- melt(MEnet)

# Perform hierarchical clustering
hc <- hclust(as.dist(diss), method = "average")

# Generate dendrogram.
p1 <- ggdendrogram(hc, rotate=FALSE)
p1

# Prepare a df for generating colored bars.
df2 <- data.frame(
  cluster = cutreeDynamic(hc, distM = diss, method="tree", minClusterSize = 2, verbose = 0),
  module  = factor(hc$labels, levels=hc$labels[hc$order]))
df2$order <- match(df2$module,dendro_data(hc)$labels$label)
df2 <- df2[order(df2$order),]
df2$module <- factor(df2$module,levels = df2$module)
df2$module <- gsub("ME","", df2$module)
head(df2)

# Number of meta modules.
length(unique(df2$cluster))

# Meta modules.
meta_modules <- data.frame(
  protein = rownames(cleanDat),
  module = net$colors,
  meta_module = df2$cluster[match(net$colors,df2$module)])
# Convert NA to zero.
meta_modules$meta_module[is.na(meta_modules$meta_module)] <- 0

# Generate colored bars.
p2 <- ggplot(df2,aes(module, y = 1, fill=factor(cluster))) + geom_tile() +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        legend.position="none")
p2

# Combine.
gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)  

maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)

fig <- as.ggplot(grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5)))
fig

# Save as tiff.
file <- paste0(outputfigsdir,"/",outputMatName,"Meta_Modules.tiff")
ggsave(file,fig)


# Modularity of the meta module network.
collectGarbage()
r <- bicor(t(cleanDat))
adjm <- ((1+r)/2)^power #signed.

# Create igraph object. 
graph <- graph_from_adjacency_matrix(
  adjmatrix = adjm, 
  mode = c("undirected"), 
  weighted = TRUE, 
  diag = FALSE)

# Calculate modularity, q.
membership <- as.numeric(as.factor(meta_modules$meta_module))
q1 <- modularity(graph, membership, weights = edge_attr(graph, "weight"))
q1

# Without grey modules.
v <- rownames(cleanDat)[!meta_modules$meta_module==0]
subg <- induced_subgraph(graph,v)
membership <- as.numeric(as.factor(meta_modules$meta_module))
membership <- membership[!meta_modules$meta_module==0]
q2 <- modularity(subg, membership, weights = edge_attr(subg, "weight"))
q2

#-------------------------------------------------------------------------------
#' ## GO Enrichment analysis of Modules (colors).
#-------------------------------------------------------------------------------
#+ eval = FALSE

## Prepare a matrix of class labels (colors) to pass to enrichmentAnalysis(). 
geneNames <- rownames(cleanDat)
dynamicColors <- net$colors
results_modules <- as.data.frame(cbind(geneNames,dynamicColors))

# Reshape the module colors data.
colors <- as.data.frame(results_modules)
mytable <- table(colors)
colors <- as.data.frame.matrix(mytable)
head(colors)

# Convert 0 and 1 to column names. 
logic <- colors == 1 # 1 will become TRUE, and 0 will become FALSE.
# Loop through each column to replace 1 with column header (color).
for (i in 1:ncol(logic)){
  col_header <- colnames(colors)[i]
  colors[logic[,i],i] <- col_header
}

## Map Uniprot IDs to Entrez.
# Get Uniprot IDs and gene symbols from rownames
#Uniprot_IDs <- as.character(rownames(colors))
Uniprot_IDs <- as.character(colsplit(rownames(colors),"\\|",c("Symbol","UniprotID"))[,2])
# Map Uniprot IDs to Entrez IDs:
entrez <- mapIds(org.Mm.eg.db, 
                 keys = Uniprot_IDs, 
                 column = "ENTREZID", 
                 keytype = "UNIPROT", 
                 multiVals = "first")

# Insure that colors is a matrix.
colors <- as.matrix(colors)
head(colors)

# look at the number of genes assigned to each cluster. 
table(colors)

# The colors matrix and vector of cooresponding entrez IDs 
# will be passed to enrichmentAnalysis().

# Build a GO annotation collection:
if (!exists(deparse(substitute(musGOcollection)))){
  musGOcollection <- buildGOcollection(organism = "mouse")
}

# Creates some space by clearing some memory.
collectGarbage()

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
  ignoreLabels = "0")

# Create some space by clearing some memory.
collectGarbage()

# Collect the results. 
results_GOenrichment <- list()
for (i in 1:length(GOenrichment$setResults)){
  results_GOenrichment[[i]] <- GOenrichment$setResults[[i]]$enrichmentTable
}

# Combined result
GO_result <- do.call(rbind,results_GOenrichment)
idx <- GO_result$Bonferroni<0.05
nsig <- length(idx[idx==TRUE])
print(paste("There are", nsig, "GO terms that exhibit significant enrichment among all modules."))

# Top Term for each module.
topGO <- subset(GO_result,rank==1)

# Number of modules with sig. GO enrichment.
table(topGO$FDR<0.05)[2]         
table(topGO$Bonferroni<0.05)[2]
sigGO <- subset(topGO,FDR<0.05)

# These are the unique sig. GO terms.
unique(sigGO$shortDataSetName)

# Add to results_GOenrichment list()
results_GOenrichment <- c(list(topGO),results_GOenrichment)

# Add names to list of results. 
names(results_GOenrichment) <- c("topGO",gsub(" FDR","",colnames(colors)))

# Write results to file. 
file <- paste0(outputtabsdir,"/",tissue,"_WGCNA_Analysis_Module_GOenrichment_Results.xlsx")
write.excel(results_GOenrichment,file)

#-------------------------------------------------------------------------------
#' ## Module summary.
#-------------------------------------------------------------------------------

# Module summary data frame. 
module_summary <- KW_results

# Add module sizes
df <- moduleColors %>% group_by(module) %>% dplyr::summarize(size = length(module))
idx <- match(module_summary$Module,df$module)
size <- df$size[idx]
module_summary <- add_column(module_summary,size,.after = 1)

# Add PVE
idx <- match(module_summary$Module,pve$Module)
PVE <- pve$PVE[idx]
module_summary <- add_column(module_summary,PVE,.after = 2)

# Add topGO annotation and pvalue. 
idx <- match(module_summary$Module,topGO$class)
module_summary$topGO <- topGO$shortDataSetName[idx]
module_summary$GO.p.val <- topGO$pValue[idx]
module_summary$GO.p.adj <- topGO$Bonferroni[idx]
module_summary$GO.FDR <- topGO$FDR[idx]

# Genes
genes <- do.call(rbind,
                 lapply(split(rownames(cleanDat),net$colors),
                        function(x) paste(x,collapse=",")))
idx <- match(module_summary$Module,rownames(genes))
module_summary$Proteins <- genes[idx]

# Hubs
hubs <- do.call(rbind,
                lapply(split(HubProteins$Protein,HubProteins$Module),
                       function(x) paste(x,collapse=",")))
idx <- match(module_summary$Module,rownames(hubs))
module_summary$Hubs <- hubs[idx]

# Write to excel.
file <- paste0(outputtabsdir,"/",tissue,"_WGCNA_Analysis_Module_Summary.xlsx")
write.excel(module_summary, file)

#-------------------------------------------------------------------------------
#' ## Examine module GO enrichment.
#-------------------------------------------------------------------------------

# Singple plot.
# Specify the top percent of terms to print with topN. 
ggplotGOscatter(results_GOenrichment,color = "darkred", topN = 1)

# All plots. 
colors <- as.list(net$colors)
plot_list <- lapply(colors,function(x) ggplotGOscatter(results_GOenrichment,x,topN=0.5))
names(plot_list) <- colors
goplots <- plot_list

# Examine a plot.
plot_list$lightyellow

# Function for visualizing GO terms.
ggplotGOscatter <- function(results_GOenrichment,color,topN=1.0){
  # Collect data in df.
  GOres <- results_GOenrichment[[color]]
  x <- GOres$enrichmentRatio
  y <- -log(GOres$pValue)
  FDR <- as.numeric(GOres$Bonferroni)
  nGenes <- GOres$nCommonGenes
  label <- GOres$shortDataSetName
  df <- data.frame(x,y,FDR,nGenes,label)
  df <- df[order(df$FDR),]
  
  # Hide some of the labels.
  df$label[seq(round(topN*nrow(df)),nrow(df))] <- ""
  
  # Generate plot. 
  plot <- ggplot(df,aes(x = x,y = y, colour = FDR, size = nGenes, label=label)) + 
    geom_point() +  geom_text_repel(colour = "black",alpha = 0.85) + 
    scale_colour_gradient(low = color, high = "white") + 
    xlab("Fold Enrichment") +
    ylab("-Log(P-value)") + 
    ggtitle("Go Enrichment") + 
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"), 
      axis.title.y = element_text(color = "black", size = 11, face = "bold"))
  return(plot)
}

# Save top modules as tiff.
idx <- module_summary$Module[module_summary$p.adj<0.1]
for (i in 1:length(idx)){
  color = idx[i]
  file <- paste0(outputfigsdir,"/",outputMatName,"GO_Scatter_",color,".tiff")
  ggsave(file,plot_list[[color]])
}

#-------------------------------------------------------------------------------
# EBD::Generate Global Network Plots.
#-------------------------------------------------------------------------------

# Insure tissue is specified correctly. 
tissue <- c("Cortex","Striatum")[type]

# Name for output.
FileBaseName <- paste0(outputMatName)

# Initialize PDF
CairoPDF(file = paste0(outputfigsdir,"/",FileBaseName,"_GlobalNetPlots",".pdf"), 
         width = 16, height = 12)

## Plot dendrogram with module colors and trait correlations
MEs <- tmpMEs <- data.frame()
MEList <- moduleEigengenes(t(cleanDat), colors = net$colors)
MEs <- orderMEs(MEList$eigengenes)
# let's be consistent in case prefix was added, remove it.
colnames(MEs) <- gsub("ME", "", colnames(MEs)) 
rownames(MEs) <- rownames(numericMeta)

# Warnings OK; This determines which traits are numeric and if forced to numeric values, 
# non-NA values do not sum to 0
numericIndices <- unique(c(which(!is.na(apply(numericMeta, 2, function(x) sum(as.numeric(x))))), 
                           which(!(apply(numericMeta, 2, function(x) sum(as.numeric(x), na.rm = T))) == 0)))

# Protein-wise Correlations with traits.
geneSignificance <- cor(sapply(numericMeta[, numericIndices], as.numeric), 
                        t(cleanDat), use = "pairwise.complete.obs")
rownames(geneSignificance) <- colnames(numericMeta)[numericIndices]
geneSigColors <- t(numbers2colors(t(geneSignificance), signed = TRUE, 
                                  lim = c(-1, 1), naColor = "black"))
rownames(geneSigColors) <- colnames(numericMeta)[numericIndices]

plotDendroAndColors(
  dendro = net$dendrograms[[1]],
  colors = t(rbind(net$colors, geneSigColors)),
  cex.dendroLabels = 1.2, addGuide = TRUE,
  dendroLabels = FALSE,
  groupLabels = c("Module Colors", colnames(numericMeta)[numericIndices])
)

## Plot eigengene dendrogram/heatmap - using bicor.
tmpMEs <- MEs # net$MEs
colnames(tmpMEs) <- paste("ME", colnames(MEs), sep = "")
MEs[, "grey"] <- NULL
tmpMEs[, "MEgrey"] <- NULL

# Relationships between eigenProteins.
plotEigengeneNetworks(tmpMEs, "Eigengene Network", 
                      marHeatmap = c(3, 4, 2, 2), 
                      marDendro = c(0, 4, 2, 0), 
                      plotDendrograms = TRUE, 
                      xLabelsAngle = 90, 
                      heatmapColors = blueWhiteRed(50))

## Evaluate ME association with groups (numericMeta) using a linear model or ANOVA.
pvec <- lm1 <- regvars <- list()

# First for ANOVA over all Sample.Model combinations/groups
regvars[["overall"]] <- data.frame(as.factor(numericMeta[, "Sample.Model"]), 
                                   as.numeric(numericMeta[, "Age"]), 
                                   as.factor(numericMeta[, "Sex"]))
## data frame with covariates in case we want to try multivariate regression.
colnames(regvars[["overall"]]) <- c("Group", "Age", "Sex") 
lm1[["overall"]] <- lm(data.matrix(MEs) ~ Group, data = regvars[["overall"]])
pvec[["overall"]] <- rep(NA, ncol(MEs))

# Extract F-statistic and model p-value.
for (i in 1:ncol(MEs)) {
  f <- summary(lm1[["overall"]])[[i]]$fstatistic
  pvec[["overall"]][i] <- pf(f[1], f[2], f[3], lower.tail = F)
}
names(pvec[["overall"]]) <- colnames(MEs)

# Then ANOVA on subsets of samples, if appropriate # Region is genetic background!
# Groupings: "Syngap1_KO" "Ube3a_KO"   "Shank2_KO"  "Shank3_KO" vs else????
for (region in unique(numericMeta$Grouping)){
  regvars[[region]] <- data.frame(as.factor(numericMeta[, "Sample.Model"]), 
                                  as.numeric(numericMeta[, "Age"]), 
                                  as.factor(numericMeta[, "Sex"]))[which(numericMeta$Grouping == region), ]
  colnames(regvars[[region]]) <- c("Group", "Age", "Sex")
  ## ANOVA framework yields same results:
  # aov1 <- aov(data.matrix(MEs[which(numericMeta$Group==region),])~Group,data=regvars) 
  lm1[[region]] <- lm(data.matrix(MEs[which(numericMeta$Grouping == region), ]) ~ numericMeta$Sample.Model[which(numericMeta$Grouping==region)],
                      data = regvars[[region]])
  
  pvec[[region]] <- rep(NA, ncol(MEs))
  # Extract F-statistic and p-value corresponding to the  model. 
  for (i in 1:ncol(MEs)) {
    f <- summary(lm1[[region]])[[i]]$fstatistic 
    pvec[[region]][i] <- pf(f[1], f[2], f[3], lower.tail = F)
  }
  names(pvec[[region]]) <- colnames(MEs)
}

## Get sigend kME values
# kME- connectivity module eigengenes. Aka module membership. 
# The protein-wise correlation with module Eigenproteins. 
# tmpMEs is Module EigenProteins with the gray module removed. 
kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc = "bicor")

######################
## Plot eigengene-trait correlations - using p value of bicor for heatmap scale
# bicorAndPvalues seems really similar to Lm/Anova...

library(RColorBrewer)
MEcors <- bicorAndPvalue(MEs, numericMeta[, numericIndices])
# Warning is okay. 
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p

# Which modules exhibit significant correlation with traits?
foo <- melt(moduleTraitCor)
colnames(foo) <- c("Module","Trait","Cor")
man <- melt(moduleTraitPvalue)
colnames(man) <- c("Module","Trait","Pvalue")
modulesTraits <- merge(foo,man,by="Module")


#write.csv(moduleTraitPvalue,"moduleTraitPvalue.csv")
#write.csv(module_pvec,"module_pvec.csv")

textMatrix <- apply(moduleTraitCor, 2, function(x) signif(x, 2))
par(mfrow = c(1, 1))
par(mar = c(6, 8.5, 3, 3))
## Display the correlation values within a heatmap plot
cexy <- if (nModules > 75) {
  0.8
} else {
  1
}
colvec <- rep("white", 1500)
colvec[1:500] <- colorRampPalette(rev(brewer.pal(8, "BuPu")[2:8]))(500)
colvec[501:1000] <- colorRampPalette(c("white", brewer.pal(8, "BuPu")[2]))(3)[2] # interpolated color for 0.05-0.1 p
labeledHeatmap(
  Matrix = apply(moduleTraitPvalue, 2, as.numeric),
  xLabels = colnames(numericMeta)[numericIndices],
  yLabels = paste0("ME", names(MEs)),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = colvec,
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  cex.lab.y = cexy,
  zlim = c(0, 0.15),
  main = paste("Module-trait relationships\n bicor r-value shown as text\nHeatmap scale: Student correlation p value"),
  cex.main = 0.8
)


######################
## Plot eigengene-trait heatmap custom - using bicor color scale

numericMetaCustom <- numericMeta[, numericIndices]
colnames(numericMetaCustom)[which(colnames(numericMetaCustom) == "Treatment.Grouping")] <- "Treatment (P=0 M=0.5 F=3)"
MEcors <- bicorAndPvalue(MEs, numericMetaCustom)
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p

moduleTraitPvalue <- signif(moduleTraitPvalue, 1)
moduleTraitPvalue[moduleTraitPvalue > as.numeric(0.05)] <- as.character("")

textMatrix <- moduleTraitPvalue
# paste(signif(moduleTraitCor, 2), " / (", moduleTraitPvalue, ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
# textMatrix = gsub("()", "", textMatrix,fixed=TRUE)

labelMat <- matrix(nrow = (length(names(MEs))), ncol = 2, data = c(rep(1:(length(names(MEs)))), labels2colors(1:(length(names(MEs))))))
labelMat <- labelMat[match(names(MEs), labelMat[, 2]), ]
for (i in 1:(length(names(MEs)))) {
  labelMat[i, 1] <- paste("M", labelMat[i, 1], sep = "")
}
for (i in 1:length(names(MEs))) {
  labelMat[i, 2] <- paste("ME", labelMat[i, 2], sep = "")
}

# rowMin(moduleTraitPvalue) # if we want to resort rows by min P value in the row
xlabAngle <- if (nModules > 75) {
  90
} else {
  45
}

par(mar = c(16, 12, 3, 3))
par(mfrow = c(1, 1))

bw <- colorRampPalette(c("#0058CC", "white"))
wr <- colorRampPalette(c("white", "#CC3300"))

colvec <- c(bw(50), wr(50))

labeledHeatmap(
  Matrix = t(moduleTraitCor)[, ],
  yLabels = colnames(numericMetaCustom),
  xLabels = labelMat[, 2],
  xSymbols = labelMat[, 1],
  xColorLabels = TRUE,
  colors = colvec,
  textMatrix = t(textMatrix)[, ],
  setStdMargins = FALSE,
  cex.text = 0.5,
  cex.lab.x = cexy,
  xLabelsAngle = xlabAngle,
  verticalSeparator.x = c(rep(c(1:length(colnames(MEs))), as.numeric(ncol(MEs)))),
  verticalSeparator.col = 1,
  verticalSeparator.lty = 1,
  verticalSeparator.lwd = 1,
  verticalSeparator.ext = 0,
  horizontalSeparator.y = c(rep(c(1:ncol(numericMetaCustom)), ncol(numericMetaCustom))),
  horizontalSeparator.col = 1,
  horizontalSeparator.lty = 1,
  horizontalSeparator.lwd = 1,
  horizontalSeparator.ext = 0,
  zlim = c(-1, 1),
  main = "Module-trait Relationships\n Heatmap scale: signed bicor r-value", # \n (Signif. p-values shown as text)"),
  cex.main = 0.8
)

## Plot annotated heatmap - annotate all the metadata, plot the eigengenes!
toplot <- MEs
Grouping <- numericMeta$Grouping
Grouping[numericMeta$SampleType=="WT"] <- "WT"
Gender <- numericMeta$Sex
Gender[Gender == 0] <- "Female"
Gender[Gender == 1] <- "Male"
Gender <- factor(Gender)

metdat <- data.frame(WT_vs_other = Grouping, Age = numericMeta[, "Age"], 
                     Gender = Gender, Genotype = numericMeta[, "Model"])

colnames(toplot) <- colnames(MEs)
rownames(toplot) <- rownames(MEs)
toplot <- t(toplot) # have to transpose

for (pvector in names(pvec)) {
  pvec[[pvector]] <- pvec[[pvector]][match(names(pvec[[pvector]]), rownames(toplot))]
}
rownames(toplot) <- paste(orderedModules[match(colnames(MEs), orderedModules[, 2]), 1], " ", 
                          rownames(toplot), "  |  K-W P.overall=", signif(pvec[["overall"]], 2), sep = "")

## Plot heatmaps of eigenproteins
treatmentColorVec <- statusColorVec <- vector()
for (Genotype in 1:length(unique(numericMeta$Model))){
  treatmentColorVec <- c(treatmentColorVec, 
                         numericMeta$Model[which(numericMeta$Model == levels(factor(numericMeta$Model))[Genotype])[Genotype]])
}

treatmentColorVec[treatmentColorVec == " Shank2"] <- "yellow"
treatmentColorVec[treatmentColorVec == " Shank3"] <- "blue"
treatmentColorVec[treatmentColorVec == " Syngap1"] <- "green"
treatmentColorVec[treatmentColorVec == " Ube3a"] <- "purple"

statusColorVec <- levels(as.factor(unique(numericMeta$Color)))

heatmapLegendColors <- list(
  "WT_vs_other" = statusColorVec,
  "Age" = c("white", "seagreen3"),
  "Gender" = c("pink", "dodgerblue"),
  "Genotype" = treatmentColorVec,
  "Modules" = sort(colnames(MEs))
)

# Requires NMF
par(mfrow = c(1, 1))
aheatmap(
  x = toplot, ## Numeric Matrix
  main = "Plot of Eigengene-Trait Relationships - SAMPLES IN BATCH (REGION) ORDER",
  annCol = metdat,
  annRow = data.frame(Modules = colnames(MEs)),
  annColors = heatmapLegendColors,
  border = list(matrix = TRUE),
  scale = "row",
  distfun = "correlation", hclustfun = "average", ## Clustering options
  cexRow = 0.8, ## Character sizes
  cexCol = 0.8,
  col = blueWhiteRed(100), ## Color map scheme
  treeheight = 80,
  Rowv = TRUE, Colv = NA
) ## Do not cluster columns - keep given order


aheatmap(
  x = toplot, ## Numeric Matrix
  main = "Plot of Eigengene-Trait Relationships - SAMPLES CLUSTERED",
  annCol = metdat,
  annRow = data.frame(Modules = colnames(MEs)),
  annColors = heatmapLegendColors,
  border = list(matrix = TRUE),
  scale = "row",
  distfun = "correlation", hclustfun = "average", ## Clustering options
  cexRow = 0.8, ## Character sizes
  cexCol = 0.8,
  col = blueWhiteRed(100), ## Color map scheme
  treeheight = 80,
  Rowv = TRUE, Colv = TRUE
) ## Cluster columns


dev.off()
## Close PDF 1.

#-------------------------------------------------------------------------------
## Extract module significance data.
#-------------------------------------------------------------------------------
# Get Module p-values (p-value for ME lm fit to numericMeta groupings)
module_pvec <- t(do.call(rbind,pvec))
# Rename columns
colnames(module_pvec) <- paste(colnames(module_pvec),"p_value")
head(module_pvec)
dim(module_pvec)

# Number of significant associations.
print(paste(table(module_pvec[,1]<0.05)[2],"of", nModules, 
            "modules exhibit a significant association to traits."))
print(paste(table(module_pvec[,2:ncol(module_pvec)]<0.05)[2],
            "instances of significant association to traits overall."))

#-------------------------------------------------------------------------------
## Output PDF of verboseScatterPlots and boxplots for MEs factored into your groups of interest...
#-------------------------------------------------------------------------------
CairoPDF(file = paste0(outputfigs, "/2.GlobalNetPlots(BoxPlots)_", FileBaseName, ".pdf"), 
         width = 18, height = 11.25)

# control where plots created in the below for loop will be displayed [out of order from the order of creation in the inner for (region... loop]
layout(
  mat = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20), nrow = 4, ncol = 5, byrow = TRUE),
  heights = c(0.9, 0.9, 0.9, 0.9), # Heights of the four rows
  widths = c(1.75, 0.9, 0.9, 0.9, 0.9)
) # Widths of the 5 columns
# sapply(c(1:20),layout.show)


for (i in 1:(nrow(toplot))) { # grey already excluded when setting MEs data frame, -1 unnecessary
  boxplot(toplot[i, ] ~ factor(numericMeta$Sample.Model, 
                               levels = rev(levels(factor(numericMeta$Sample.Model)))), 
          col = colnames(MEs)[i], ylab = "Eigenprotein Value", 
          main = paste0(gsub("K-W", "Kruskal-Wallis", rownames(toplot)[i])), xlab = NULL, las = 2, outline = FALSE)
  for (region in unique(numericMeta$Group)) {
    titlecolor <- if (signif(pvec[[region]], 2)[i] < 0.05) {
      "red"
    } else {
      "black"
    }
    boxplot(toplot[i, which(numericMeta$Group == region)] ~ factor(numericMeta$Sample.Model[which(numericMeta$Group == region)], 
                                                                   levels = rev(levels(factor(numericMeta$Sample.Model[which(numericMeta$Group == region)])))), 
            col = colnames(MEs)[i], ylab = "Eigenprotein Value", main = paste0(orderedModules[match(colnames(MEs)[i], orderedModules[, 2]), 1], "\nK-W P.", region, " = ", signif(pvec[[region]], 2)[i]), xlab = NULL, las = 2, col.main = titlecolor, outline = TRUE) # outline=T, show outliers
  }
}

dev.off()
## Close PDF 2.