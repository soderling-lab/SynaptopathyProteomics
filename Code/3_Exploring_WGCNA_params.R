#' ---
#' title: Exploring WGCNA Parameters. 
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
#' 
#-------------------------------------------------------------------------------
#' ## Prepare the workspace.
#-------------------------------------------------------------------------------
#+ eval = TRUE, echo = FALSE, error = FALSE

# Use ctl+alt+T to execute a code chunk.

# Run this chunk before doing anything!
rm(list = ls())
dev.off()
cat("\014") # alternative is cat("\f")
options(stringsAsFactors = FALSE)

#  Load required packages.
suppressPackageStartupMessages({
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

# To install TBmiscr:
#library(devtools)
#devtools::install_github("twesleyb/TBmiscr")

# Define version of the code.
CodeVersion <- "Exploring_Params"

# Define tisue type: cortex = 1; striatum = 2.
type <- 3
tissue <- c("Cortex", "Striatum", "Combined")[type]

# Set the working directory.
rootdir <- "D:/Documents/R/Synaptopathy-Proteomics"
#rootdir <- "C:/Users/User/Documents/Tyler Bradshaw/Synaptopathy-Proteomics"
setwd(rootdir)

# Set any other directories.
functiondir <- paste(rootdir, "Code", sep = "/")
datadir <- paste(rootdir, "Input", sep = "/")
Rdatadir <- paste(rootdir,"RData", sep = "/")

# Load required custom functions.
my_functions <- paste(functiondir, "0_TMT_Preprocess_Functions.R", sep = "/")
source(my_functions)

# Define prefix for output figures and tables.
outputMatName <- paste(tissue, "_WGCNA_Analysis_", sep = "")

# Globally set ggplots theme.
ggplot2::theme_set(theme_gray())

#-------------------------------------------------------------------------------
#' ## Start WGCNA. Choosing a soft thresholding power, Beta.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Estimate powers?
estimatePower <- TRUE

# Data is...
# Load TAMPOR cleanDat from file: #2918 of 2918
datafile <- paste(Rdatadir,tissue,"TAMPOR_data_outliersRemoved.Rds",sep="/")
cleanDat <- readRDS(datafile)
cleanDat <- log2(cleanDat)
cleanDat[1:5,1:5] # Data should be log transformed. 
dim(cleanDat)

# Load combined sample info.
traitsfile <- paste(Rdatadir,tissue,"Combined_Cortex_Striatum_traits.Rds",sep="/")
sample_info <- readRDS(traitsfile)
sample_info[1:5,1:5]
dim(sample_info)

# Allow parallel WGCNA calculations:
allowWGCNAThreads()
parallelThreads <- 11
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
  #file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology_Table.tiff")
  #ggsave(file,table)
  
  # Figure. ggplotScaleFreeFit() generates three plots.
  plots <- ggplotScaleFreeFit(sft)
  plots$Grid
  
  # Save as tiff.
  #file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology.tiff")
  #ggsave(file,plots$Grid)
  
  # Save plots and table as PDF.
  #plot_list <- list(table,plots$ScaleFreeFit, plots$MeanConnectivity)
  #file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology.pdf")
  #ggsavePDF(plot_list,file)
}

#-------------------------------------------------------------------------------
#' ## Prepare to sample WGCNA parameters.
#-------------------------------------------------------------------------------

# Allow parallel WGCNA calculations:
allowWGCNAThreads()
parallelThreads <- 11
clusterLocal <- makeCluster(c(rep("localhost", parallelThreads)), type = "SOCK")
registerDoParallel(clusterLocal)

# Main network building parameters. 
power <- 12
corType <- "bicor"
networkType <- "signed"

# Other key parameters for optimization:
minModSize <- seq(3,50,by=1)           # Minimum module size. seq(1,50,by=1) 
deepSplit <- seq(0,4,by=1)             # Sensitivity for module splitting [0-4]. 4 is most sensitive. Increasing results in more modules.  
mergeCutHeight <- seq(0.01,0.2,0.01)      # Cut height for module detection. Was 0.07. Increasing results in more modules.  
reassignThresh <- seq(0.01,0.1,0.01)     # pvalue threshold for reassigning nodes to modules. 
minKMEtoStay <- seq(0.1,0.7,0.1)        # minimum module connectivity score for assigning to a module. 
minCoreKMESize <- seq(3,15,1)
pamStage <- TRUE
maxBlockSize <- 12000
detectCutHeight <- 0.995

# Calculate the adjacency network.
r <- bicor(t(cleanDat))
adjm <- ((1+r)/2)^power #signed network.
#adjm <- abs(r)^power     #un-signed.

# Create igraph object. 
graph <- graph_from_adjacency_matrix(
  adjmatrix = adjm, 
  mode = c("undirected"), 
  weighted = TRUE, 
  diag = FALSE)

# Combine into a list. 
params_list <- list(minModSize,deepSplit,mergeCutHeight,reassignThresh,minKMEtoStay,minCoreKMESize,pamStage)
names(params_list) <- c("minModSize","deepSplit","mergeCutHeight","reassignThresh","minKMEtoStay","minCoreKMESize","pamStage")

# Use epand.grid() to generate matrix of all posible combinations of params.
params_grid <- expand.grid(params_list)

# Sample parameters, nboot iterations.
nboot <- 1000
rand_params <- sample(1:nrow(params_grid),nboot)
out <- list()

#-------------------------------------------------------------------------------
#' ## Loop to sample parameters.
#-------------------------------------------------------------------------------

for (i in 1:nboot){
  print(paste("Sampling network building parameters, iteration",i))
  
  # Random sampling of params.
  idx <- rand_params[i]
  params <- params_grid[idx,]
  params$idx <- idx
  
  # Call blockwiseModules to build WGCNA network. 
  # Setting saveTOM = FALSE will really slow thing down.
  net <- blockwiseModules(t(cleanDat),
                          power = power, 
                          deepSplit = params$deepSplit, 
                          minModuleSize = params$minModSize,
                          mergeCutHeight = params$mergeCutHeight, 
                          TOMDenom = "mean",
                          detectCutHeight=detectCutHeight, 
                          corType = corType,
                          networkType = networkType, 
                          pamStage = params$pamStage,
                          pamRespectsDendro = TRUE, 
                          reassignThresh = params$reassignThresh, 
                          minCoreKMESize = params$minCoreKMESize, 
                          minKMEtoStay = params$minKMEtoStay,
                          verbose = 0, 
                          saveTOMs = FALSE, 
                          maxBlockSize = maxBlockSize)
  
  ## Enforce module preservation. .
  # Input for NetRep:
  data_list <- list(data = t(cleanDat))
  correlation_list <- list(data = r)       
  network_list <- list(data = adjm)       
  module_labels <- net$colors     
  names(module_labels) <- rownames(cleanDat)
  
  # Progress report. 
  print(paste("Performing permutation tests to calculate module preservation, iter",i))
  
  # Calculate module preservation statistics.
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
    nThreads = parallelThreads, 
    #nPerm = 1000, # nPerm will be determined by the function. 
    null = "overlap", 
    alternative = "greater", 
    simplify = TRUE,
    verbose = FALSE)
  
  # Extract preservation statistics and p-values. 
  preservation <- preservation[c("observed","p.values")]
  
  # Get the maximum permutation test p-value.
  maxPval <- apply(preservation$p.values, 1, function(x) max(x,na.rm=TRUE))
  nModules <- length(unique(net$colors))-1
  params$nModules <- nModules
  
  # Modules removed if adjusted pvalue is greater than alpha.
  alpha = 0.05
  modules_out <- names(maxPval)[maxPval>alpha/nModules]
  #names(net$MEs)[grepl(paste(modules_out,collapse="|"),names(net$MEs))]<-"MEgrey"
  
  # Drop nsig Modules (Set color to grey)
  net$colors[net$colors %in% modules_out] <- "grey"
  nModules <- length(unique(net$colors))-1
  
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
  
  # Fraction of un-assigned proteins. 
  params$PercentGrayNodes <- 100*table(net$colors)["grey"]/length(net$goodGenes)
  params$nSigModules <- nModules
  
  # Rogdi with Wdr7?
  params$Rogdi.Wdr7 <- net$colors[rownames(cleanDat) == "Rogdi|Q3TDK6"]==net$colors[rownames(cleanDat) == "Wdr7|Q920I9"]
  
  # PSMD Complex
  idx <- match(rownames(cleanDat)[grepl("Psm",rownames(cleanDat))],rownames(cleanDat))
  tab <- as.data.frame(table(net$colors[idx]))
  row <- grepl("grey",tab$Var1)
  tab <- tab[!row,]
  if (dim(tab)[1] ==0){
    params$MaxPSMD <- 0
  }else{
    params$MaxPSMD <- max(tab$Freq)
  }
  
  # Calculate modularity, q.
  membership <- as.numeric(as.factor(net$colors))
  q1 <- modularity(graph, membership, weights = edge_attr(graph, "weight"))
  params$q1 <- q1
  
  # Without "grey" nodes.
  v <- rownames(cleanDat)[!net$colors=="grey"]
  subg <- induced_subgraph(graph,v)
  membership <- as.numeric(as.factor(net$colors))
  membership <- membership[!net$colors=="grey"]
  q2 <- modularity(subg, membership, weights = edge_attr(subg, "weight"))
  params$q2 <- q2
  
  # Module coherence (aka PVE)
  modc <- as.matrix(
    propVarExplained(datExpr=t(cleanDat), colors=net$colors, MEs=net$MEs)
  )
  modc <- modc[!(rownames(modc) =="PVEgrey"),]
  
  # Insure numericMeta rows and cleanDat columns match
  numericMeta <- sample_info
  numericMeta <- numericMeta[match(colnames(cleanDat), rownames(numericMeta)), ]
  
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
  
  # Determine numerical indices. The columns of numericMeta with numerical data.
  # Warnings OK; This determines which traits are numeric and if forced to numeric values, 
  # non-NA values do not sum to 0.
  numericIndices <- unique(c(which(!is.na(apply(numericMeta, 2, function(x) sum(as.numeric(x))))), 
                             which(!(apply(numericMeta, 2, function(x) sum(as.numeric(x), na.rm = T))) == 0)))
  ## Calc MEs
  # Calculate Module EigenProteins (MEs)
  MEs <- tmpMEs <- data.frame()
  MEList <- moduleEigengenes(t(cleanDat), colors = net$colors)
  MEs <- orderMEs(MEList$eigengenes)
  
  # Remove prefix.
  colnames(MEs) <- gsub("ME", "", colnames(MEs)) 
  rownames(MEs) <- rownames(numericMeta)
  
  ## Determine the cost of removing grey nodes.
  geneSignificance <- abs(
    cor(sapply(numericMeta[, numericIndices], as.numeric), 
        t(cleanDat), use = "pairwise.complete.obs"))
  #sum(geneSignificance)
  
  # Add rownames. 
  rownames(geneSignificance) <- colnames(numericMeta)[numericIndices]
  idx <- colnames(geneSignificance) %in% rownames(cleanDat)[net$colors=="grey"]
  GSsub <- geneSignificance[,idx]
  costGrey <- sum(GSsub)/sum(geneSignificance)
  
  ## Calc ME cor Traits.
  
  # Warnings are okay.
  MEcorTraits <- bicorAndPvalue(MEs, numericMeta[, numericIndices])
  
  # Convergence... a module that is highly correlated with several traits.
  # Also, the degree of a module node in the module-trait network.
  MSdm <- MEcorTraits$bicor
  
  # Get comparisons of interest. 
  cols <- c(
    "Cortex.KO.Shank2","Striatum.KO.Shank2",
    "Cortex.KO.Shank3","Striatum.KO.Shank3",
    "Cortex.HET.Syngap1","Striatum.HET.Syngap1",
    "Cortex.KO.Ube3a","Striatum.KO.Ube3a")
  idx <- match(cols,colnames(MSdm))
  dm <- MSdm[,idx]
  
  MSdm <- as.data.frame(dm)
  MSdm$Degree <- rowSums(abs(MSdm))
  
  params$MSmax <- max(MSdm$Degree)
  params$MStotal <- sum(MSdm$Degree)
  params$MSmean <- mean(MSdm$Degree)
  params$MSmedian <- median(MSdm$Degree)
  params$MSgray <-MSdm$Degree[rownames(MSdm)=="grey"]
  params$costGrey <- costGrey
  params$medianModCoherence <- median(modc)
  
  ##KW stats.
  # Insure that ME data is correct format.
  traits <- numericMeta
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
  
  # Pool WT.Cortex and WT.Striatum.
  groups <- paste(numericMeta$TissueType,numericMeta$Sample.Model,sep=".")
  groups[grepl("Cortex.WT",groups)] <- "WT.Cortex"
  groups[grepl("Striatum.WT",groups)] <- "WT.Striatum"
  
  # Perform KW test. 
  KWtest <- apply(MEs,2,function(x) kruskal.test(x,as.factor(groups)))
  
  # Extract pvalues from the list of KW tests. 
  # The pvalue is the 3rd element of each list.
  KW_results <- as.data.frame(do.call(rbind, sapply(KWtest,"[",3)))
  
  # Clean up the result. Drop Grey!
  colnames(KW_results) <- "p.value"
  rownames(KW_results) <- gsub(".p.value","",rownames(KW_results))
  KW_results$p.adj <- NA
  KW_results <- KW_results[!rownames(KW_results)=="grey",]
  KW_results$p.adj <- p.adjust(KW_results$p.value,method="bonferroni")
  
  # Number of significant modules. 
  params$KWsigModules <- sum(KW_results$p.adj<0.05)
  sigModules <- rownames(KW_results)[KW_results$p.adj<0.05]
  
  # Post-hoc comparisons with DunnettTest (comparison to control)
  
  #x = as.numeric(MEs[[1]])
  g <- as.factor(groups)
  cox_subset <- grepl("Cortex",g)
  str_subset <- grepl("Striatum",g)
  #DunnettTest(x[cox_subset],g[cox_subset],control="WT.Cortex")
  #DunnettTest(x[str_subset],g[str_subset],control="WT.Striatum")
  
  list1 <- lapply(MEs,function(x) 
    DunnettTest(as.numeric(x)[cox_subset], as.factor(groups)[cox_subset], control = "WT.Cortex"))
  list1 <- sapply(list1,"[",1)
  
  list2 <- lapply(MEs,function(x) 
    DunnettTest(as.numeric(x)[str_subset], as.factor(groups)[str_subset], control = "WT.Striatum"))
  list2 <- sapply(list2,"[",1)
  
  new_list <- list()
  for (k in 1:length(list1)){
    new_list[[k]] <- rbind(list1[[k]],list2[[k]])
  }
  names(new_list) <- names(MEs)
  new_list <- lapply(new_list,function(x) as.data.frame(x))
  nsig <- sum(unlist(lapply(new_list,function(x) sum(x$pval<0.05))))
  params$DTsigTests <- nsig

  # Add p.adj.
  new_list <- lapply(new_list,function(x) add_column(x, p.adj = p.adjust(x$pval,method="bonferroni")))
  
  # Dunnett sig tests for KW sig modules.
  dm <- do.call(rbind,lapply(new_list,function(x) sum(x$p.adj<0.05)))
  idx <- match(rownames(KW_results),rownames(dm))
  KW_results$nSigDunnett <- dm[idx,]
  params$nModConvergence <- sum(KW_results$nSigDunnett>1)
  
  # Store in list. 
  out[[i]] <- list(NetworkParams = params,
                   ModuleColors = net$colors, 
                   ModulePreservation = preservation)
  
  # Save output every 10 iters.
  if (any((i==seq(1,nboot,10)))){
    collectGarbage()
    print(paste("Saving progress at:",Sys.time()))
    # Save results to file.
    file <- paste0(Rdatadir,"/","Sample_blockwiseModules_Stats.RDS")
    saveRDS(out,file)
    # Push to GitHub every 100.
    if (any((i==seq(1,nboot,100)))){
      print("Updating GitHub with progress.")
      gitstatus()
      gitadd()
      gitcommit(msg= paste("Exploring params iter",i,"complete."))
      gitpush()
    }
  }
}

## END LOOP.

#-------------------------------------------------------------------------------
# Load the results.
#-------------------------------------------------------------------------------

# Load Parameters and network stats.
files <- list.files(Rdatadir, pattern = "Stats")
file <- paste0(Rdatadir,"/",files)
file
out <- readRDS(file)
length(out)

# Load Dunnett tests.
#files <- list.files(Rdatadir)
#files
#file <- paste0(Rdatadir,"/",files[10])
#Dtests <- readRDS(file)
#length(Dtests)

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

# Rank stats.
result$coherenceRank <- rank(result$medianModCoherence)/nrow(result)
result$costRank <- rank(result$costGrey)/nrow(result)
result$modRank <- rank(result$q2)/nrow(result)
result$nModRank <- rank(result$nSigModules)/nrow(result)
result$nSigModRank <- rank(result$KWsigModules)/nrow(result)
result$maxPSMDRank <- rank(result$MaxPSMD)/nrow(result)
result$score <- (result$coherenceRank + result$modRank)/2 - result$costRank
#result$score <- (result$coherenceRank + result$nSigModRank + result$modRank)/3 - result$costRank
#result$score2 <- result$maxPSMDRank + result$nSigModRank - result$costRank

# Examine Dtests results.
#result$nSigDunnettTests <- do.call(rbind,
#                                   lapply(Dtests,function(x) sum(x$nSigDunnettTest)))
# If KW pval is <0.05, then number of sig Dunnett tests.
#Dtests <- lapply(Dtests,function(x) 
#  add_column(x,logic = x$nSigDunnettTest*as.numeric(x$p.adj<0.05),.after=5))
# Convergence, the number of modules with more than 2 sig changes.
#result$nModsConvergance <- do.call(rbind,
#                                   lapply(Dtests,function(x) sum(x$logic>1)))
# Write to file. 
#write.csv(result,"params_result.csv")

#-------------------------------------------------------------------------------
#' ## Loop to recalculate Dunnett test results with p.adj, exploring convergance.
#-------------------------------------------------------------------------------

# CleanDat
file <- paste(Rdatadir,tissue,"TAMPOR_data_outliersRemoved.RDS",sep="/")
cleanDat <- readRDS(file)
cleanDat <- log2(cleanDat)
dim(cleanDat)

# Sample info.
file <- paste(Rdatadir,tissue,"sample_info.RDS",sep="/")
sample_info <- numericMeta <- readRDS(file)

# Load Sampled Network parameters and stats.
files <- list.files(Rdatadir)
file <- paste0(Rdatadir,"/",files[17])
out <- readRDS(file)

# Network stats
result <- do.call(rbind,sapply(out,"[",1))
result$iter <- c(1:nrow(result))
rownames(result) <- paste0("params_",1:nrow(result))

# Module colors.
modcolors <- sapply(out,"[",2)
names(modcolors) <- rownames(result)

#### LOOP to perform Dunnett test.
output <- list()
n_iter <- 1

for (k in 94){
  print(paste0("iteration ",k,"..."))
  # Calculate Module EigenProteins.
  MEs <- tmpMEs <- data.frame()
  MEList <- moduleEigengenes(t(cleanDat), colors = modcolors[[k]])
  MEs <- orderMEs(MEList$eigengenes)
  # let's be consistent in case prefix was added, remove it.
  colnames(MEs) <- gsub("ME", "", colnames(MEs)) 
  rownames(MEs) <- colnames(cleanDat)
  
  # Insure traits are in matching order.
  idx <- match(rownames(MEs),rownames(numericMeta))
  traits <- numericMeta[idx,]
  all(rownames(traits)==rownames(MEs))
  
  # Define groups, the biological groups of interest. 
  groups <- paste(traits$Tissue,traits$Sample.Model,sep=".")
  groups[grepl("Cortex.WT",groups)] <- "WT.Cortex"
  groups[grepl("Striatum.WT",groups)] <- "WT.Striatum"
  #unique(groups)
  
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
  #print(paste("nModules with p.adj < 0.05 =", length(KWsigModules)))
  #print(paste("nModules with FDR < 0.05 =", length(KW_results$Module[KW_results$FDR<0.05])))
  
  # Split Module EigenProtein (MEs) dm into a list of column vectors for lapply. 
  ME_list <- split(as.matrix(MEs), rep(1:ncol(MEs), each = nrow(MEs)))
  names(ME_list) <- colnames(MEs)
  
  # Add vector of groups
  ME_list <- lapply(ME_list,function(x) data.frame(x = x, groups = groups))
  
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
  
  # Loop through ME_list and generate verboseBoxPlot.
  # lapply wont work here because the name is not preserved when you call lapply()...
  # method is the p.adj method for the Dunn's test p-value.
  plot_data <- list()
  for (i in 1:dim(MEs)[2]){
    x <- ME_list[[i]]$x
    g <- ME_list[[i]]$groups
    color <- names(ME_list)[[i]]
    plot <- ggplotVerboseBoxplot(x,g,levels,contrasts,color,stats=TRUE,method="dunnett")
    plot_data[[i]] <- plot
    names(plot_data)[[i]] <- color
  }
  
  plots <- sapply(plot_data,"[",1)
  
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
  
  # Just KW overall significant modules. 
  res <- do.call(rbind,lapply(Dtest_stats,function(x) sum(x$P.adj<0.05)))
  idx <- match(KW_results$Module,rownames(res))
  KW_results$nSigDunnettTest <- res[idx,]
  
  output[[k]] <- KW_results
  
  # Save output.
  if (k==n_iter){
    print(paste("Saving progress at:",Sys.time()))
    # Save results to file.
    file <- paste0(Rdatadir,"/","Params_DunnettTests.RDS")
    saveRDS(output,file)
    gitstatus()
    gitadd()
    gitcommit(msg= paste("Dunnett tests for sampled params"))
    gitpush()
    
  }
}
