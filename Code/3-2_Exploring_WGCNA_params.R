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

# Use ctl+alt+T to execute a code chunk.

# Run this chunk before doing anything!
rm(list = ls())
if (.Device != "null device" ) dev.off()
cat("\f")
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
# library(devtools)
# devtools::install_github("twesleyb/TBmiscr")

# Define version of the code.
CodeVersion <- "Exploring_Params"

# Define tisue type: cortex = 1; striatum = 2.
type <- 3
tissue <- c("Cortex", "Striatum", "Combined")[type]

# Set the working directory.
rootdir <- "D:/projects/Synaptopathy-Proteomics"
rootdir <- "C:/Users/User/Documents/Tyler/Synaptopathy-Proteomics"
setwd(rootdir)

# Set any other directories.
functiondir <- paste(rootdir, "Code", sep = "/")
datadir <- paste(rootdir, "Input", sep = "/")
Rdatadir <- paste(rootdir, "RData", sep = "/")

# Load required custom functions.
my_functions <- paste(functiondir, "0_TMT_Preprocess_Functions.R", sep = "/")
source(my_functions)

# Define prefix for output figures and tables.
outputMatName <- paste(tissue, "_WGCNA_Analysis_", sep = "")

# Globally set ggplots theme.
ggplot2::theme_set(theme_gray())

#-------------------------------------------------------------------------------
#' ## Build a  dictionary-like object for mapping gene identifiers.
#-------------------------------------------------------------------------------

# Load WGCNA network and meta data.
file <- paste(Rdatadir, "Network_and_metaModules.Rds", sep = "/")
data <- readRDS(file)
net <- data$net
meta <- data$meta

# Create a dictionary like object for mapping entrez to gene|uniprot.
entrez2protein <- as.list(meta$protein)
names(entrez2protein) <- meta$entrez

#-------------------------------------------------------------------------------
#' ## Start WGCNA. Choosing a soft thresholding power, Beta.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Estimate powers?
estimatePower <- TRUE

# Data is...
# Load TAMPOR cleanDat from file: #2918 of 2918
datafile <- paste(Rdatadir, tissue, "TAMPOR_data_outliersRemoved.Rds", sep = "/")
cleanDat <- readRDS(datafile)
cleanDat <- log2(cleanDat)
cleanDat[1:5, 1:5] # Data should be log transformed.
dim(cleanDat)

################################################################################
## Run this chunk if subsetting the data based on DEP communities.

# Load the DEP protein communities. 
file <- paste0(Rdatadir,"/","DEP_Communities.RDS")
community_results <- readRDS(file)

# Pick a group/genotype.
n <- 5
group <- c("Shank2", "Shank3", "Syngap1", "Ube3a", "All")[n]

# Define proteins of interst.
if (group == "All"){
  v <- unique(as.vector(unlist(sapply(community_results,"[",4))))
  prots <- unlist(entrez2protein[v])
}else{
  prots <- unlist(entrez2protein[community_results[[group]][[4]]])
}
subg_name <- group
subDat <- subset(cleanDat, rownames(cleanDat) %in% prots)

# Write over cleanDat.
cleanDat <- subDat
dim(cleanDat)

################################################################################
## Run this chunk if splitting data into WT and KO!

# Load combined sample info.
traitsfile <- paste(Rdatadir, tissue, "Combined_Cortex_Striatum_traits.Rds", sep = "/")
sample_info <- readRDS(traitsfile)
sample_info[1:5, 1:5]
dim(sample_info)

wt_samples <- subset(sample_info$SampleID, sample_info$SampleType == "WT")
ko_samples <- subset(sample_info$SampleID, 
                     sample_info$SampleType == "HET" | sample_info$SampleType == "KO")

wt_dat <- cleanDat[,colnames(cleanDat) %in% wt_samples]
ko_dat <- cleanDat[,colnames(cleanDat) %in% ko_samples]
dim(wt_dat)
dim(ko_dat)
subdat <- list(wt = wt_dat,
               ko = ko_dat)

# Allow parallel WGCNA calculations:
allowWGCNAThreads()
parallelThreads <- 8
clusterLocal <- makeCluster(c(rep("localhost", parallelThreads)), type = "SOCK")
registerDoParallel(clusterLocal)

## Determine soft power, beta.
# Vector of powers to test:
powers <- seq(4, 20, by = 1.0)

# Soft Power selection
sft <- lapply(subdat, function(x) {
  pickSoftThreshold(t(x), 
  powerVector = powers, 
  corFnc = "bicor",
  blockSize = 15000, 
  verbose = 3, 
  networkType = "signed")})

# Figure. ggplotScaleFreeFit() generates three plots.
plots <- lapply(sft, function(x) ggplotScaleFreeFit(x))

# Choose minimum power to achieve scale free fit > 0.8.
power_beta <- lapply(sft, function(x) x$fitIndices$Power[x$fitIndices$SFT.R.sq > 0.8][1])
power_beta

################################################################################
## Estimate scale free fit.

# Load combined sample info.
traitsfile <- paste(Rdatadir, tissue, "Combined_Cortex_Striatum_traits.Rds", sep = "/")
sample_info <- readRDS(traitsfile)
sample_info[1:5, 1:5]
dim(sample_info)

# Allow parallel WGCNA calculations:
allowWGCNAThreads()
parallelThreads <- 8
clusterLocal <- makeCluster(c(rep("localhost", parallelThreads)), type = "SOCK")
registerDoParallel(clusterLocal)

## Determine soft power, beta.
# Vector of powers to test:
powers <- seq(4, 20, by = 0.5)

# Soft Power selection
if (estimatePower == TRUE) {
  sft <- pickSoftThreshold(t(cleanDat),
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
  # ggsave(file,table)

  # Figure. ggplotScaleFreeFit() generates three plots.
  plots <- ggplotScaleFreeFit(sft)
  plots$Grid

  # Save as tiff.
  # file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology.tiff")
  # ggsave(file,plots$Grid)

  # Save plots and table as PDF.
  # plot_list <- list(table,plots$ScaleFreeFit, plots$MeanConnectivity)
  # file <- paste0(outputfigsdir,"/",outputMatName,"ScaleFreeTopology.pdf")
  # ggsavePDF(plot_list,file)
}

# Choose minimum power to achieve scale free fit > 0.8.
power <- sft$fitIndices$Power[sft$fitIndices$SFT.R.sq > 0.8][1]
print(paste("Power (beta):",power))

#-------------------------------------------------------------------------------
#' ## Prepare to sample WGCNA parameters.
#-------------------------------------------------------------------------------

################################################################################
## Parameter optimization will be performed with the WT network.
print("Using the WT data for WGCNA parameter optimization!")
power_beta <- power_beta$wt
cleanDat <- subdat$wt
dim(cleanDat)
cleanDat[1:5,1:5] # Data should be log transformed!

################################################################################

# Allow parallel WGCNA calculations:
allowWGCNAThreads()
parallelThreads <- 8
clusterLocal <- makeCluster(c(rep("localhost", parallelThreads)), type = "SOCK")
registerDoParallel(clusterLocal)

# Main network building parameters.
power_beta <- power_beta
corType <- "bicor"
networkType <- "signed"

# Other key parameters for optimization:
minModSize <- seq(3, 50, by = 1) # Minimum module size. seq(1,50,by=1)
deepSplit <- seq(0, 4, by = 1) # Sensitivity for module splitting [0-4]. 4 is most sensitive. Increasing results in more modules.
mergeCutHeight <- seq(0.01, 0.2, 0.01) # Cut height for module detection. Was 0.07. Increasing results in more modules.
reassignThresh <- seq(0.01, 0.1, 0.01) # pvalue threshold for reassigning nodes to modules.
minKMEtoStay <- seq(0.1, 0.7, 0.1) # minimum module connectivity score for assigning to a module.
minCoreKMESize <- seq(3, 15, 1)
pamStage <- TRUE
maxBlockSize <- 12000
detectCutHeight <- 0.995

# Calculate the adjacency network.
r <- bicor(t(cleanDat))
adjm <- ((1 + r) / 2)^power_beta # signed network.

# Create igraph object.
graph <- graph_from_adjacency_matrix(
  adjmatrix = adjm,
  mode = c("undirected"),
  weighted = TRUE,
  diag = FALSE
)

# Combine into a list.
params_list <- list(minModSize, deepSplit, mergeCutHeight, reassignThresh, minKMEtoStay, minCoreKMESize, pamStage)
names(params_list) <- c("minModSize", "deepSplit", "mergeCutHeight", "reassignThresh", "minKMEtoStay", "minCoreKMESize", "pamStage")

# Use epand.grid() to generate matrix of all posible combinations of params.
params_grid <- expand.grid(params_list)

# Sample parameters, nboot iterations.
nboot <- 1000
rand_params <- sample(1:nrow(params_grid), nboot)
out <- list()

# File for saving output.
file <- paste0(Rdatadir, "/", "Sample_blockwiseModules_WT_Only_Stats.RDS")

# Should progress be pushed to git?
git_push <- FALSE

# Clean up.
collectGarbage()

#-------------------------------------------------------------------------------
#' ## Loop to sample parameters.
#-------------------------------------------------------------------------------
for (i in 1:nboot) {
  print(paste("Sampling network building parameters, iteration", i))

  # Random sampling of params.
  idx <- rand_params[i]
  params <- params_grid[idx, ]
  params$idx <- idx

  # Call blockwiseModules to build WGCNA network.
  # Setting saveTOM = FALSE will really slow thing down.
  net <- blockwiseModules(t(cleanDat),
    power = power_beta,
    deepSplit = params$deepSplit,
    minModuleSize = params$minModSize,
    mergeCutHeight = params$mergeCutHeight,
    TOMDenom = "mean",
    detectCutHeight = detectCutHeight,
    corType = corType,
    networkType = networkType,
    pamStage = params$pamStage,
    pamRespectsDendro = TRUE,
    reassignThresh = params$reassignThresh,
    minCoreKMESize = params$minCoreKMESize,
    minKMEtoStay = params$minKMEtoStay,
    verbose = 0,
    saveTOMs = FALSE,
    maxBlockSize = maxBlockSize
  )

  ## Enforce module preservation. .
  # Input for NetRep:
  data_list <- list(data = t(cleanDat))
  correlation_list <- list(data = r)
  network_list <- list(data = adjm)
  module_labels <- net$colors
  names(module_labels) <- rownames(cleanDat)

  # Progress report.
  print(paste("Performing permutation tests to calculate module preservation, iter", i))

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
    # nPerm = 1000, # nPerm will be determined by the function.
    null = "overlap",
    alternative = "greater",
    simplify = TRUE,
    verbose = FALSE
  )

  # Extract preservation statistics and p-values.
  preservation <- preservation[c("observed", "p.values")]

  # Get the maximum permutation test p-value.
  maxPval <- apply(preservation$p.values, 1, function(x) max(x, na.rm = TRUE))
  nModules <- length(unique(net$colors)) - 1
  params$nModules <- nModules

  # Modules removed if adjusted pvalue is greater than alpha.
  alpha <- 0.05
  modules_out <- names(maxPval)[maxPval > alpha / nModules]
  # names(net$MEs)[grepl(paste(modules_out,collapse="|"),names(net$MEs))]<-"MEgrey"

  # Drop nsig Modules (Set color to grey)
  net$colors[net$colors %in% modules_out] <- "grey"
  nModules <- length(unique(net$colors)) - 1

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
    softPower = power_beta,
    scale = TRUE,
    verbose = 0, indent = 0
  )
  net$MEs <- MEs$eigengenes

  # Fraction of un-assigned proteins.
  params$PercentGrayNodes <- 100 * table(net$colors)["grey"] / length(net$goodGenes)
  params$nSigModules <- nModules

  # Calculate modularity, q.
  membership <- as.numeric(as.factor(net$colors))
  q1 <- modularity(graph, membership, weights = edge_attr(graph, "weight"))
  params$q1 <- q1

  # Without "grey" nodes.
  v <- rownames(cleanDat)[!net$colors == "grey"]
  subg <- induced_subgraph(graph, v)
  membership <- as.numeric(as.factor(net$colors))
  membership <- membership[!net$colors == "grey"]
  q2 <- modularity(subg, membership, weights = edge_attr(subg, "weight"))
  params$q2 <- q2

  # Module coherence (aka PVE)
  modc <- as.matrix(
    propVarExplained(datExpr = t(cleanDat), colors = net$colors, MEs = net$MEs)
  )
  modc <- modc[!(rownames(modc) == "PVEgrey"), ]

  # Store in list.
  out[[i]] <- list(
    NetworkParams = params,
    ModuleColors = net$colors,
    ModulePreservation = preservation
  )

  # Save output every 10 iters.
  if (any((i == seq(1, nboot, 10)))) {
    collectGarbage()
    print(paste("Saving progress at:", Sys.time()))
    # Save results to file.
    saveRDS(out, file)
    # Push to GitHub every 100 iterations.
    if (git_push == TRUE & any((i == seq(1, nboot, 100)))) {
      print("Updating GitHub with progress.")
      gitstatus()
      gitadd()
      gitcommit(msg = paste("Exploring params iter", i, "complete."))
      gitpush()
    }
  }
} ## END LOOP.

#-------------------------------------------------------------------------------
# Load the results.
#-------------------------------------------------------------------------------

# Load Parameters and network stats.
files <- list.files(Rdatadir, pattern = "Stats")
files
file <- paste0(Rdatadir, "/", files[8])
file
out <- readRDS(file)
length(out)

# Load Dunnett tests.
# files <- list.files(Rdatadir)
# files
# file <- paste0(Rdatadir,"/",files[10])
# Dtests <- readRDS(file)
# length(Dtests)

# Network stats
result <- do.call(rbind, sapply(out, "[", 1))
result$iter <- c(1:nrow(result))
rownames(result) <- paste0("params_", 1:nrow(result))

# Module colors.
modcolors <- sapply(out, "[", 2)
names(modcolors) <- rownames(result)

# Preservation stats.
modstats <- sapply(out, "[", 3)
names(modstats) <- rownames(result)

# Rank stats.
result$score <- result$q2/(result$PercentGrayNodes/100)

# Examine Dtests results.
# result$nSigDunnettTests <- do.call(rbind,
#                                   lapply(Dtests,function(x) sum(x$nSigDunnettTest)))
# If KW pval is <0.05, then number of sig Dunnett tests.
# Dtests <- lapply(Dtests,function(x)
#  add_column(x,logic = x$nSigDunnettTest*as.numeric(x$p.adj<0.05),.after=5))
# Convergence, the number of modules with more than 2 sig changes.
# result$nModsConvergance <- do.call(rbind,
#                                   lapply(Dtests,function(x) sum(x$logic>1)))
# Write to file.
# write.csv(result,"params_result.csv")

#-------------------------------------------------------------------------------
#' ## Loop to recalculate Dunnett test results with p.adj, exploring convergance.
#-------------------------------------------------------------------------------

# CleanDat
file <- paste(Rdatadir, tissue, "TAMPOR_data_outliersRemoved.RDS", sep = "/")
cleanDat <- readRDS(file)
cleanDat <- log2(cleanDat)
dim(cleanDat)

# Sample info.
file <- paste(Rdatadir, tissue, "sample_info.RDS", sep = "/")
sample_info <- numericMeta <- readRDS(file)

# Load Sampled Network parameters and stats.
files <- list.files(Rdatadir)
file <- paste0(Rdatadir, "/", files[17])
out <- readRDS(file)

# Network stats
result <- do.call(rbind, sapply(out, "[", 1))
result$iter <- c(1:nrow(result))
rownames(result) <- paste0("params_", 1:nrow(result))

# Module colors.
modcolors <- sapply(out, "[", 2)
names(modcolors) <- rownames(result)

#### LOOP to perform Dunnett test.
output <- list()
n_iter <- 1

for (k in 94) {
  print(paste0("iteration ", k, "..."))
  # Calculate Module EigenProteins.
  MEs <- tmpMEs <- data.frame()
  MEList <- moduleEigengenes(t(cleanDat), colors = modcolors[[k]])
  MEs <- orderMEs(MEList$eigengenes)
  # let's be consistent in case prefix was added, remove it.
  colnames(MEs) <- gsub("ME", "", colnames(MEs))
  rownames(MEs) <- colnames(cleanDat)

  # Insure traits are in matching order.
  idx <- match(rownames(MEs), rownames(numericMeta))
  traits <- numericMeta[idx, ]
  all(rownames(traits) == rownames(MEs))

  # Define groups, the biological groups of interest.
  groups <- paste(traits$Tissue, traits$Sample.Model, sep = ".")
  groups[grepl("Cortex.WT", groups)] <- "WT.Cortex"
  groups[grepl("Striatum.WT", groups)] <- "WT.Striatum"
  # unique(groups)

  # Calculate Kruskal-Wallis pvalues for all modules (columns of MEs df).
  KWtest <- apply(MEs, 2, function(x) kruskal.test(x, as.factor(groups)))

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
  KWsigModules <- KW_results$Module[KW_results$p.adj < 0.05]
  # print(paste("nModules with p.adj < 0.05 =", length(KWsigModules)))
  # print(paste("nModules with FDR < 0.05 =", length(KW_results$Module[KW_results$FDR<0.05])))

  # Split Module EigenProtein (MEs) dm into a list of column vectors for lapply.
  ME_list <- split(as.matrix(MEs), rep(1:ncol(MEs), each = nrow(MEs)))
  names(ME_list) <- colnames(MEs)

  # Add vector of groups
  ME_list <- lapply(ME_list, function(x) data.frame(x = x, groups = groups))

  # Define levels for order of bars in plot.
  # levels <- c("WT.Cortex","ASD.Cortex","WT.Striatum","ASD.Striatum")
  levels <- c(
    "WT.Cortex", "WT.Striatum",
    "Cortex.KO.Shank2", "Striatum.KO.Shank2",
    "Cortex.KO.Shank3", "Striatum.KO.Shank3",
    "Cortex.HET.Syngap1", "Striatum.HET.Syngap1",
    "Cortex.KO.Ube3a", "Striatum.KO.Ube3a"
  )

  # Generate contrasts matrix for comparisons of interest.
  # contrasts <- makePairwiseContrasts(list("WT.Cortex","WT.Striatum"),list("ASD.Cortex","ASD.Striatum"))
  # contrasts
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
    plot <- ggplotVerboseBoxplot(x, g, levels, contrasts, color, stats = TRUE, method = "dunnett")
    plot_data[[i]] <- plot
    names(plot_data)[[i]] <- color
  }

  plots <- sapply(plot_data, "[", 1)

  # Extract Dunnett test stats.
  Dtest_stats <- sapply(plot_data, "[", 3)
  names(Dtest_stats) <- sapply(strsplit(names(Dtest_stats), "\\."), "[", 1)

  # Loop to add module column.
  for (i in 1:length(Dtest_stats)) {
    df <- Dtest_stats[[i]]
    namen <- names(Dtest_stats)[i]
    df <- add_column(df, Module = namen, .before = 1)
    Dtest_stats[[i]] <- df
  }

  # Just KW overall significant modules.
  res <- do.call(rbind, lapply(Dtest_stats, function(x) sum(x$P.adj < 0.05)))
  idx <- match(KW_results$Module, rownames(res))
  KW_results$nSigDunnettTest <- res[idx, ]

  output[[k]] <- KW_results

  # Save output.
  if (k == n_iter) {
    print(paste("Saving progress at:", Sys.time()))
    # Save results to file.
    file <- paste0(Rdatadir, "/", "Params_DunnettTests.RDS")
    saveRDS(output, file)
    gitstatus()
    gitadd()
    gitcommit(msg = paste("Dunnett tests for sampled params"))
    gitpush()
  }
}
