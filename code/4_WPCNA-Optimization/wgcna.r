#!/usr/bin/env Rscript

#' ---
#' title: wgcna.r
#' authors: Tyler W Bradshaw
#' description: Performs WGCNA given some data and parameters (optional).
#' usage: ./wgcna.R data.Rds [--parameters parameters.txt] 

#' This script performs Weighted Gene (or Protein) Co-expression Analysis (WGCNA).
#' The analysis is broken into three chunks:
#' * Parsing the users input.
#' * Prepaaring a named list of parameters for WGCNA.
#' * Performing WGCNA and identifying appropriate soft-thresholding power, beta.
#' * Analysis of the quality of the partition.

#' The wgcna() function is a wrapper around WGCNA::BlockWiseModules(). This funciton
#' performs WGCNA using the user's provided data and (optional) hyperparameters.
#' If parameters are not provided, then the defaults will be used. 
#' If a power for weighting the network is not provided,
#' then it will be chosen such that the overall topology of the weighted 
#' network is approximately scale free (R^2 > 0.8).

#' Note: the default power is 9, and this is not appropriate for all datasets.
#' Set power to NULL, if you wish to calculate the best soft-thresholding power.
#' ---

###############################################################################
## # Load data for testing if no command line arguments passed.
# if (length(commandArgs(trailingOnly=TRUE)) == 0){
#   rm(list = ls())
#   if (.Device != "null device") { dev.off() }
#   cat("\f") 
#   options(stringsAsFactors = FALSE)
#   message("Using wtDat and saved params!")
#   if (as.character(Sys.info()[1]) == "Linux"){
#	   message("Working on Linux OS!")
#	   dir <- getwd()
#	   setwd(dir)
#   }else if (as.character(Sys.info()[1]) == "Windows"){
#	   dir <- "D:/Projects/Synaptopathy-Proteomics/code/4_WPCNA-Optimization/EI"
#	   message("Working on Windows OS!")
#	   setwd(dir)
#   }
#   data_file  <- paste(dir, "wtDat.Rds", sep="/")
#   params_file <- paste(dir, "default_parameters.txt", sep="/")
#   exprDat <- readRDS(data_file)
#   temp_params <- read.delim(params_file, header = FALSE, col.names = c("Parameter","Value"))
#   temp_params$Value[temp_params$Value == "True"] <- TRUE
#   temp_params$Value[temp_params$Value == "False"] <- FALSE
#   temp_params <- temp_params[!(temp_params$Value == "None"),]
#   user_params <- as.list(temp_params$Value)
#   names(user_params) <- temp_params$Parameter
# }
#
###############################################################################

#------------------------------------------------------------------------------
# ## Parse the command line arguments.
#------------------------------------------------------------------------------

# Global options:
options(stringsAsFactors = FALSE)
require(argparser, quietly = TRUE)

# Check for input data before proceeding.
if(length(commandArgs(trailingOnly=TRUE)) == 0) {
  msg <- paste("Please provide normalized n x m matrix of expression data as input for WGCNA!",
               "Use ./wgcna.r --help for help.", sep = "\n")
  stop(msg)
}

# Parse the command line arguments.
p <- arg_parser("Perform WGCNA given a normalized n x m matrix of protein or gene expression data. ")
p <- add_argument(p, "data", 
                  help = paste("normalized n (sample) x m (gene) expression matrix",
                               "provided as a .Rds file in the same directory as this script"),
                  default = NULL)
p <- add_argument(p, "--parameters", short = "-p", 
                  help = "optional parameters for WGCNA algorithm", 
                  default = NULL)
args <- parse_args(p)

# Load data as n x m normalized expression matrix. 
dir <- getwd()
project_dir <- dirname(dirname(dir))
data_file <- paste(dir, args$data, sep="/")
exprDat <- readRDS(data_file)

# If provided, parse the user's hyperparameters.
if (!is.na(args$parameters)) {
  params_file <- paste(dir, args$parameters, sep="/")
  temp_params <- read.delim(params_file, header = FALSE, col.names = c("Parameter","Value"))
  unlink(params_file)

  # Replace True/False with TRUE/FALSE
  temp_params$Value[temp_params$Value == "True"] <- TRUE
  temp_params$Value[temp_params$Value == "False"] <- FALSE
  
  # Remove 'None' type arguments. These will be replaced with defaults.
  temp_params <- temp_params[!(temp_params$Value == "None"),]
  
  # Format as list. 
  user_params <- as.list(temp_params$Value)
  names(user_params) <- temp_params$Parameter
  
  # If no parameters passed, defaults will be generated by get_wgcna_params().
  } else if (is.na(args$parameters)) {
    user_params <- args$parameters     # NA
  } else {
    stop("Unable to parse user provided parameters.")
}

#------------------------------------------------------------------------------
# ## Define default WGCNA parameters.
#------------------------------------------------------------------------------

## Define a function that can suppress unwanted messages from a function. 
silently <- function(func, ...) {
  sink(tempfile())
  out <- func(...)
  sink(NULL)
  return(out)
}

## Given expression data matrix and optional user defined parameters, 
#  define parameters for WGCNA analysis. 

get_wgcna_params <- function(exprDat, overrides = NULL){
  # Global options:
  options(stringsAsFactors = FALSE)
  # WGCNA Parameters:
  keys <- c("weights", "checkMissingData", "blocks", "maxBlockSize", "blockSizePenaltyPower", 
            "nPreclusteringCenters", "randomSeed", 'loadTOM', 'corType', 'maxPOutliers', 
            "quickCor", "pearsonFallback", "cosineCorrelation", "power", "networkType", 
            "replaceMissingAdjacencies", "suppressTOMForZeroAdjacencies", 'TOMType', 
            'TOMDenom', 'getTOMs', 'saveTOMs', "saveTOMFileBase", "deepSplit", 
            "detectCutHeight", "minModuleSize", "maxCoreScatter", 'minGap', 'maxAbsCoreScatter',
            "minAbsGap", 'minSplitHeight', 'minAbsSplitHeight', 'useBranchEigennodeDissim', 
            'minBranchEigennodeDissim', 'stabilityLabels', 'stabilityCriterion', 
            'minStabilityDissim', 'pamStage', 'pamRespectsDendro', 'reassignThreshold', 
            'minCoreKME', 'minCoreKMESize', 'minKMEtoStay', 'mergeCutHeight', 'impute', 
            'trapErrors', 'numericLabels', 'nThreads', 'useInternalMatrixAlgebra', 
            'useCorOptionsThroughout', 'verbose', 'indent')
  default_params <- vector("list", length(keys))
  names(default_params) <- keys
  ## Define defaults.
  # Input data
  #default_params$exprDat <- exprDat
  #default_params$weights <- NULL
  # Data checking options
  default_params$checkMissingData <- FALSE
  # Options for splitting data into blocks
  #default_params$blocks                <- NULL
  default_params$maxBlockSize          <- 15000
  default_params$blockSizePenaltyPower <- 5
  default_params$nPreclusteringCenters <- as.integer(min(ncol(exprDat)/20, 100*ncol(exprDat)/default_params$maxBlockSize))
  default_params$randomSeed            <- 12345
  # load TOM from previously saved file?
  default_params$loadTOM <- FALSE
  # Network construction arguments: correlation options
  default_params$corType           <- "bicor"
  default_params$maxPOutliers      <- 1
  default_params$quickCor          <- 0
  default_params$pearsonFallback   <- "individual"
  default_params$cosineCorrelation <- FALSE
  # Adjacency function options
  default_params$power                         <- 13
  default_params$networkType                   <- "signed"
  default_params$replaceMissingAdjacencies     <- FALSE
  default_params$suppressTOMForZeroAdjacencies <- FALSE
  # Topological overlap options
  default_params$TOMType  <- "signed"  # c()
  default_params$TOMDenom <- "min"     # c("min","mean")
  # Saving or returning TOM
  #default_params$getTOMs         <- NULL
  default_params$saveTOMs        <- FALSE 
  default_params$saveTOMFileBase <- "blockwiseTOM"
  # Basic tree cut options
  default_params$deepSplit       <- 2
  default_params$detectCutHeight <- 0.995 
  default_params$minModuleSize   <- min(20, ncol(exprDat)/2)
  # Advanced tree cut options
  # default_params$maxCoreScatter           <- NULL 
  # default_params$minGap                   <- NULL
  # default_params$maxAbsCoreScatter        <- NULL
  # default_params$minAbsGap                <- NULL
  # default_params$minSplitHeight           <- NULL
  # default_params$minAbsSplitHeight        <- NULL
  default_params$useBranchEigennodeDissim <- FALSE
  default_params$minBranchEigennodeDissim <- 0.15
  # default_params$stabilityLabels          <- NULL
  default_params$stabilityCriterion       <- "Individual fraction" # c("Individual fraction", "Common fraction")
  # default_params$minStabilityDissim       <- NULL
  default_params$pamStage                 <- TRUE
  default_params$pamRespectsDendro        <- TRUE
  # Gene reassignment, module trimming, and module "significance" criteria
  default_params$reassignThreshold <- 1e-6
  default_params$minCoreKME        <- 0.5
  default_params$minCoreKMESize    <- round(min(20, ncol(exprDat)/2)/3)
  default_params$minKMEtoStay      <- 0.3
  # Module merging options
  default_params$mergeCutHeight <- 0.15
  default_params$impute         <- TRUE
  default_params$trapErrors     <- FALSE
  # Output options
  default_params$numericLabels <- FALSE
  # Options controlling behaviour
  default_params$nThreads                 <- 8
  default_params$useInternalMatrixAlgebra <- FALSE
  default_params$useCorOptionsThroughout  <- TRUE
  default_params$verbose                  <- 0 
  default_params$indent                   <- 0
  # # Save defaults to file.
  # if (!"default_parameters.txt" %in% list.files(dir)) {
  #   message("Saving default parameters!")
  #   write.table(as.matrix(unlist(default_params)),
  #               "default_parameters.txt", quote = FALSE, sep ="\t", col.names = FALSE)
  # }
  # If provided, overwrite default parameters with user provided ones.
  if (inherits(overrides, "list")) {
    idx <- match(names(overrides), names(default_params))
    parameters <- default_params
    parameters[idx] <- overrides
    # Insure that params are the correct data type. Ignore NULLS.
    parameters <- lapply(parameters, function(x) if(!is.null(x)) { type.convert(x) })
    # Make sure factors are converted back to characters!
    idx <- c(1:length(parameters))[unlist(lapply(parameters, function(x) is.factor(x)))]
    parameters[idx] <- unlist(lapply(parameters[idx], function(x) as.character(x)))
    return(parameters)
    # Otherwise, use the defaults.
    } else if (is.null(user_params) | is.na(user_params)) {
      message("No user defined parameters. Using default parameters!")
      parameters <- default_params
      return(parameters)
      } else {
        stop("please provide a list of parameters, or use the defaults.")
      }
}

# Get parameters for WGCNA and save them to file.
parameters <- get_wgcna_params(exprDat, overrides = user_params)

#------------------------------------------------------------------------------
# ## Define a function to perform WGCNA.
#------------------------------------------------------------------------------

wgcna <- function(exprDat, parameters) {
  # Use tryCatch to handle errors caused by parameters which return a single module.
  tryCatch(
    {
      ## Global options and imports. 
      options(stringsAsFactors=FALSE)
      suppressPackageStartupMessages({
        require(WGCNA)
        require(doParallel)
        require(parallel)
      })
      # Allow parallel WGCNA calculations if nThreads is > 0.
      if (parameters$nThreads > 0) {
        # Use sink to supress unwanted output.
        temp <- tempfile()
        sink(temp) 
        allowWGCNAThreads(nThreads = parameters$nThreads)
        sink(NULL)
        unlink(temp)
        clusterLocal <- makeCluster(c(rep("localhost", parameters$nThreads)), type = "SOCK")
        registerDoParallel(clusterLocal)
      }
      # Function to get best power, suppressing unwanted output with sink.
      pickPower <- function(...){
        temp <- tempfile()
        sink(temp) 
        output <- WGCNA::pickSoftThreshold(...)
        sink(NULL)
        unlink(temp)
        return(output)
      }
      # If no power was provided, then determine the best power to achieve ~scale free toplogy. 
      if (is.null(parameters$power)) {
        message("Picking the best soft thresholding power to achieve scale free fit > 0.8!")
        sft <- pickPower(
          data        = exprDat,
          dataIsExpr  = TRUE,
          weights     = parameters$weights,
          RsquaredCut = 0.80,
          powerVector = seq(2, 20, by = 1.0),
          removeFirst = FALSE, 
          nBreaks     = 10, 
          corFnc      = parameters$corType,
          corOptions  = list(use = 'p'),
          blockSize   = parameters$maxBlockSize,
          networkType = parameters$networkType,
          moreNetworkConcepts = FALSE,
          gcInterval  = NULL,
          verbose     = parameters$verbose,
          indent      = parameters$indent
        )
        # Calculate the best power_beta.
        parameters$power <- sft$powerEstimate
        message(paste("Best soft-threshold power, beta:", sft$powerEstimate))
      }
      ## Perform WGCNA by calling the blockwiseModules() function.
      # Function to supress unwanted output from the blockwiseModules with sink().
      blockwiseWGCNA <- function(...){
        temp <- tempfile()
        sink(temp) 
        output <- WGCNA::blockwiseModules(...)
        sink(NULL)
        unlink(temp)
        return(output)
      }
      net <- blockwiseWGCNA(
        # Input data
        datExpr               = exprDat, 
        weights               = parameters$weights,
        # Data checking options
        checkMissingData      = parameters$checkMissingData,
        # Options for splitting data into blocks
        blocks                = parameters$blocks,
        maxBlockSize          = parameters$maxBlockSize,
        blockSizePenaltyPower = parameters$blockSizePenaltyPower,
        nPreclusteringCenters = parameters$nPreclusteringCenters,
        randomSeed            = parameters$randomSeed,
        # load TOM from previously saved file? This will speed things up. 
        loadTOM               = parameters$loadTOM,
        # Network construction arguments: correlation options
        corType               = parameters$corType,
        maxPOutliers          = parameters$maxPOutliers, 
        quickCor              = parameters$quickCor,
        pearsonFallback       = parameters$pearsonFallback,
        cosineCorrelation     = parameters$cosineCorrelation,
        # Adjacency function options
        power                     = parameters$power,
        networkType               = parameters$networkType,
        replaceMissingAdjacencies = parameters$replaceMissingAdjacencies,
        # Topological overlap options
        TOMType                       = parameters$TOMType,
        TOMDenom                      = parameters$TOMDenom,
        suppressTOMForZeroAdjacencies = parameters$suppressTOMForZeroAdjacencies,
        suppressNegativeTOM           = parameters$suppressNegativeTOM,
        # Saving or returning TOM
        getTOMs         = parameters$getTOMs,
        saveTOMs        = parameters$saveTOMs, 
        saveTOMFileBase = parameters$saveTOMFileBase,
        # Basic tree cut options
        deepSplit       = parameters$deepSplit,
        detectCutHeight = parameters$detectCutHeight,
        minModuleSize   = parameters$minModuleSize,
        # Advanced tree cut options
        maxCoreScatter           = parameters$maxCoreScatter, 
        minGap                   = parameters$minGap,
        maxAbsCoreScatter        = parameters$maxAbsCoreScatter, 
        minAbsGap                = parameters$minAbsGap,
        minSplitHeight           = parameters$minSplitHeight, 
        minAbsSplitHeight        = parameters$minAbsSplitHeight,
        useBranchEigennodeDissim = parameters$useBranchEigennodeDissim,
        minBranchEigennodeDissim = parameters$mergeCutHeight,
        stabilityLabels          = parameters$stabilityLabels,
        stabilityCriterion       = parameters$stabilityCriterion,
        minStabilityDissim       = parameters$minStabilityDissim,
        pamStage                 = parameters$pamStage, 
        pamRespectsDendro        = parameters$pamRespectsDendro,
        # Gene reassignment, module trimming, and module "significance" criteria
        reassignThreshold  = parameters$reassignThreshold,
        minCoreKME         = parameters$minCoreKME, 
        minCoreKMESize     = round(parameters$minModuleSize/3),
        minKMEtoStay       = parameters$minKMEtoStay,
        # Module merging options
        mergeCutHeight = parameters$mergeCutHeight, 
        impute         = parameters$impute, 
        trapErrors     = parameters$trapErrors, 
        # Output options
        numericLabels = parameters$numericLabels,
        # Options controlling behaviour
        nThreads                 = parameters$nThreads,
        useInternalMatrixAlgebra = parameters$useInternalMatrixAlgebra,
        useCorOptionsThroughout  = parameters$useCorOptionsThroughout,
        verbose                  = parameters$verbose, 
        indent                   = parameters$indent)
      return(list("data" = exprDat, "network" = net, "hyperparameters" = parameters))
    }, 
    # Catch error messages. 
    error = function(cond) {
      stop("Error: unable to perform WGCNA!")
    },
    # Catch warning messages.
    # TryCatch() will abort code prematurely, so we need to close sink opened above.
    # Otherwise stdout will be sopped up by sink.
    warning = function(cond) {
      sink(NULL) 
      unlink(temp)
      message("Warning: unable to complete analysis! Likely cause: WGCNA returned 0 or 1 module.")
      print(100) # a bad score
      quit()
    }
  )
} 
# Ends wgcna()

## Perform WGCNA!
results <- wgcna(exprDat, parameters)

#------------------------------------------------------------------------------
# ## Exctract key WGCNA results. 
#------------------------------------------------------------------------------

suppressPackageStartupMessages({
  require(data.table, quietly = TRUE)
})

# Extract WGCNA results.
exprDat <- results$data
net <- results$network
params <- results$hyperparameters

# The number of modules and percent grey.
nmodules <- length(unique(results$network$colors)) - 1 # exclude grey

# Stop if nModules < 2, as we cannot compute quality statistics with fewer than 2 modules. 
if (nmodules < 2) {
  message("Error: cannot compute quality indices if nModules < 2.0")
  print(100) # A bad score
  quit()
}

# Progress report.
msg <- paste("... Total nModules:", nmodules)
message(msg)

# Calculate percent grey (unclustered nodes):
is_grey <- net$colors == "grey"
percent_grey <- sum(is_grey)/length(is_grey)
message(paste("... Percent grey  :", round(percent_grey,3)))

# Calculate median percent variance explained, module cohesivness:
pve <- WGCNA::propVarExplained(exprDat, net$colors, net$MEs, corFnc = params$corType)
pve <- pve[!names(pve) == "PVEgrey"]
message(paste("... Cluster medPVE:", round(median(pve),3)))

# Calculate percent total variation explained by the clustering.
adjm <- silently(bicor, exprDat)
subadjm <- silently(bicor, exprDat[,!is_grey])
ptve <- sum(subadjm^2) / sum(adjm^2)
loss_ptve <- 1 - ptve
message(paste("... Partition PVE :", round(ptve,3)))

# # Write WGCNA partition to file.
# v <- as.numeric(as.factor(net$colors))
# file <- paste(dir,"wgcna-partition.csv")
# fwrite(as.data.table(t(v)),"wgcna-partition.txt", append = TRUE) 

#------------------------------------------------------------------------------
# ## Calculate Modularity, the quality of the partition.
#------------------------------------------------------------------------------
# Modularity typically only applies to unsigned graphs.
# Modularity for signed graphs: Gomez et al., 2018 
# REF: (https://arxiv.org/abs/0812.3030)

suppressPackageStartupMessages({
  require(reshape2, quietly = TRUE)
  require(igraph, quietly = TRUE)
})

# Function to check if value is even.
# If params$power is even then we will need to enforce sign of adjm.
is_even <- function(x){
  return((x %% 2) == 0)
}

# Calculate signed, weighted adjacency matrix.
# Weight (power) affects sign and modularity. If power is even, then ensure that
# sign of interaction is enforced. (negative value ^ even power = positive)
if (is_even(params$power)) {
	r <- adjm
	r[r<0] <- -1
	r[r>0] <- 1
        signed_adjm <- r * (adjm^params$power)
} else {
	signed_adjm <- adjm^params$power
}

# Write cluster info to file.
script_dir <- paste(project_dir,"bin","radalib", sep="/")
cluster_file <- paste(script_dir,"clusters.clu", sep="/")
cl <- as.matrix(as.numeric(as.factor(net$colors[match(colnames(adjm), names(net$colors))])))
colnames(cl) <- paste("*Vertices",ncol(signed_adjm))
write.table(cl, quote = FALSE, file = cluster_file, row.names = FALSE, col.names = TRUE)

# Write network to file in Pajak format: use fwrite for faster performance!
network_file <- paste(script_dir, "network.net", sep="/")
n <- signed_adjm
colnames(n) <- rownames(n) <- c(1:ncol(n))
edge_list <- as.data.table(na.omit(melt(n)))
colnames(edge_list) <- c("protA","protB","weight")
v <- as.data.table(paste(seq(1,ncol(n)), " \"", seq(1,ncol(n)), "\"", sep = ""))
write.table(paste("*Vertices", dim(n)[1]), file = network_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
fwrite(v, file = network_file, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE, append = TRUE)
write.table("*Edges", file = network_file, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
fwrite(edge_list, file = network_file, sep = " ", col.names = FALSE, append = TRUE)

# Build a command to be passed to radalib.
script <- "./modularity_calculation.exe"
type <- "WS" # Weighted signed network. 
cmd <- paste(script, "network.net", "clusters.clu", type)

# Evaluate modularity of the partition.
# Need to be in radalib/tools directory!
setwd(script_dir)
result <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
setwd(dir)

# Parse the result.
x <- trimws(result[grep(" Q = ", result)])
message(paste("...", x))
Q <- as.numeric(unlist(strsplit(x,"\ "))[4])

# Return score as the inverse of modularity. Bigger is better = smaller.
score <- 1/Q
print(score)

# Remove .net and .clu files.
unlink(network_file)
unlink(cluster_file)

quit()

#------------------------------------------------------------------------------
# ## Calculate modified Davies-Bouldin quality index
#------------------------------------------------------------------------------

# DB Index - quality of network partition:

# Rij = si + sj / dij    # si and sj should be small, dij should be big!

# For cluster i and its most similar cluster j compute
# the average distance between its nodes and its center
# as si. For cluster j do the same. Normalize the sum
# of si and sj by the distance between the center of
# clusters i and j.

# DB = 1/k * sumi,k(maxRij) 

# For all modules find max Rij (closest cluster) and sum this quantity.

# WGCNA interpretation.
# Most similar cluster = Closest ME.
# Similarity between grey for all clusters should be minimized! 
# si and sj are average of distance to module center (1-kME).

suppressPackageStartupMessages({
  require(reshape2, quietly = TRUE)
  require(dplyr, quietly = TRUE)
})

# Calculate percent total variation explained by the clustering.
adjm <- silently(bicor, exprDat)
subadjm <- silently(bicor, exprDat[,!is_grey])
ptve <- sum(subadjm^2) / sum(adjm^2)
loss_ptve <- 1 - ptve
message(paste("... Partition PVE :", round(ptve,3)))

# KME is module membership. Bicor between protein and modules ME (its center).
kME <- signedKME(exprDat, datME=net$MEs, outputColumnName = "kME", corFnc = params$corType)

# For all clusters find average distance to center of the cluster (si = 1 - kME).
# small kME then larger distance to center of cluster.
si <- list()
color_vec <- unique(net$colors)
for (color in color_vec) {
  v <- names(net$colors[net$colors == color])
  subkME <- subset(kME, rownames(kME) %in% v)
  idy <- match(color, substring(colnames(subkME),4))
  si[color] <- mean(1-subkME[,idy])
}

# Calculate ME distance matrix in order to find closest ME for each cluster.
MEs <- net$MEs
adjm <- silently(bicor, MEs)
diag(adjm) <- NA
datME <- na.omit(reshape2::melt(adjm))
colnames(datME) <- c("MEi","MEj","bicor")

# Calculate dij, the distance between MEs (dij = 1-bicor).
datME$dij <- 1 - datME$bicor

# Get average distance to center of cluster for MEi and MEj
datME$si <- unlist(lapply(as.list(datME$MEi), function(x) si[substring(x,3)]))
datME$sj <- unlist(lapply(as.list(datME$MEj), function(x) si[substring(x,3)]))

# Calculate similarity metric Rij.
datME$Rij <- (datME$si + datME$sj) / datME$dij

# Modify DB index by also calculating distance between all modules and grey.
# For all clusters, similarity to grey should be minimized:
datMEgrey <- subset(datME, datME$MEi=="MEgrey")

# For remaining clusters (excluding grey), similarity to nearest neighbor 
# should be minimized (better clusering = greater cluster seperation).
# Get top (max Rij) for each MEi.
datMEij <- subset(datME, !(datME$MEi == "MEgrey" | datME$MEj=="MEgrey"))
datMEij <- as.data.frame(datMEij %>% dplyr::group_by(MEi) %>% dplyr::top_n(1, Rij))

# Calculate modified DB index as Rij + Rgrey (similarity between all clusters and grey).
Rgrey <- sum(datMEgrey$Rij) # Similarity to grey should bi minimized!
Rij <- sum(datMEij$Rij)
k <- nmodules
modDB <- (Rij + Rgrey)/k
message(paste("... Modified DB   :", round(modDB,3))) # smaller is better!

# Calculate quality score as modDB * loss_ptve (the percent unexplained variance).
score <- modDB/ptve # smaller is better!
print(score)

quit()

###############################################################################
## ENDOFILE ##
###############################################################################

###############################################################################
## Additional quality metrics:

## NOTES:
# W = within cluster variance; if W is small, then PVEk ~ 100. As k increases, so does PVEk. 
# T = total cluster variance --> does not change based on partition Pi.

# Measure heterogenity of partition:
# Diameter of partition is the maximum diameter of its clusters. This can be calculated with any dissimilarity matrix.

# The total number of partitions for n items is given as:
# K^n/K!

#------------------------------------------------------------------------------
# ## Evaluate quality of partition... Modified Calinski-Harabasz Index.
#------------------------------------------------------------------------------

# Calinksi-Harabasz definition of partition quality:

#  Variance Ratio * (N-k)/(k-1)

# The variance Ratio is the ratio between within cluster variance and between 
# cluster variance: 

# CH = TBk/TWk * (N-k)/(k-1) 

# When:
# TWk is the sum of squared distances between the nodes in a cluster and its center.
# TBk is the sum squared distances between the center of all clusters and the 
# center of the network. 

# Here, Twk is the sum of the distances between proteins and their ME. This is 
# equivalent to: sum((1-kme)^2)
# As kme is a proteins module membership (the correlation between its expression 
# vector and its ME). 

# Between cluster variance, TBk is the sum of squared distances between module 
# centers and either the:
#     * center of the network
#     * the center of all clusters
# Define the center of the network as the average of all MEs.
# Define distance between clusters as distance between MEs.

# A higher CH score indicates a better partition.

# KME is module membership. Bicor between protein and modules ME (its center).
kME <- signedKME(exprDat, datME=net$MEs, corFnc = params$corType)

# For each cluster compute sum of squared distance to ME.
wk <- list()
color_vec <- unique(net$colors)
for (color in color_vec) {
  v <- names(net$colors[net$colors == color])
  subkME <- subset(kME, rownames(kME) %in% v)
  idy <- match(color, substring(colnames(subkME),4))
  wk[color] <- sum(1-subkME[,idy])^2
}

# Total within cluster variance.
wk <- unlist(wk)
idx <- names(wk) == "grey"
TWk <- as.vector(sum(wk[!idx]) - wk[idx]) # Minimize this quantity! More negative is better.

## Calculate total between cluster variance:
method = c("centroid", "all")[2]
if (method == "centroid") {
  # Total Between Cluster variance is sum of squared distances between modules 
  # centers (their MEs) and the center of the network, the average ME. 
  MEs <- net$MEs
  MEcentroid <- apply(MEs,1,mean)
  adjm <- silently(bicor,(cbind(MEs,MEcentroid)))
  Bk <- (1 - adjm[,"MEcentroid"])^2 # Distance beween all clusters and network centroid.
  TBk <- sum(Bk) # Should be maximized!
  
} else if (method == "all") {
  # Total Between Cluster variance is sum of squared distances between modules MEs.
  MEs <- net$MEs
  adjm <- silently(bicor, MEs)
  Bk <- (1 - adjm)^2 # Distance beween all clusters centers.
  TBk <- sum(Bk) # Should be maximized!
}

# Return score. 
N <- dim(exprDat)[2] # All nodes!!
k <- nmodules
score <- (TBk/TWk)^-1 * (N-k)/(k-1) # Smaller is better!
print(score)

quit()


#------------------------------------------------------------------------------
## Calculate simple quality stats.
#------------------------------------------------------------------------------

# Variance ratio and PVE will not penalize small number of clustered nodes.

# Given a dendrogram, the total variance is sum of all heights (distances).
# Within cluster variance is sum of all intra-cluster heights distances.
# The Between cluster variance can then be computed from:
# T = W + B; 
# B = T - W

# Calculate dissimilarity matrix. 
r <- silently(bicor,exprDat)
adjm <- ((1 + r) / 2)^params$power
diss <- 1 - TOMsimilarity(adjm,
                          TOMType  = params$TOMType,
                          TOMDenom = params$TOMDenom,
                          suppressTOMForZeroAdjacencies = params$suppressTOMForZeroAdjacencies,
                          useInternalMatrixAlgebra      = params$useInternalMatrixAlgebra,
                          verbose = params$verbose,
                          indent  = params$indent)
colnames(diss) <- colnames(adjm)
rownames(diss) <- rownames(adjm)

# Calculate distance matrix, D and T, W, and B:
D <- as.matrix(as.dist(diss))
Total <- sum(D)
out <- rownames(D) %in% names(net$colors[is_grey])
W <- sum(D[!out, !out])
B <- Total - W

# Variance Ratio:
R <- W/B
message(paste("... Variance Ratio:", R)) # Smaller is better!

# Percent variance explained:
PVE <- 100*(1 - W/Total) 
message(paste("... Total VarEx(%):", PVE))

#------------------------------------------------------------------------------
# ## Calculate quality as... Total variance explained by clustering.
#------------------------------------------------------------------------------

# Sum of total within cluster variance explained.
# Does not account for distances between clusters!

suppressPackageStartupMessages({
  require(igraph, quietly = TRUE)
})

# Calculate total amount of information in the network as the sum of its adjm.
# Raise to the power of two so that all adjacencies are positive.
adjm <- silently(bicor,exprDat)
total_information <- sum(adjm^2)

# Generate igraph of this adjacency matrix.
g <- graph_from_adjacency_matrix(adjm, mode = "directed", weighted = TRUE, diag = FALSE)

# Calculate total variance explained within a cluster.
# The variance in a cluster is the sum of the squared adjacency matrix.
# The total variance explained is the percent_variance_explained * total_variance
ve <- list()
for (color in unique(net$colors)){
  v <- names(net$colors[net$colors == color])
  subg <- induced_subgraph(g, vids = v)
  subadjm <- as.matrix(as_adjacency_matrix(subg, attr = "weight"))
  total_var <- sum(subadjm^2) # Square such that all positive.
  ve[color] <- total_var * pve[match(color,substring(names(pve),4))]
}

# The sum of all variance explained by the clustering:
tve <- sum(unlist(ve), na.rm = TRUE)     # Grey is NA.
cluster_quality <- tve # bigger is better!

# total information is always the same. We want to maximize the total
# variation explained by the network partition. This should maximize
# cluster quality (cohesiveness) while requiring that percent grey is
# minimized.
quality <- 1/k * ((total_information - tve)/total_information) # Percent variance explained normalized by k.
score = 1/quality # Minimize!
print(score)

quit()

# In addition to module cohesiveness, modules should be well sperated.
# Calculate distance between module eigengenes.
r <- silently(bicor, net$MEs)
diss <- 2 - r
diag(diss) <- 0
cluster_dispersion <- sum(diss)  # Bigger is better!

# Calculate quality score:
N = length(net$colors[!is_grey])
k = nmodules
score <- cluster_quality * cluster_dispersion * (N-k)/(k-1) # Bigger is better!
print(10000000000/score)

quit()



