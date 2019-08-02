#!/usr/bin/env Rscript

## Usage:
# wgcna.r [data.Rds] [--parameters parameters.txt] 

# This script performs Weighted Gene (or Protein) Co-expression Analysis (WGCNA).
# The analysis is broken into three chunks:
# * Parsing the users input.
# * Prepaaring a named list of parameters for WGCNA.
# * Performing WGCNA and identifying appropriate soft-thresholding power, beta.
# * Analysis of the quality of the partition.

# The wgcna() function is a wrapper around WGCNA::BlockWiseModules(). This funciton
# performs WGCNA using the user's provided data and (optional) hyperparameters.
# If parameters are not provided, then the defaults will be used. 
# If a power for weighting the network is not provided,
# then it will be chosen such that the overall topology of the weighted 
# network is approximately scale free (R^2 > 0.8).

# Note: the default power is 9, and this is not appropriate for all datasets.
# Set power to NULL, if you wish to calculate the best soft-thresholding power.

###############################################################################
# # Load data for testing if no command line arguments passed.
# if (length(commandArgs(trailingOnly=TRUE)) == 0){
#   rm(list = ls())
#   if (.Device != "null device") { dev.off() }
#   cat("\f") # alternative is cat("\f")
#   options(stringsAsFactors = FALSE)
#   print("Using exprDat and saved params!")
#   dir <- "D:/projects/Synaptopathy-Proteomics/code/wgcna-hpo"
#   setwd(dir)
#   data_file  <- paste(dir, "exprDat.Rds", sep="/")
#   params_file <- paste(dir, "default_parameters.txt", sep="/")
#   exprDat <- readRDS(data_file)
#   temp_params <- read.delim(params_file, header = FALSE, col.names = c("Parameter","Value"))
#   temp_params$Value[temp_params$Value == "True"] <- TRUE
#   temp_params$Value[temp_params$Value == "False"] <- FALSE
#   temp_params <- temp_params[!(temp_params$Value == "None"),]
#   user_params <- as.list(temp_params$Value)
#   names(user_params) <- temp_params$Parameter
# }

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
data_file <- paste(dir, args$data, sep="/")
exprDat <- readRDS(data_file)

# If provided, parse the user's hyperparameters.
if (!is.na(args$parameters)) {
  params_file <- paste(dir, args$parameters, sep="/")
  temp_params <- read.delim(params_file, header = FALSE, col.names = c("Parameter","Value"))
  
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
  
  # Save defaults to file.
  if (!"default_parameters.txt" %in% list.files(dir)) {
    message("Saving default parameters!")
    write.table(as.matrix(unlist(default_params)),
                "default_parameters.txt", quote = FALSE, sep ="\t", col.names = FALSE)
  }
  
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
      message("No user defined parameters. Using default WGCNA parameters!")
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
    }, 
    finally = {
      # wrap up code if you want. 
    }
  )
} # Ends wgcna()

## Perform WGCNA!
results <- wgcna(exprDat, parameters)

#------------------------------------------------------------------------------
# ## Evaluate quality of the WGCNA partition.
#------------------------------------------------------------------------------

# FIXME: catch xlaunch error and tell the user to use bin/xlaunch

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
message(paste("... Median PVE    :", round(median(pve),3)))

#------------------------------------------------------------------------------
# Calculate quality as...
#------------------------------------------------------------------------------

# Sum of total within cluster variance explained.
# Sum of distance between eigengenes.
suppressPackageStartupMessages({
  require(igraph, quietly = TRUE)
})

#CH = (TWk/TBk) * (N-k)/(k-1)
# score = 1/CH
#kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc = "bicor")

## Caclulate total variance explained within a cluster.

# Calculate weighted signed adjacency matrix.
sink(tempfile())
adjm <- bicor(exprDat)^2
total_information <- sum(adjm)
sink(NULL)

# Complete graph.
g <- graph_from_adjacency_matrix(adjm, mode = "directed", weighted = TRUE, diag = FALSE)

# Get total variance explained for all clusters.
tve <- list()
for (color in unique(net$colors)){
  v <- names(net$colors[net$colors == color])
  subg <- induced_subgraph(g, vids = v)
  subadjm <- as.matrix(as_adjacency_matrix(subg, attr = "weight"))
  tv <- sum(subadjm)
  tve[color] <- tv * pve[match(color,substring(names(pve),4))]
}

tve <- sum(unlist(tve), na.rm = TRUE)

# total information is always the same. We want to maximize the total 
# variation explained by the network partition. This should maximize
# cluster quality (cohesiveness) while requiring that percent grey is 
# minimized. 

# Distance between module eigengenes.
sink(tempfile())
r <- bicor(net$MEs[!colnames(net$MEs) == "MEgrey"])
diss <- 1 - r
diag(diss) <- 0
sink(NULL)

# Maximize speration of clusters. 
N = dim(exprDat)[2]
k = nmodules
cluster_dispersion <- sum(diss)
score <- tve * cluster_dispersion * (N-k)/(k-1) * 1/100000000000
print(1/score)
