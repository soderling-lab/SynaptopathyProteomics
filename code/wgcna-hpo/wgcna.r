#!/usr/bin/env Rscript

## wgcna.R [exprDat.Rds] [parameters.txt]

#------------------------------------------------------------------------------
# ## Parse the command line input.
#------------------------------------------------------------------------------

# If provided, parse the user provided hyperparemeters.
args <- commandArgs(trailingOnly=TRUE)
nargs <- length(args)

if (nargs == 0) { stop("Please provide input expression data!")
  } else if (nargs > 0) {
    
    # Load the expression data and parameters. 
    dir <- getwd()
    data_file  <- paste(dir, args[1], sep="/")
    params_file <- paste(dir, args[2], sep="/")
    
    # Load data as n x m normalized expression matrix. 
    exprDat <- log2(t(readRDS(data_file)))
    temp_params <- read.delim(params_file, header = FALSE, col.names = c("Parameter","Value"))
    
    # Replace True/False with TRUE/FALSE
    temp_params$Value[temp_params$Value == "True"] <- TRUE
    temp_params$Value[temp_params$Value == "False"] <- FALSE
    
    # Remove 'None' type arguments. These will be replaced with defaults.
    temp_params <- temp_params[!(temp_params$Value == "None"),]
    
    # Format as list. 
    user_params <- as.list(temp_params$Value)
    names(user_params) <- temp_params$Parameter
  }

##### LOAD INPUT IF JUST TESTING ##############################################
dir <- "D:/projects/Synaptopathy-Proteomics/code/wgcna-hpo"
setwd(dir)
data_file  <- paste(dir, "exprDat.Rds", sep="/")
params_file <- paste(dir, "parameters.txt", sep="/")
exprDat <- log2(t(readRDS(data_file)))
temp_params <- read.delim(params_file, header = FALSE, col.names = c("Parameter","Value"))

# Replace True/False with TRUE/FALSE.
temp_params$Value[temp_params$Value == "True"] <- TRUE
temp_params$Value[temp_params$Value == "False"] <- FALSE

# Remove 'None' type arguments. These will be replaced with defaults.
temp_params <- temp_params[!(temp_params$Value == "None"),]

# Format as list.
user_params <- as.list(temp_params$Value)
names(user_params) <- temp_params$Parameter

###############################################################################

#------------------------------------------------------------------------------
# ## Define a function to perform WGCNA.
#------------------------------------------------------------------------------

# This is a wrapper function around WGCNA::BlockWiseModules() that performs
# WGCNA using the user's provided data and hyperparameters.If parameters are 
# not provided, then the defaults will be used. If a power_beta is not provided,
# then it will be chosen such that the overall topology of the weighted 
# network is approximately scale free.

wgcna <- function(exprDat, parameters=NULL){

  ## Global options and imports. 
  options(stringsAsFactors = FALSE)
  suppressPackageStartupMessages({
    require(WGCNA)
    require(doParallel)
    require(parallel)
  })
  
  ## Load WGCNA defaults stored in Bin/.
  # These will be used if the user does not provide any.
  source("./defaults.r")
  
  ## If provided, parse the user provided parameters.
  if (!exists("parameters")) {
    params <- default_params
  } else if (inherits(parameters,"list")) {
    user_params <- parameters
  } else {
    stop("please provide a list of parameters, or use the defaults.")
  }
  
  # Overwrite default parameters with user provided ones.  
  idx <- match(names(user_params), names(default_params))
  params <- default_params
  params[idx] <- user_params
  
  # Insure that params are the correct data type. Ignore NULLS.
  params <- lapply(params, function(x) if(!is.null(x)) { type.convert(x) })
  # Make sure factors are converted back to characters!
  idx <- c(1:length(params))[unlist(lapply(params, function(x) is.factor(x)))]
  params[idx] <- unlist(lapply(params[idx], function(x) as.character(x)))
  
  # Allow parallel WGCNA calculations if nThreads is > 0.
  if (params$nThreads > 0) {
    # Use sink to supress unwanted output.
    temp <- tempfile()
    sink(temp) 
    allowWGCNAThreads(nThreads = params$nThreads)
    sink(NULL)
    unlink(temp)
    clusterLocal <- makeCluster(c(rep("localhost", params$nThreads)), type = "SOCK")
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
  if (is.null(params$power)) {
    sft <- pickPower(
      data        = exprDat,
      dataIsExpr  = TRUE,
      weights     = params$weights,
      RsquaredCut = 0.80,
      powerVector = seq(2, 20, by = 1.0),
      removeFirst = FALSE, 
      nBreaks     = 10, 
      corFnc      = params$corType,
      corOptions  = list(use = 'p'),
      blockSize   = params$maxBlockSize,
      networkType = params$networkType,
      moreNetworkConcepts = FALSE,
      gcInterval  = NULL,
      verbose     = params$verbose,
      indent      = params$indent
    )
  }
  
    # Calculate the best power_beta.
    params$power <- sft$powerEstimate
    
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
    
    net <-blockwiseWGCNA(
      # Input data
      datExpr               = exprDat, 
      weights               = params$weights,
      
      # Data checking options
      checkMissingData      = params$checkMissingData,
      
      # Options for splitting data into blocks
      blocks                = params$blocks,
      maxBlockSize          = params$maxBlockSize,
      blockSizePenaltyPower = params$blockSizePenaltyPower,
      nPreclusteringCenters = params$nPreclusteringCenters,
      randomSeed            = params$randomSeed,
      
      # load TOM from previously saved file? This will speed things up. 
      loadTOM               = params$loadTOM,
      
      # Network construction arguments: correlation options
      corType               = params$corType,
      maxPOutliers          = params$maxPOutliers, 
      quickCor              = params$quickCor,
      pearsonFallback       = params$pearsonFallback,
      cosineCorrelation     = params$cosineCorrelation,
      
      # Adjacency function options
      power                     = params$power,
      networkType               = params$networkType,
      replaceMissingAdjacencies = params$replaceMissingAdjacencies,
      
      # Topological overlap options
      TOMType                       = params$TOMType,
      TOMDenom                      = params$TOMDenom,
      suppressTOMForZeroAdjacencies = params$suppressTOMForZeroAdjacencies,
      suppressNegativeTOM           = params$suppressNegativeTOM,
      
      # Saving or returning TOM
      getTOMs         = params$getTOMs,
      saveTOMs        = params$saveTOMs, 
      saveTOMFileBase = params$saveTOMFileBase,
      
      # Basic tree cut options
      deepSplit       = params$deepSplit,
      detectCutHeight = params$detectCutHeight,
      minModuleSize   = params$minModuleSize,
      
      # Advanced tree cut options
      maxCoreScatter           = params$maxCoreScatter, 
      minGap                   = params$minGap,
      maxAbsCoreScatter        = params$maxAbsCoreScatter, 
      minAbsGap                = params$minAbsGap,
      minSplitHeight           = params$minSplitHeight, 
      minAbsSplitHeight        = params$minAbsSplitHeight,
      useBranchEigennodeDissim = params$useBranchEigennodeDissim,
      minBranchEigennodeDissim = params$mergeCutHeight,
      stabilityLabels          = params$stabilityLabels,
      stabilityCriterion       = params$stabilityCriterion,
      minStabilityDissim       = params$minStabilityDissim,
      pamStage                 = params$pamStage, 
      pamRespectsDendro        = params$pamRespectsDendro,
      
      # Gene reassignment, module trimming, and module "significance" criteria
      reassignThreshold  = params$reassignThreshold,
      minCoreKME         = params$minCoreKME, 
      minCoreKMESize     = round(params$minModuleSize/3),
      minKMEtoStay       = params$minKMEtoStay,
      
      # Module merging options
      mergeCutHeight = params$mergeCutHeight, 
      impute         = params$impute, 
      trapErrors     = params$trapErrors, 
      
      # Output options
      numericLabels = params$numericLabels,
      
      # Options controlling behaviour
      nThreads                 = params$nThreads,
      useInternalMatrixAlgebra = params$useInternalMatrixAlgebra,
      useCorOptionsThroughout  = params$useCorOptionsThroughout,
      verbose                  = params$verbose, 
      indent                   = params$indent)
  
    # Output:
    return(list("network" = net, "hyperparameters" = params))
    } 
# END FUNCTION.

#------------------------------------------------------------------------------
# ## Perform WGCNA!
#------------------------------------------------------------------------------

results <- wgcna(exprDat, parameters = user_params)

#------------------------------------------------------------------------------
# ## Evaluate the quality of the WGCNA partition.
#------------------------------------------------------------------------------

suppressPackageStartupMessages({
  require(clusterSim)
})

# Extract WGCNA results.
net <- results$network
params <- results$parameters

# The number of modules.
nmodules <- length(unique(net$colors)) - 1 # exclude grey

# Calculate percent grey (unclustered nodes):
is_grey <- net$colors == "grey"
percent_grey <- sum(is_grey)/length(is_grey)

# Calculate median percent variance explained, module cohesivness:
pve <- WGCNA::propVarExplained(exprDat, net$colors, net$MEs, corFnc = params$corType)
median_pve <- median(pve[!names(pve) == "PVEgrey"])

# Remove grey nodes from data.
subDat <- exprDat[,!is_grey]
  
# Calculate weighted signed adjacency matrix.
sink(tempfile())
adjm <- ((1 + bicor(subDat))/ 2)^params$power
sink(NULL)

# Calculate TOM dissimilarity. 
diss <- 1 - TOMsimilarity(
  adjm, 
  TOMType  = params$TOMType,
  TOMDenom = params$TOMDenom,
  suppressTOMForZeroAdjacencies = params$suppressTOMForZeroAdjacencies,
  useInternalMatrixAlgebra      = params$useInternalMatrixAlgebra,
  verbose  = params$verbose,
  indent   = params$indent
)

# Calcualte network cluster quality indices.
cl <- as.integer(as.numeric(as.factor(result$network$colors[!is_grey])))
db <- index.DB(adjm, cl)$DB
ch <- index.G1(adjm, cl)
sc <- index.S(as.dist(diss), cl)

# Return results.
quality_results <- list(
  n_Modules         = nmodules,
  Percent_grey      = percent_grey,
  Median_PVE        = median_pve, # exluding grey module.
  Davies_Bouldin    = db, # DB index evaluates intra-cluster similarity and inter-cluster differences
  Calinski_Harabasz = ch, # Silhouette Index measure the distance between each data point, the centroid of the cluster it was assigned to and the closest centroid belonging to another cluster
  Silhouett_index   = sc)


print("Done!")





