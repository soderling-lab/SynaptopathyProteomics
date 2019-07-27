#!/usr/bin/env Rscript

## wgcna.R [exprDat.Rds] [parameters.txt]

#------------------------------------------------------------------------------
# ## Parse the command line input.
#------------------------------------------------------------------------------

#FIXME: save wgcna output.
# Need to account for percent grey in quality score.

# Global options
options(stringsAsFactors = FALSE)

# If provided, parse the user provided hyperparemeters.
args <- commandArgs(trailingOnly=TRUE)
nargs <- length(args)

if (nargs == 0) { stop("Please provide input expression data!")
  } else if (nargs == 1) {
    # Load the expression data.
	  dir <- getwd()
	  data_file <- paste(dir, args[1], sep="/")
	  exprDat <- log2(t(readRDS(data_file)))
	  user_params <- NULL
	  } else if (nargs == 2) {
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
      } else {
        stop("Unable to parse user input.")
      }

###############################################################################
# Load data for testing if no command line arguments passed.
if (nargs == 0){
  dir <- "D:/projects/Synaptopathy-Proteomics/code/wgcna-hpo"
  setwd(dir)
  data_file  <- paste(dir, "exprDat.Rds", sep="/")
  params_file <- paste(dir, "parameters.txt", sep="/")
  exprDat <- log2(t(readRDS(data_file)))
  temp_params <- read.delim(params_file, header = FALSE, col.names = c("Parameter","Value"))
  temp_params$Value[temp_params$Value == "True"] <- TRUE
  temp_params$Value[temp_params$Value == "False"] <- FALSE
  temp_params <- temp_params[!(temp_params$Value == "None"),]
  user_params <- as.list(temp_params$Value)
  names(user_params) <- temp_params$Parameter
}
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
  suppressPackageStartupMessages({
    require(WGCNA)
    require(doParallel)
    require(parallel)
  })
  
  ## Load WGCNA defaults stored in Bin/.
  # These will be used if the user does not provide any.
  source("./defaults.r")
  
  ## If provided, parse the user provided parameters.
  if (!exists("parameters") | length(parameters) == 0) {
	  message("No user defined parameters. Using the saved default parameters!")
	  params <- default_params
  } else if (inherits(parameters,"list") & length(parameters) > 0) {
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
	  message("Picking the best soft thresholding power to achieve scale free fit > 0.8!")
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
    
    # Calculate the best power_beta.
    params$power <- sft$powerEstimate
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
    return(list("data" = exprDat, "network" = net, "hyperparameters" = params))
    } 
# EOF

#------------------------------------------------------------------------------
# ## Perform WGCNA!
#------------------------------------------------------------------------------

results <- wgcna(exprDat, parameters = user_params)

# Extract WGCNA results.
exprDat <- results$data
net <- results$network
params <- results$hyperparameters

# The number of modules and percent grey.
nmodules <- length(unique(results$network$colors)) - 1 # exclude grey

# Catch errors caused by 0 or 1 clusters.
if (nmodules < 2){
  message("Cannot compute quality statistics for 0 or 1 modules!")
  print(10) # <- A bad quality score
  quit()
}

msg <- paste("WGCNA identified", nmodules, "modules!")
message(msg)

# Calculate percent grey (unclustered nodes):
is_grey <- net$colors == "grey"
percent_grey <- sum(is_grey)/length(is_grey)

# Calculate median percent variance explained, module cohesivness:
pve <- WGCNA::propVarExplained(exprDat, net$colors, net$MEs, corFnc = params$corType)
median_pve <- median(pve[!names(pve) == "PVEgrey"])

# Remove grey nodes from data.
subDat <- exprDat[,!is_grey]

#------------------------------------------------------------------------------
# ## Evaluate the quality of the WGCNA partition.
#------------------------------------------------------------------------------
# If using WSL, you need to initialize x11 display. Before calling clusterSim,
# insure you have launched x11 display first! 

# FIXME: catch this error and tell the user to use bin/xlaunch
# Warning messages:                                                                                  1: In rgl.init(initValue, onlyNULL) : RGL: unable to open X11 display                              2: 'rgl_init' failed, running with rgl.useNULL = TRUE

suppressPackageStartupMessages({
  library(clusterSim)
})

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
cl <- as.integer(as.numeric(as.factor(net$colors[!is_grey])))
db <- index.DB(adjm, cl)$DB
ch <- index.G1(adjm, cl)
sc <- index.S(as.dist(diss), cl)

# quality results
quality <- list(
  n_Modules         = nmodules,
  Percent_grey      = percent_grey,
  Median_PVE        = median_pve, # exluding grey module.
  Davies_Bouldin    = db,         # Smaller is better, min(0): DB index evaluates intra-cluster similarity and inter-cluster differences
  Calinski_Harabasz = ch,         # CH index: larger is better (ratio of between cluster variance and within cluster variance)
  Silhouette_Coeff   = sc)        # S Index measure the distance between each data point, the centroid of the cluster it was assigned to and the closest centroid belonging to another cluster

# Which metric to use?
metric = "Davies_Bouldin"

if (metric == "Calinski_Harabasz") {
  score = 1/ch
} else if (metric == "Davies_Bouldin") {
  score = db
} else if (metric == "Silhouette_Coeff") {
  score = -1*sc
}

# Return score
print(score * percent_grey)

# END
