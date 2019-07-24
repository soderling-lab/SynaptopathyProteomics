## Default WGCNA Parameters

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

## Define defaults. Ignore NULL defaults.

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
#default_params$power                         <- NULL
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
default_params$minCoreKMESize    <- (min(20, ncol(exprDat)/2))/3
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
