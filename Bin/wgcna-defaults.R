#!/usr/bin/env Rscript

## Define Default WGCNA Hyperparameters.
params <- list(); {
  # Input data
  #params$datExpr <- datExpr
  params$weights <- NULL
  
  # Data checking options
  params$checkMissingData <- TRUE
  
  # Options for splitting data into blocks
  params$blocks                <- NULL
  params$maxBlockSize          <- 5000
  params$blockSizePenaltyPower <- 5
  params$nPreclusteringCenters <- as.integer(min(ncol(datExpr)/20, 100*ncol(datExpr)/params$maxBlockSize))
  params$randomSeed            <- 12345
  
  # load TOM from previously saved file?
  params$loadTOM <- FALSE
  
  # Network construction arguments: correlation options
  params$corType           <- "bicor"
  params$maxPOutliers      <- 1
  params$quickCor          <- 0
  params$pearsonFallback   <- "individual"
  params$cosineCorrelation <- FALSE
  
  # Adjacency function options
  params$power                         <- NULL
  params$networkType                   <- "unsigned"
  params$replaceMissingAdjacencies     <- FALSE
  params$suppressTOMForZeroAdjacencies <- FALSE
  
  # Topological overlap options
  params$TOMType  <- "signed"  # c()
  params$TOMDenom <- "min"     # c("min","mean")
  
  # Saving or returning TOM
  params$getTOMs         <- NULL
  params$saveTOMs        <- FALSE 
  params$saveTOMFileBase <- "blockwiseTOM"
  
  # Basic tree cut options
  params$deepSplit       <- 2
  params$detectCutHeight <- 0.995 
  params$minModuleSize   <- min(20, ncol(datExpr)/2)
  
  # Advanced tree cut options
  params$maxCoreScatter           <- NULL 
  params$minGap                   <- NULL
  params$maxAbsCoreScatter        <- NULL
  params$minAbsGap                <- NULL
  params$minSplitHeight           <- NULL
  params$minAbsSplitHeight        <- NULL
  params$useBranchEigennodeDissim <- FALSE
  params$minBranchEigennodeDissim <- 0.15
  params$stabilityLabels          <- NULL
  params$stabilityCriterion       <- "Individual fraction" # c("Individual fraction", "Common fraction")
  params$minStabilityDissim       <- NULL
  params$pamStage                 <- TRUE
  params$pamRespectsDendro        <- TRUE
  
  # Gene reassignment, module trimming, and module "significance" criteria
  params$reassignThreshold <- 1e-6
  params$minCoreKME        <- 0.5
  params$minCoreKMESize    <- (min(20, ncol(datExpr)/2))/3
  params$minKMEtoStay      <- 0.3
  
  # Module merging options
  params$mergeCutHeight <- 0.15
  params$impute         <- TRUE
  params$trapErrors     <- FALSE
  
  # Output options
  params$numericLabels <- FALSE
  
  # Options controlling behaviour
  params$nThreads                 <- 0
  params$useInternalMatrixAlgebra <- FALSE
  params$useCorOptionsThroughout  <- TRUE
  params$verbose                  <- verbose 
  params$indent                   <- 0
}