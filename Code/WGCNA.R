#!/usr/bin/env Rscript

# Wrapper function around WGCNA::BlockWiseModules().
wgcna <- function(datExpr, powerBeta = "scale free", hyperparameters="Default", verbose = 0){
  
  ## Global options and imports. 
  options(stringsAsFactors = FALSE)
  suppressPackageStartupMessages({
    library(WGCNA)
    library(doParallel)
    library(parallel)
  })
  
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
  
  ## Parse the users input.
  
  # Check that hyperparameters is character or list.
  if (!(inherits(hyperparameters,"list") | inherits(hyperparameters, "character"))) {
    stop("please provide a named list of hyperaparameters, or use the defaults.")
    
    # If list of user defined parameters is defined, then use them. 
    } else if (inherits(hyperparameters, "list")) {
      msg <- paste("Using", length(hyperparameters), "user defined parameters!")
      if (verbose > 0) { print(msg)}
      idx <- match(names(hyperparameters), names(params))
      params <- params
      params[idx] <- hyperparameters
      
      # Else, use default parameters.
      } else if (tolower(hyperparameters) == "default") {
        msg <- c("Using default parameters!")
        if (verbose > 0) { print(msg) }
        params <- params # The defaults defined above. 
        
        # Otherwise, quit. 
        } else {
          stop("unable to parse the hyperparameters input.")
          }
  
  # Allow parallel WGCNA calculations if nThreads is >0.
  if (params$nThreads > 0) {
    allowWGCNAThreads(nThreads = params$nThreads)
    clusterLocal <- makeCluster(c(rep("localhost", params$nThreads)), type = "SOCK")
    registerDoParallel(clusterLocal)
    }
  
    # If specified, use the user's powerBeta. 
    if (inherits(powerBeta, "numeric")) {
    params$power <- powerBeta
    msg <- paste("Using user provided powerBeta:", powerBeta)
    if (verbose > 0) { print(msg) }
    
    # If not specified, determine soft power to achieve scale free toplogy. 
    } else if (tolower(powerBeta) == "scale free") {
      sft <- capture.output(
        WGCNA::pickSoftThreshold(datExpr,
                          weights     = params$weights,
                          powerVector = seq(2, 20, by = 1.0),
                          corFnc      = params$corType,
                          blockSize   = params$maxBlockSize,
                          networkType = params$networkType,
                          verbose     = params$verbose,
                          indent      = params$indent
        )
      )
      # Parse output of pickSoftThreshold().
      out <- sft[c((match("$fitIndices", sft)+1):(length(sft)-1))]
      rows <- unlist(strsplit(out,"\t"))
      header <- unlist(strsplit(rows[1], "\\s+"))[c(2:8)]
      fit <- do.call(rbind,strsplit(rows[2:length(rows)], "\\s+"))
      fit <- as.data.frame(apply(fit,2,function(x) as.numeric(x)))
      fit[,1] <- NULL
      colnames(fit) <- header
      # Calculate powerBeta
      powerBeta <- fit$Power[fit$SFT.R.sq > 0.8][1]
      params$power <- powerBeta
      fit <- subset(fit, fit$Power == powerBeta)
      msg <- paste0("Using a power of ", powerBeta, " to achieve a scale free fit of ", fit$SFT.R.sq,".")
      if (params$verbose > 0) { print(msg) }
      # else quit.    
      } else {
  stop("error parsing powerBeta.")
    }
  
  ## Perform WGCNA by calling the blockwiseModules() function. 
  net <- WGCNA::blockwiseModules(
    # Input data
    datExpr               = datExpr, 
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
    corType                   = params$corType,
    maxPOutliers              = params$maxPOutliers, 
    quickCor                  = params$quickCor,
    pearsonFallback           = params$pearsonFallback,
    cosineCorrelation         = params$cosineCorrelation,
    
    # Adjacency function options
    power                     = params$power,
    networkType               = params$networkType,
    replaceMissingAdjacencies = params$replaceMissingAdjacencies,
    
    # Topological overlap options
    TOMType                   = params$TOMType,
    TOMDenom                  = params$TOMDenom,
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
    minCoreKMESize     = params$minModuleSize/3,
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
  
  # Status message.
  if (params$verbose > 0) { print("Done!") }
  
  # Output:
  return(list("network" = net, "hyperparameters" = params))
} # END FUNCTION.

#------------------------------------------------------------------------------
# Perform WGCNA.

net <- wgcna(datExpr, powerBeta = "Scale free", hyperparameters, verbose = 0)  

  

