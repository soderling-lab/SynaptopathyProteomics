#------------------------------------------------------------------------------
# ## wgcna-db.R
#------------------------------------------------------------------------------

# Wrapper function around WGCNA::BlockWiseModules() to perform WGCNA using
# user provided hyperparameters.If parameters are not provided, then the 
# defaults will be used.

wgcna <- function(exprDat, powerBeta = NULL, parameters=NULL, verbose = 0){

  ## Load WGCNA defaults stored in Bin/.
  # These will be used if the user does not provide any.
  source("D:/projects/Synaptopathy-Proteomics/Bin/wgcna-defaults.R")
  
  ## If provided, parse the user provided parameters.
  if (!exists("parameters")) {
    params <- default_params
  } else if (inherits(parameters,"list")) {
    user_params <- parameters
  } else {
    stop("please provide a list of parameters, or use the defaults.")
  }
    
  idx <- match(names(user_params), names(default_params))
  params <- default_params
  params[idx] <- user_params
  
  # Insure that params are the correct data type.
  # None type should be removed. What about factors?
  params <- lapply(params, function(x) type.convert(x))
  
  # Allow parallel WGCNA calculations if nThreads is >0.
  if (params$nThreads > 0) {
    allowWGCNAThreads(nThreads = params$nThreads)
    clusterLocal <- makeCluster(c(rep("localhost", params$nThreads)), type = "SOCK")
    registerDoParallel(clusterLocal)
  }
  
  # If specified, use the user's powerBeta. 
  if (inherits(params$power, "numeric")) {
    params$power <- powerBeta
    # If not specified, determine soft power to achieve scale free toplogy. 
  } else {
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
    checkMissingData      = as.logical(params$checkMissingData),
    
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
    minBranchEigennodeDissim = as.numeric(params$mergeCutHeight),
    stabilityLabels          = params$stabilityLabels,
    stabilityCriterion       = params$stabilityCriterion,
    minStabilityDissim       = params$minStabilityDissim,
    pamStage                 = as.logical(params$pamStage), 
    pamRespectsDendro        = as.logical(params$pamRespectsDendro),
    
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
  
  # Output:
  return(list("network" = net, "hyperparameters" = params))
} # END FUNCTION.

## END




