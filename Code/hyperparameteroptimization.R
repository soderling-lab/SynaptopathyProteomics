#------------------------------------------------------------------------------
# Clean up.
rm(list = ls())
if (.Device != "null device" ) dev.off()
cat("\f")

#------------------------------------------------------------------------------
# Trying to figure our bayesianhyperparameter optimization...

library(WGCNA)
options(stringsAsFactors = FALSE)

# Generate some test data.
set.seed(2)
nGenes <- 1000
nSamples <- 8
data <- matrix(rnorm(n = nSamples * nGenes, mean = 100), nSamples, nGenes)

# WGCNA function to be optimized:
score_network <- function(data, hyperparameters){
  net <- blockwiseModules(datExpr = data,
                          # Hyperparameters:
                          deepSplit      = hyperparameters$deepSplit,
                          minModuleSize  = hyperparameters$minModuleSize,
                          mergeCutHeight = hyperparameters$mergeCutHeight,
                          reassignThresh = hyperparameters$reassignThresh,
                          minCoreKMESize = hyperparameters$minCoreKMESize,
                          minKMEtoStay   = hyperparameters$minKMEtoStay,
                          # Defaults:
                          power             = 5,
                          TOMDenom          = "mean",
                          detectCutHeight   = 0.995,
                          corType           = "bicor",
                          networkType       = "signed",
                          pamStage          = TRUE,
                          pamRespectsDendro = TRUE,
                          saveTOMs          = FALSE,
                          maxBlockSize      = 12000,
                          verbose           = 0)
  
  # Calculate and return median module coherence. 
  pve <- propVarExplained(datExpr = data, colors = net$colors, MEs = net$MEs)
  return(list(Score = median(pve)))
}
  
#------------------------------------------------------------------------------3
# Test the function.

# Combine into a list.
hyperparameters <- list(
  minModuleSize = 34,
  deepSplit = 4,
  mergeCutHeight = 0.18069189726375,
  reassignThresh = 0.0173420400521718,
  minKMEtoStay = 0.40141052021645,
  minCoreKMESize = 4)
  
test <- score_network(data, hyperparameters)
test # 0.482824

#------------------------------------------------------------------------------

library(ParBayesianOptimization)

hyperparameters = list(
  minModuleSize = c(3L, 50L),  
  deepSplit =  c(0L, 4L),
  mergeCutHeight = c(0.01, 0.2),
  reassignThresh = c(0.01, 0.1),
  minKMEtoStay = c(0.1,0.7), 
  minCoreKMESize = c(3L,15L)
)

hpo <- BayesianOptimization(score_network, 
                     bounds = hyperparameters, 
                     saveIntermediate = NULL,
                     leftOff = NULL, 
                     parallel = FALSE, 
                     packages = "WGCNA", 
                     export = "data",
                     initialize = TRUE, 
                     initGrid = NULL, 
                     initPoints = 1, 
                     bulkNew = parallelThreads,
                     nIters = 100, 
                     kern = "Matern52", 
                     beta = 0, 
                     acq = "ucb",
                     stopImpatient = list(newAcq = "ucb", rounds = Inf), 
                     kappa = 2.576,
                     eps = 0, 
                     gsPoints = 100, 
                     convThresh = 1e+07,
                     minClusterUtility = NULL, 
                     noiseAdd = 0.25, 
                     verbose = 1)



  
  
  