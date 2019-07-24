#!/usr/bin/env python

## Performing hyperparameter optimization of the WGCNA function.

#------------------------------------------------------------------------------
# ## WGCNA Hyperparameter optimization.
#------------------------------------------------------------------------------

# Bayesian parameter optimization using the skopt-optimizer module.
# See: https://scikit-optimize.github.io/optimizer/index.html
from skopt.space import Real, Integer, Categorical
from skopt.utils import use_named_args
import subprocess

# Define WGCNA parameter space for optimization (17):
hyperparameters = [
        # Network construction arguments: correlation options
        Real(name = "maxPOutliers", low=0.0, high=1.0), #(0,1) pearson -> bicor
        # Basic tree cut options
        Integer(name ="deepSplit",       low = 1, high = 4),  #(0,4)
        Real(name    ="detectCutHeight", low=0.0, high=1.0),  #(0.995,0,1)
        Integer(name ="minModuleSize",   low = 2, high=3022), #(2,m)
        # Advanced tree cut options
        Real(name        ="maxCoreScatter",           low=0.0, high=1.0), #(None,0,1)
        Real(name        ="minGap",                   low=0.0, high=1.0), #(None,0,1)
        Real(name        ="minSplitHeight",           low=0.0, high=1.0), #(None,0,1)
        Real(name        ="minBranchEigennodeDissim", low=0.0, high=1.0), #(0.15,0,1)
        Real(name        ="minStabilityDissim",       low=0.0, high=1.0), #(None,0,1)
        Categorical(name ="useBranchEigennodeDissim", categories=[False,True]), #(False,True)
        Categorical(name ="pamStage",                 categories=[True,False]), #(True,False)
        Categorical(name ="pamRespectsDendro",        categories=[True,False]), #(True, False
        # Gene reassignment, module trimming, and module "significance" criteria
        Real(name        ="reassignThreshold", low=0.0, high=1.0), #(1e-6,0,1)
        Real(name        ="minCoreKME",        low=0.0, high=1.0), #(0,1)
        Real(name        ="minKMEtoStay",      low=0.0, high=1.0), #(0,1)
        Integer(name     ="minCoreKMESize",    low=0,   high=5),   #(0,minModuleSize)
        # Module merging options
        Real(name ="mergeCutHeight", low=0.0, high=1.0) #(0.15,0,1)
        ]

# Define a function, decorated with named arguments.
#@use_named_args(dimensions = hyperparameters)
def wgcna_evaluation(
        maxPOutliers,
        deepSplit,
        detectCutHeight,
        minModuleSize,
        maxCoreScatter,
        minGap,
        minSplitHeight,
        useBranchEigennodeDissim,
        minStabilityDissim,
        pamStage,
        pamRespectsDendro,
        reassignThreshold,
        minCoreKME,
        minCoreKMESize,
        minKMEtoStay,
        mergeCutHeight
        ):
    # Get function arguments and write these to file.
    args = locals()
    with open('parameters.txt', 'w') as f:
        for key in sorted(args):
            f.write(key + "\t" + str(args[key]) + "\n")
    f.close()
    # Call wgcna.r to perform WGCNA.
    cmd = ["./wgcna.r", "exprDat.Rds", "parameters.txt"]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    out = process.communicate()
    # Return output.
    return out
# EOF

#------------------------------------------------------------------------------
# Test the function.
#------------------------------------------------------------------------------

# None arguments will be replaced with defaults!
wgcna_evaluation(
        maxPOutliers = None,
        deepSplit    = None,
        detectCutHeight = None,
        minModuleSize   = None,
        maxCoreScatter  = None,
        minGap          = None,
        minSplitHeight = None, 
        useBranchEigennodeDissim = None,
        minStabilityDissim = None,
        pamStage = None,
        pamRespectsDendro = None,
        reassignThreshold = None,
        minCoreKME = None,
        minCoreKMESize = None,
        minKMEtoStay = None,
        mergeCutHeight = None
        )
print(out)

#------------------------------------------------------------------------------
# Perform optimization
#------------------------------------------------------------------------------

from skopt import gp_minimize
result = gp_minimize(func=wgcna_evaluation, dimensions=yperparameters,
        n_calls=20, base_estimator="ET", random_state=4)

print("Best fitness:", result.fun)
print("Best parameters:", result.x)

