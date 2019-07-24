#!/usr/bin/env python

## Performing hyperparameter optimization of the WGCNA function.

#-----------------------------------------------------------------------------
# ## Define WGCNA parameters.
#-----------------------------------------------------------------------------

# Parameters to be optimized (17):
# Please specify all parameters that you want to optimized or don't include them.
# Don't use None! This won't be interpreted correctly by R.
hyperparameters = {
        # Network construction arguments: correlation options
        "maxPOutliers"             : 1,     #(0,1) pearson -> bicor behavior
        # Basic tree cut options
        "deepSplit"                : 2,     #(0,4)
        "detectCutHeight"          : 0.995, #(0,1)
        "minModuleSize"            : 2,     #(2,m) | m = Total number of genes.
        # Advanced tree cut options
        "maxCoreScatter"           : 0,     #(None,0,1)
        "minGap"                   : 0,     #(None,0,1)
        "minSplitHeight"           : 0,     #(None,0,1)
        "useBranchEigennodeDissim" : False, #(True,False)
        "minBranchEigennodeDissim" : 0.15,  #(0,1)
        "minStabilityDissim"       : 0,     #(None,0,1)
        "pamStage"                 : True,  # TRUE,FALSE
        "pamRespectsDendro"        : True,  # True, False
        # Gene reassignment, module trimming, and module "significance" criteria
        "reassignThreshold"        : 1e-6,  #(0,1)
        "minCoreKME"               : 0.5,   # (0,1)
        "minCoreKMESize"           : 5,     #(0,minModuleSize)
        "minKMEtoStay"             : 0.3,   #(0,1)
        # Module merging options
        "mergeCutHeight"           : 0.15   #(0,1)
        }

# Save parameters.txt. These will be passed as input to the WGCNA function.
import json  
with open('parameters.txt', 'w') as f:
    for key in sorted(hyperparameters):
        f.write(key + "\t" + str(hyperparameters[key]) + "\n")
f.close()

#------------------------------------------------------------------------------
# ## Optimization of WGCNA Parameters.
#------------------------------------------------------------------------------

# Call wgcna.r to perform WGCNA and evaluate partition quality.
import subprocess
process = subprocess.Popen(["./wgcna.r", "exprDat.Rds", "parameters.txt"], stdout=subprocess.PIPE)
out = process.communicate()

# Check output.
print(out)

#------------------------------------------------------------------------------
# Testing the skopt module.

## Using multiple named arguments!
from skopt.space import Real

# Build a list of parameters.
params = [
        Real(name='foo', low=0.0, high=1.0),
        Real(name='bar', low=0.0, high=1.0),
        Real(name='baz', low=0.0, high=1.0),
        ]

# Define a function.
from skopt.utils import use_named_args
@use_named_args(dimensions=params)
def f(foo, bar, baz):
    out = foo ** 2 + bar ** 4 + baz ** 8
    return out

# Perform optimization
from skopt import forest_minimize
result = forest_minimize(func=f, dimensions=params,
        n_calls=20, base_estimator="ET", random_state=4)

print("Best fitness:", result.fun)
print("Best parameters:", result.x)


#------------------------------------------------------------------------------
## Hyperparameter optimization.

from skopt.space import Real, Integer, Categorical

# Define parameter space for WGCNA algorithm.
hyperparameters = [
        # Network construction arguments: correlation options
        Real(name        = "maxPOutliers", low=0.0, high=1.0), #(0,1) pearson -> bicor
        # Basic tree cut options
        Integer(name     ="deepSplit", low = 1, high = 4),   #(0,4)
        Real(name        ="detectCutHeight", low=0, high=1.0),  #(0.995,0,1)
        Integer(name     ="minModuleSize", low=2, high=3022),     #(2,m)
        # Advanced tree cut options
        Real(name        ="maxCoreScatter", low=0.0, high=1.0),            #(None,0,1)
        Real(name        ="minGap", low=0.0, high=1.0),                    #(None,0,1)
        Real(name        ="minSplitHeight", low=0.0, high=1.0),            #(None,0,1)
        Categorical(name ="useBranchEigennodeDissim", categories=[False,True]), #(False,True)
        Real(name        ="minBranchEigennodeDissim", low=0.0, high=1.0),   #(0.15,0,1)
        Real(name        ="minStabilityDissim",low=0.0, high=1.0),         #(None,0,1)
        Categorical(name ="pamStage", categories=[True,False]),               #(TRUE,FALSE)
        Categorical(name ="pamRespectsDendro", categories=[True,False]),       #(True, False
        # Gene reassignment, module trimming, and module "significance" criteria
        Real(name        ="reassignThreshold", low=0.0, high=1.0),           #(1e-6,0,1)
        Real(name        ="minCoreKME", low=0.0, high=1.0),                  # (0,1)
        Integer(name     ="minCoreKMESize", low=0, high=5),     #(0,minModuleSize)
        Real(name        ="minKMEtoStay", low=0.0, high=1.0),   #(0,1)
        # Module merging options
        Real(name       ="mergeCutHeight", low=0.0, high=1.0) #(0.15,0,1)
        ]

# Define a function, decorated with named arguments.
from skopt.utils import use_named_args

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
    return None


# Call wgcna.r to perform WGCNA.
import subprocess
cmd = ["./wgcna.r", "exprDat.Rds", "parameters.txt"]
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
out = process.communicate()

# Perform optimization
from skopt import gp_minimize
result = gp_minimize(func=wgcna_evaluation, dimensions=yperparameters,
        n_calls=20, base_estimator="ET", random_state=4)

#######
# Test the function.

# None arguments will be replaced with defaults!
wgcna_evaluation(
        maxPOutliers = None,
        deepSplit = None,
        detectCutHeight = None,
        minModuleSize = None,
        maxCoreScatter = None,
        minGap = None,
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
