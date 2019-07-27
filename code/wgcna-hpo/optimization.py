#!/usr/bin/env python3

## Performing hyperparameter optimization of the WGCNA function.

#------------------------------------------------------------------------------
# ## WGCNA Hyperparameter optimization.
#------------------------------------------------------------------------------

# Bayesian parameter optimization using the skopt-optimizer module.
# The DepractionWarning which results from importing skopt is acknowledged by 
# skopt's developers and will be fixed in the next release. It cannot be 
# suppressed by `warnings.filterwarnings("ignore", category=DeprecationWarning)`.
# Just ignore it.

from skopt.space import Real, Integer, Categorical
from skopt.utils import use_named_args
    
import subprocess

# Define WGCNA hyperparameters for optimization:
# Use 'default' : None for parameters you don't want to optimize.
hyperparameters = {
        # Network construction arguments: correlation options:
        "maxPOutliers"    : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : None},
        # Basic tree cut options:
        "deepSplit"       : {'type' : Integer, 'low' : 0, 'high' : 4.0, 'default' : 2},
        "detectCutHeight" : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : 0.995},
        "minModuleSize"   : {'type' : Integer, 'low' : 2, 'high' : 302, 'default' : 12},   
        # Advanced tree cut options:
        "maxCoreScatter"     : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : None},  
        "minGap"             : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : None}, 
        "minSplitHeight"     : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : None},                    
        "pamRespectsDendro"  : {'type' : Categorical, 'categories' : [False,True], 'default' : True}, 
        "minBranchEigennodeDissim" : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : None}, 
        "minStabilityDissim"       : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : None}, 
        "useBranchEigennodeDissim" : {'type' : Categorical, 'categories' : [False,True], 'default' : False}, 
        # Gene reassignment, module trimming, and module "significance" criteria
        "reassignThreshold"   : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : 1e-6}, 
        "minCoreKME"          : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : 0.5},
        "minKMEtoStay"        : {'type' : Real, 'low' : 0.1, 'high' : 1.0, 'default' : 0.2},
        "minCoreKMESize"      : {'type' : Integer, 'low' : 1, 'high' : 100, 'default' : 4},
        # Module merging options
        "mergeCutHeight" : {'type' : Real, 'low' : 0, 'high' : 1, 'default' : 0.15}
        }
    
# Remove parameters with default of 'None'. 
out = [key for key in hyperparameters if hyperparameters.get(key).get('default') is None] 
for key in out: del hyperparameters[key]

# Loop through hyperparameter dict and create list which will be passed to skopt-optimize. 
# use if else to handle categories and integers/real numbers.
space = list()
for key in hyperparameters.keys():
    param = key
    values = hyperparameters.get(key)
    fun = hyperparameters.get(key).get('type')
    if fun is Categorical:
        space.append(fun(name = param, categories = values.get('categories')))
    else:
        lower_bound = values.get('low')
        upper_bound = values.get('high')
        space.append(fun(name = param, low = lower_bound, high = upper_bound))

# Parse the hyperparameter defaults as a starting point for the function.
defaults = [hyperparameters.get(key).get('default') for key in hyperparameters]

# Define a function, decorated with named arguments.
# NOTE: parameters need to be defined in the same order as they are above.
@use_named_args(dimensions = space)
def wgcna_evaluation(
        deepSplit,
        detectCutHeight,
        minModuleSize,
        pamRespectsDendro,
        useBranchEigennodeDissim,
        reassignThreshold,
        minCoreKME,
        minKMEtoStay,
        minCoreKMESize,
        mergeCutHeight
        ):
    # Get function arguments and write these to file.
    args = locals()
    with open('parameters.txt', 'w') as f:
        for key in sorted(args):
            f.write(key + "\t" + str(args[key]) + "\n")
    f.close()
    # Call wgcna.r to perform WGCNA and evaluate quality of partition..
    cmd = ["./wgcna.r", "exprDat.Rds", "parameters.txt"]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    out = process.communicate()
    # Parse output of wgcna.r
    score = float(out[0].decode().split()[1])
    # Return quality score.
    return score
# EOF

#------------------------------------------------------------------------------
# Perform optimization.
#------------------------------------------------------------------------------

from skopt import gp_minimize

result = gp_minimize(func = wgcna_evaluation, dimensions = space,
        base_estimator="ET", 
        n_calls=100,         # total number of evaluations
        n_random_starts=10,  # Number of random starts before the approximating the function with base_estimator.
        acq_func='LCB', # function to minimize c(LCB, EI, PI, gp_hedge) gp_hedge is a probabilistic combination of LCB, EI, and PI.  
        acq_optimizer='auto',# .
        x0=defaults, # If provided then f(x0) is evaluated, followed by n_random_starts. Finally, n_calls - len(x0) - n_random_starts are evaluated.
        kappa=1.96,
        random_state=4,      # For reproducible results use something other than None.
        verbose=True,
        )

print("Best fitness:", result.fun)
print("Best parameters:", result.x)



