#!/usr/bin/env python3

## Performing baesian hyperparameter optimization of the WGCNA function using  
#  the skopt-optimizer module.
# Usage:
# ./wgcna-optimization.py data.Rds [acq_func]

#------------------------------------------------------------------------------
# ## Parse the command line input.
#------------------------------------------------------------------------------

from argparse import ArgumentParser

# Required input:
ap = ArgumentParser(description = ''' Perform optimization of the WGCNA 
        algorithm with Baesian hyperparameter optimization.''')
ap.add_argument("data", type = str, 
        help = ''' The normalized n x m expression
        data matrix that will be clustered by WGCNA.''')
# Optional arguments passed to gp_minimize:
ap.add_argument("-e", "--estimator", type = str,
        help = ''' The gp base estimator to use for the optimization.
        Default is Matern kernel.''', default = None)
ap.add_argument("-n", "--n_calls", type = int,
        help = ''' The total number of evaluations to be performed.''',
        default = 100)
ap.add_argument("-r", "--random_starts", type = int,
        help = ''' The number of random evaluations of wgcna with
        random parameters to be performed approximating it with
        'base_estimator'.''', default = 10)
ap.add_argument("-a", "--acq_func", type = str,
        help = ''' The acquisition function used by gp_minimize to 
        minimize over the posterior distribution. One of LCB, EI, PI, 
        gp_hedge, EIps, or PIps. ''', default = 'gp_hedge')
ap.add_argument("-o", "--optimizer", type = str,
        help = ''' The Method used to optimize the 'acq_function'.
        One of 'auto', 'sampling', or 'lbfgs'.''',
        default = 'auto')
ap.add_argument("-xi", "--xi", type = float,
        help = '''Use to set how much imporvement one wants over 
        previous evaluation used if acquisition function is 'EI',
        or 'PI'.''', default = 0.01)
ap.add_argument("-k", "--kappa", type = float,
        help = ''' How much variance in expected values. Used when
        acquisition function is LCB. Higher values favor exploration
        over exploitation.''', default = 1.96)
ap.add_argument("-s", "--seed", type = int,
        help = ''' For reproducible results, use something other than
        None.''', default = 4)
ap.add_argument("-v", "--verbose", type = bool,
        help = ''' If True then progress of function will be 
        printed to stdout. ''', default = True)

# Parse input arguments.
args = vars(ap.parse_args())
data = args['data']

#------------------------------------------------------------------------------
# ## Define WGCNA hyperparameters and an optimizer function.
#------------------------------------------------------------------------------

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
        "maxPOutliers"    : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : 0.5},
        # Basic tree cut options:
        "deepSplit"       : {'type' : Integer, 'low' : 0, 'high' : 4.0, 'default' : 2},
        "detectCutHeight" : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : 0.995},
        "minModuleSize"   : {'type' : Integer, 'low' : 2, 'high' : 302, 'default' : 12},   
        # Advanced tree cut options:
        "maxCoreScatter"     : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : 0.5},  
        "minGap"             : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : 0.5}, 
        "minSplitHeight"     : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : 0.5},                    
        "pamRespectsDendro"  : {'type' : Categorical, 'categories' : [False,True], 'default' : True}, 
        "minBranchEigennodeDissim" : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : 0.5}, 
        "minStabilityDissim"       : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : 0.5}, 
        "useBranchEigennodeDissim" : {'type' : Categorical, 'categories' : [False,True], 'default' : False}, 
        # Gene reassignment, module trimming, and module "significance" criteria
        "reassignThreshold"   : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : 1e-6}, 
        "minCoreKME"          : {'type' : Real, 'low' : 0.0, 'high' : 1.0, 'default' : 0.5},
        "minKMEtoStay"        : {'type' : Real, 'low' : 0.1, 'high' : 1.0, 'default' : 0.2},
        "minCoreKMESize"      : {'type' : Integer, 'low' : 1, 'high' : 100, 'default' : 4},
        # Module merging options
        "mergeCutHeight" : {'type' : Real, 'low' : 0, 'high' : 1, 'default' : 0.15}
        }
    
# If any parameters are 'None', then remove them. 
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
def wgcna_evaluation(**space):
    # Get function arguments with locals() and write these to file.
    args = locals()
    user_params = args['space']
    with open('inputparams.txt', 'w') as f:
        for key in sorted(user_params):
            f.write(key + "\t" + str(user_params[key]) + "\n")
    f.close()
    # Call wgcna.r to perform WGCNA and evaluate quality of partition..
    cmd = ["./wgcna.r", data, "--parameters", "inputparams.txt"]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    out = process.communicate()
    # Parse output of wgcna.r
    score = float(out[0].decode().split()[1])
    # Return quality score.
    return score
# EOF

#------------------------------------------------------------------------------
# Perform baesian optimizahyerparameterotion of the WGCNA function with gp_minimize.
#------------------------------------------------------------------------------
# FIXME: It might be helpful if arguments were saved to file.

from skopt import gp_minimize

result = gp_minimize(func = wgcna_evaluation, dimensions = space,
        base_estimator  = args['estimator'],        
        n_calls         = args['n_calls'],               
        n_random_starts = args['random_starts'],         
        acq_func        = args['acq_func'], 
        acq_optimizer   = args['optimizer'],       
        x0              = defaults,        
        xi              = args['xi'],                    
        kappa           = args['kappa'],                
        random_state    = args['seed'],            
        verbose         = args['verbose']            
        )

#------------------------------------------------------------------------------
# Save the results and clean-up. 
#------------------------------------------------------------------------------

# Save the search results.
import pandas as pd

params = [param.name for param in space]
df = pd.DataFrame(data = result.x_iters, columns = params) 
df['Quality'] = result.func_vals
df.to_csv('search_space.csv', index=False)

# Save the best parameters.
values = result.x
optimized_params = dict(zip(params,values))
out_file = "optimized_params.txt"
with open(out_file, 'w') as f:
    for key in sorted(optimized_params):
        f.write(key + "\t" + str(optimized_params[key]) + "\n")
    f.close()

# Clean-up.
import os
os.remove("parameters.txt")
