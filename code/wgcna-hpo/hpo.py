#!/usr/bin/env python

#-----------------------------------------------------------------------------
# ## hpo.py
#-----------------------------------------------------------------------------

# Performs hyperparameter optimization of the WGCNA function.

#-----------------------------------------------------------------------------
# ## Define WGCNA parameters.
#-----------------------------------------------------------------------------

import json # for writing parameters to file. 

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
#with open('parameters.txt', 'w') as f:
#        print(hyperparameters, file = f)
with open('parameters.txt', 'w') as f:
    for key in?!?jedi=0,  sorted(hyperpara?!? (*_*s: AnyStr*_*) ?!?jedi?!?meters):
        f.write(key + "\t" + str(d[key]) + "\n")
f.close()

#------------------------------------------------------------------------------
# ## Optimization of WGCNA Parameters.
#------------------------------------------------------------------------------

import subprocess

# Path to R executable (R.exe or Rscript.exe). This must be a Linux path!
path2rexe = "/mnt/c/Program Files/R/R-3.6.1/bin/Rscript.exe" 

# Path to Rscript must be a Windows path!
rscript = "3b_wgcna-db.R"
path2rscript = "D:/projects/Synaptopathy-Proteomics/Code/" + rscript

# Define list of additional arguments to pass to rscript.
args = ["parameters.txt"]

# Send command to R on Windows side!
cmd = [path2rexe, path2rscript] + args
x = subprocess.check_output(cmd, universal_newlines = True)

#------------------------------------------------------------------------------
# ## Optimization of WGCNA Parameters (linux).
#------------------------------------------------------------------------------
subprocess.call(["/pathto/MyrScript.R"])

# Also try:
subprocess.run(["ls", "-l"])
