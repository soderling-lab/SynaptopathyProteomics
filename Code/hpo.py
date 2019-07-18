#------------------------------------------------------------------------------
# ## Optimization of WGCNA Parameters.
#------------------------------------------------------------------------------

import subprocess

# Path to R executable (R.exe or Rscript.exe) must be a Linux path!
path2rexe = "/mnt/c/Program Files/R/R-3.6.1/bin/Rscript.exe" 

# Path to Rscript must be a Windows path!
path2rscript = "D:/projects/Synaptopathy-Proteomics/Code/test.R"

# Send command to R on Windows side!
cmd = [path2rexe, path2rscript]
x = subprocess.check_output(cmd, universal_newlines = True)

# Need WGCNA input parameters to R!
