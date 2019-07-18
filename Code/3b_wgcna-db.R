#------------------------------------------------------------------------------
# ## wgcna-db.R
#------------------------------------------------------------------------------

# Performs wgcna given hyperparameters passed from python and calculates 
# clustering quality score.

# Windows Usage: [ Rscript.exe wgcna-db.R parameters.txt ] 
# Called by the python hpo.py script.

#------------------------------------------------------------------------------
# ## Prepare the workspace.
#------------------------------------------------------------------------------

## Global options and imports. 
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(WGCNA)
  library(doParallel)
  library(parallel)
})

# Load the WGCNA wrapper function:
source("D:/projects/Synaptopathy-Proteomics/Bin/wgcna.R")

#------------------------------------------------------------------------------
# ## Parse the command line input.
#------------------------------------------------------------------------------

# Load hyperparameters passed from Python as a dictionary of key : value pairs.
args = commandArgs(trailingOnly=TRUE)
args = gsub("[\r\n]", "", args)
#file <- "D:/projects/Synaptopathy-Proteomics/Code/parameters.txt"
file <- paste(dir, args, sep="/") # 
pydict <- read.delim(file, header = FALSE, sep=",")

# Quick function to parse the python dictionary.
f <- function(x){
  y <- gsub("\\{|\\}| |'","", as.character(x))
  v <- unlist(strsplit(y,":"))
  return(v)
}

# Build a list of hyperparameters.
dict <- apply(dict, 2, function(x) f(x))
keys <- dict[1,]
values <- dict[2,]
user_params <- as.list(values)
names(user_params) <- keys

# Remove NoneType from Python. These will be replaced by default = NULL.
user_params <- user_params[!user_params == "None"]

#------------------------------------------------------------------------------
# ## Perform WGCNA!
#------------------------------------------------------------------------------

# Load the normalized data.
file <- "D:/projects/Synaptopathy-Proteomics/RData/2_Combined_TAMPOR_cleanDat.Rds"
datExpr <- t(log2(readRDS(file)))

# WGCNA:
result <- wgcna(datExpr, powerBeta = NULL, parameters = user_params, verbose = 0)
  
