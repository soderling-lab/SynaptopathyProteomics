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

# The WGCNA wrapper:
source("D:/projects/Synaptopathy-Proteomics/Bin/wgcna.R")

#------------------------------------------------------------------------------
# ## Parse the command line input.
#------------------------------------------------------------------------------

# Load hyperparameters passed from Python as a dictionary of key : value pairs.
args = commandArgs(trailingOnly=TRUE)
args = gsub("[\r\n]", "", args)

file <- paste(dir, args, sep="/") # 
dict <- read.delim(file, header = FALSE, sep=",")

# Quick function to parse the python dictionary.
f <- function(x){
  y <- gsub("\\{|\\}| |'","", as.character(x))
  v <- unlist(strsplit(y,":"))
  return(v)
}

# Build a list of hyperparameters.
d <- apply(dict, 2, function(x) f(x))
keys <- d[1,]
values <- d[2,]
params <- as.list(values)
names(params) <- keys

#------------------------------------------------------------------------------
# ## Perform WGCNA
#------------------------------------------------------------------------------

# Load the normalized data.
file <- "D:/projects/Synaptopathy-Proteomics/RData/2_Combined_TAMPOR_cleanDat.Rds"
exprDat <- t(log2(readRDS(file)))

# WGCNA:
result <- wgcna(exprDat, hyperparameters=params, verbose = 0)
  
