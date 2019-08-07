#!/usr/bin/env Rscript


#------------------------------------------------------------------------------
# ## Parse the command line arguments.
#------------------------------------------------------------------------------

# Global options:
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  require(argparser, quietly = TRUE)
})

# Check for input data before proceeding.
if(length(commandArgs(trailingOnly=TRUE)) == 0) {
  msg <- paste("Please provide normalized n x m matrix of expression data as input for WGCNA!",
               "Use ./wgcna.r --help for help.", sep = "\n")
  stop(msg)
}

# Parse the command line arguments.
p <- arg_parser("Perform WGCNA given a normalized n x m matrix of protein or gene expression data. ")
p <- add_argument(p, "data",
                  help = paste("normalized n (sample) x m (gene) expression matrix",
                               "provided as a .Rds file in the same directory as this script"),
                  default = NULL)
args <- parse_args(p)

# Load data as n x m normalized expression matrix.
dir <- getwd()
project_dir <- dirname(dirname(dir))
data_file <- paste(dir, args$data, sep="/")
exprDat <- readRDS(data_file)

# Load custome functions.
my_functions <- paste(project_dir, "code", "0_Functions", "0_Functions.R", sep = "/")
source(my_functions)

#------------------------------------------------------------------------------
# ## Build a command to be passed to radalib.
#------------------------------------------------------------------------------

suppressPackageStartupMessages({
  require(WGCNA, quietly = TRUE)
  require(data.table, quietly = TRUE)
})

# Calculate WS adjacency matrix.
adjm <- silently(bicor, exprDat)

# Write network in pajek format to file.
script_dir <- paste(project_dir,"bin","radalib", sep="/")
network_file <- paste(script_dir,"network_WS.net", sep = "/")
write.pajek(adjm, network_file)

# Build a command to be passed to radalib.
script <- "./communities_detection.exe"
type <- "WS" # Weighted signed network.
output_file <- paste(script_dir,"best_partition.txt", sep="/")
cmd <- paste(script, "s", type, "l" , "10",basename(network_file), basename(output_file))


#------------------------------------------------------------------------------
# ## Community detection.
#------------------------------------------------------------------------------

# Really slow!!!
message("Searching for communities...")
setwd(script_dir)
result <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
setwd(dir)

print(result)
quit()

# ENDOFILE
#------------------------------------------------------------------------------
