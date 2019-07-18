# Testing the WGCNA function.

# Prepare the workspace.
rm(list = ls())
if (.Device != "null device") dev.off()
cat("\f")
options(stringsAsFactors = FALSE)

dir <- "D:/projects/Synaptopathy-Proteomics"
setwd(dir)
source(paste(dir,"Bin","wgcna.R",sep="/"))

# Load the data.
file <- paste(dir,"RData","2_Combined_TAMPOR_cleanDat.Rds", sep ="/")
datExpr <- t(log2(readRDS(file)))

# Perform WGCNA.
hyperparameters <- list(
  
)
wgcna_results <- wgcna(datExpr, verbose =  1)

# Saving or returning TOM. Will this be faster?
params$getTOMs         <- NULL
params$saveTOMs        <- FALSE 
params$saveTOMFileBase <- "blockwiseTOM"

