#!/usr/bin/env Rscript

# Load the normalized expression data.
here <- getwd()
root <- dirname(dirname(here))
datadir <- file.path(root,"data")

# Load misc functions.
library(TBmiscr)
library(data.table)

# Load the normalized protein abundance data.
myfile <- file.path(datadir,"2_Combined_TAMPOR_cleanDat.Rds")
cleanDat <- readRDS(myfile)

# Load sample traits.
myfile <- file.path(datadir,"2_Combined_traits.Rds")
traits <- readRDS(myfile)

# Remove QC data and log2 transform.
out <- traits$SampleType[match(colnames(cleanDat),rownames(traits))] == "QC"
data <- log2(cleanDat[,!out])

# Create signed adjacency matrix, coerce to data.table.
adjm <- data.table(silently(WGCNA::bicor, t(data)))
rownames(adjm) <- colnames(adjm)

## Pass this adjacency matrix as .csv to leidenalg-clustering.py script.
# Write adjm to .csv.
tmpfile <- "adjm.csv"
fwrite(adjm, tmpfile, row.names = TRUE) 

# Create system command.
script <- "leidenalg-clustering.py"
cmd <- paste(file.path(here,script), tmpfile)

# Call leidenalg.py and remove temporary file.
system(cmd, intern = TRUE, ignore.stderr = TRUE)

# Remove temporary file.
unlink(tmpfile)

#------------------------------------------------------------------------------

# ENDOFILE
#------------------------------------------------------------------------------
