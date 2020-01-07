#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

## User parameters to change.
# Which network/partitions to analyze?
net <- "Cortex" 
mypart <- c(Cortex = "10773682",Striatum = "10781799")[net] # relaxed criterion

# Global options and imports.
suppressPackageStartupMessages({
	library(org.Mm.eg.db)
	library(anRichment)
	library(data.table)
	library(getPPIs)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load expression data.
myfile <- file.path(rdatdir, paste0("3_",net,"_cleanDat.RData"))
data <- readRDS(myfile)

# Load network partitions--self-preservation enforced.
myfile <- list.files(rdatdir,pattern=mypart,full.names=TRUE)
partitions <- readRDS(myfile)

# What about disease enrichment...
list.files(rdatdir,pattern="DisGeneNet")
myfile <- file.path(rdatdir,"DisGeneNet_Curated_Variants_mouse.RData")
GOcollection <- readRDS(myfile)

# Mouse GO collection.
#if (!exists("musGO")) { GOcollection <- buildGOcollection(organism="mouse") }

# Perform GO analysis.
res = 93
GOresults <- moduleGOenrichment(partitions,res, protmap,GOcollection)

p = sapply(GOresults,function(x) min(x$pValue))
c(1:length(p))[p==min(p)]

sum(sapply(GOresults,function(x) any(x$FDR<0.1)))

write_excel(GOresults,"temp.xlsx")

