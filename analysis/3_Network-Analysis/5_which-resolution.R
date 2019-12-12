#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change.
net <- "Striatum" 
mypart <- c(Cortex = "10360847",Striatum = "10342568")[net]

# Global options and imports.
suppressPackageStartupMessages({
	library(org.Mm.eg.db)
	library(anRichment)
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

# Load network partitions-- self-preservation enforced.
myfile <- list.files(rdatdir,pattern=mypart,full.names=TRUE)
partitions <- readRDS(myfile)

#------------------------------------------------------------------------------
## Which resolution? Perform GO analysis of modules at every resolution.
#------------------------------------------------------------------------------

# Build mouse go collection:
message(paste("Preparing to analyze modules identified in the",
	       net,"network for GO enrichment..."))
if (!exists("musGOcollection")) {
	musGOcollection <- buildGOcollection(organism="mouse")
}

# Loop to perform GO enrichment for modules at every resolution.
message("Performing GO enrichment analysis...")
GOresults <- list()
for (i in seq_along(partitions)) {
  # Initialize progress bar.
  if (i == 1) {
    pb <- txtProgressBar(min = 0, max = length(partitions), style = 3)
  }
  # Perform GO analysis.
  GOresults[[i]] <- moduleGOenrichment(partitions,i, protmap,musGOcollection)
  # Update progress bar.
  setTxtProgressBar(pb, i)
  if (i == length(partitions)) {
    # Close pb.
	  close(pb)
	  message("\n")
  }
} # Ends loop.

# Save results.
myfile <- file.path(rdatdir,paste0("3_",net,"_Module_GO_Results.RData")) 
saveRDS(GOresults, myfile)

# Load results.
GOresults <- readRDS(myfile)

#------------------------------------------------------------------------------
## Examine GO results in order to define ~best resolution.
#------------------------------------------------------------------------------

# Examine biological enrichment of modules at every resolution.
# Summarize the biological significance of a resolution as the sum of 
# -log(GO pvalues) for all modules.
modSig <- lapply(GOresults, function(x) 
		 sapply(x, function(y) sum(-log(y$pValue))))

modSig <- lapply(GOresults, function(x) 
		 sapply(x, function(y) mean(-log(y$pValue))))

# The code above is confusing, this is what it does:
#x = results[[1]] # list of go enrichment for all modules at res 1.
#y = x[[1]] # go enrichment df of module 1.
#y$pValue
#-log(y$pValue)
#sum(-log(y$pValue))
#sapply(x,function(y) sum(-log(y$pValue)))
#lapply(results, function(x) sapply(x,function(y) sum(-log(y$pValue))))
#out = lapply(results, function(x) sapply(x,function(y) sum(-log(y$pValue))))

# Summarize every resolution.
resSum <- sapply(modSig, sum)
best_res <- c(1:length(resSum))[resSum == max(resSum)]
names(best_res) <- net
# Status report. 
message(paste("Best resolution based on GO enrichment:",best_res))

# Save results.
myfile <- file.path(rdatdir,paste0("3_",net,"_Best_Resolution.RData"))
saveRDS(best_res, myfile)
