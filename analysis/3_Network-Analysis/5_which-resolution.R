#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change.
net <- "Cortex" 
mypart <- c(Cortex = "10360847") # c(Striatum = "10342568")

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
musGOcollection <- buildGOcollection(organism="mouse")

# Loop to perform GO enrichment for modules at every resolution.
message((paste("Analyzing GO enrichment of modules identified in",
	       net,"network..."))

results <- list()
for (i in seq_along(partitions)) {
  # Initialize progress bar.
  if (i == 1) {
    pb <- txtProgressBar(min = 0, max = length(partitions), style = 3)
  }
  # Perform GO analysis.
  results[[i]] <- moduleGOenrichment(partitions,i, protmap,musGOcollection)
  # Update progress bar.
  setTxtProgressBar(pb, i)
  if (i == length(partitions)) {
    # Close pb.
    close(pb)
    message("GO enrichment analysis complete!")
  }
} # Ends loop.

# Save results.
myfile <- file.path(rdatdir,paste0("3_",net,"Module_GO_Results.RData")
saveRDS(results, myfile)

#------------------------------------------------------------------------------
## Examine GO results in order to define ~best resolution.
#------------------------------------------------------------------------------

# Remove M0 results.
moduleGO <- lapply(moduleGO, function(x) x[-grep("M0", names(x))])

# Examine biological enrichment of modules at every resolution.
# Summarize the biological significance of a resolution as the sum of 
# -log(GO pvalues) for all modules.
modSig <- lapply(moduleGO, function(x) 
		 sapply(x, function(y) sum(-log(y$pValue))))

# Remove resolutions that have no sig GO terms.
out <- c(1:length(modSig))[sapply(modSig,function(x) length(x)==0)]
modSig <- modSig[-out]

# Summarize every resolution.
resSum <- sapply(modSig, sum)
best_res <- c(1:length(resSum))[resSum == max(resSum)]
names(best_res) <- net

# Status report. 
message(paste("Best resolution based on GO enrichment:",best_res))

# Save results.
myfile <- file.path(rdatdir,paste0("3_",net,"best_resolution.RData")
saveRDS(best_res, myfile)
