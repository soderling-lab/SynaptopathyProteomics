#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(dendextend)
  library(getPPIs)
  library(anRichment)
  library(org.Mm.eg.db)
  library(utils)
  library(GOSemSim)
})

# User params.
save <- TRUE # Save (TRUE) or load (FALSE) the module GO data?

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein map.
protmap <- readRDS(file.path(rdatdir, "2_Prot_Map.RData"))

# Load network comparison results.
myfile <- list.files(rdatdir, pattern = "6171865", full.names = TRUE)
comparisons <- readRDS(myfile)

# Calculate GO similarity among gene products.
# Load GO data.
# Use Molecular function ontology.
msGO <- godata("org.Mm.eg.db", ont = "MF")

#-------------------------------------------------------------------------------
# Declare a function for module go semantic similarity.
#-------------------------------------------------------------------------------
# Function that:
# Loops through every resolution of WT or KO graph.
# Gets genes for each module.
# Calculate GO sim for those genes.
# Return list containing GO sim scores for every module at every resolution.
module_GOsemSim <- function(comparisons, myProts, myPartition) {
  # Empty list for output
  out <- list()
  for (i in 1:1) {
    message(paste("Working on resolution", i))
    # Get data for resolution i.
    resdat <- comparisons[[i]]
    # Split entrez ids by module partition.
    entrez <- protmap$entrez[match(names(resdat[[myProts]]), protmap$ids)]
    modules <- split(entrez, resdat[[myPartition]])
    # Remove module 0.
    modules <- modules[names(modules)[-1]]
    # For every module in the resolution,
    # compute its GOsemantic similarity.
    # This is a matrix of GO sem sim scores for every gene in the module.
    # summarize a module as the mean of its edges.
    module_gs <- sapply(modules, function(x) {
      mean(mgeneSim(entrez,
        semData = msGO,
        measure = "Wang",
        verbose = FALSE
      ))
    })
    # Return a list.
    names(module_gs) <- names(modules)
    out[[i]] <- module_gs
  } # ENDS LOOP.
  return(out)
}

#------------------------------------------------------------------------------
# Deploy the function.
#------------------------------------------------------------------------------

# For every resolution in WT graph, get module go semantic similarity.
#message("Working on WT partitions!")
#wtModuleGOsim <- module_GOsemSim(comparisons, "wtProts", "wtPartition")

# For every resolution in KO graph, get module go semantic similarity.
message("Working on KO partitions!")
koModuleGOsim <- module_GOsemSim(comparisons, "koProts", "koPartition")

# Save data.
#saveRDS(wtModuleGOsim, file.path(rdatdir, "3_WT_Module_GO_SemSim.RData"))
saveRDS(koModuleGOsim, file.path(rdatdir, "3_KO_Module_GO_SemSim.RData"))
