#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
## Prepare the workspace.
#------------------------------------------------------------------------------

# Directories.
here <- "D:/projects/SynaptopathyProteomics/analysis/3_Network-Analysis"
rootdir <- dirname(dirname(here))
datadir <- file.path(rootdir, "data")
rdatdir <- file.path(rootdir, "rdata")
tabsdir <- file.path(rootdir, "tables")
figsdir <- file.path(rootdir, "figures")
funcdir <- file.path(rootdir, "functions")

# Global options and imports.
suppressPackageStartupMessages({
  library(RCy3)
  library(getPPIs)
  })

#------------------------------------------------------------------------------
## Using the RCy3 module to interact with Cytoscape. 
#------------------------------------------------------------------------------

# Load the module subgraphs.
myfiles <- list.files(rdatdir,"Module_Graphs",full.names = TRUE)
wtGraphs <- readRDS(myfiles[2])
koGraphs <- readRDS(myfiles[1])

# Send to Cytoscape.
cytoscapePing()

# Need to add sigProt annotation...



# Send graph and node attributes to cytoscape.
for (i in 1:length(wtGraphs)) {
  namen <- paste0("WT",names(wtGraphs)[i])
  g <- wtGraphs[[i]]
  createNetworkFromIgraph(g,namen)
  # Apply FD layout.
  layoutNetwork(layout.name = "force-directed")
}

# Send graph and node attributes to cytoscape.
for (i in 1:length(koGraphs)) {
  namen <- paste0("KO",names(koGraphs)[i])
  g <- koGraphs[[i]]
  createNetworkFromIgraph(g,namen)
  # Apply FD layout.
  layoutNetwork(layout.name = "force-directed")
}
