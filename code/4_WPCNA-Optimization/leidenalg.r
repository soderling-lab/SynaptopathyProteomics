#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# Examine the optimized partition of the WGCNA graph.
rootdir <- dirname(dirname(getwd))

# Load partition profile.
file <- paste(rootdir,tables,"wpcgraph_partition_profile.csv", sep = "/")
file <- paste(rootdir,tables,"ppi_partition_profile.csv", sep = "/")

# Examine relationship between:
# * Resolution and modularity
# * Resolution and number of clusters (k)

# Are their plateaus in this graph that suggest which resolution is ~best or
# Which resolutions are more stable? Which resolution optimizes modularity?

# Compare partition of WT graph with KO graph. Are any modules preserved or 
# are any modules not preserved (i.e. different)? 

# How can we better address which clusters are different? 


