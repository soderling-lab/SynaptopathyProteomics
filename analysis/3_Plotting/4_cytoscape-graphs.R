#!/usr/bin/env Rscript

#' ---
#' title: Network Analysis
#' description: Protein co-expression network analysis
#' authors: Tyler W Bradshaw
#' ---

# User parameters to change:

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

# Parse command line input:
args <- commandArgs(trailingOnly = TRUE)
msg <- c("Please specify a tissue type to be analyzed:\n",
	 "       Choose either 'Cortex' or 'Striatum'.")
if (!length(args == 1)) { 
	stop(msg) 
} else { 
	type <- match(args[1],c("Cortex","Striatum"))
	tissue <- c("Cortex", "Striatum")[type]
	start <- Sys.time()
	message(paste("\nStarting analysis at:", start))
}

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global imports.
suppressPackageStartupMessages({
  library(RCy3)
  library(dplyr)
  library(igraph)
  library(data.table)
})

# Directories.
figsdir <- file.path(root, "figs")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
netwdir <- file.path(root, "networks")

# Load project's functions.
suppressWarnings({ devtools::load_all() })

# Load the data from root/data.
# Contains data, adjm, netw, and metadata.
dataset <- tolower(paste0(tissue,"_data"))
data(list=dataset)
eval(parse(text=paste0("data_list=",dataset)))

# Load the graph partition in root/data.
dataset <- tolower(paste0(tissue,"_partition"))
data(list=dataset)
eval(parse(text=paste0("partition=",dataset)))

# Load gene map.
data(gene_map)

# Load PPI adjm.
myfile <- file.path(rdatdir,"PPI_Adjm.csv")
ppin <- fread(myfile,drop=1) %>% as.matrix()
rownames(ppin) <- colnames(ppin)

#---------------------------------------------------------------------
## Prepare the data.
#---------------------------------------------------------------------

# Extract adjm and network from data_list.
adjm <- data_list$Adjm
netw <- data_list$Netw

# Insure that matrices are in matching order.
colNames <- colnames(adjm)
netw <- netw[colNames,colNames]
adjm <- netw[colNames,colNames]
ppin <- netw[colNames,colNames]

# Check:
check <- all(colnames(netw) == colnames(adjm) & colnames(adjm) == colnames(ppin))
if (!check) { stop() }

# Get modules from partition.
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))
modules <- modules[-which(names(modules)=="M0")]

# Create igraph graphs from adjacency matrices.
adjm_g <- graph_from_adjacency_matrix(adjm,mode="undirected",diag=FALSE,weighted=TRUE)
netw_g <- graph_from_adjacency_matrix(netw,mode="undirected",diag=FALSE,weighted=TRUE)
ppin_g <- graph_from_adjacency_matrix(ppin,mode="undirected",diag=FALSE)

# We have three graphs:
graphs <- list(adjm=adjm_g,
	       netw=netw_g,
	       ppin=ppin_g)

# Add gene symbols to graphs in list.
graphs <- lapply(graphs, function(x) {
			 symbols <- gene_map$gene[match(names(V(x)),gene_map$ids)]
		         x <- set_vertex_attr(x,"symbol",value = symbols)
		         return(x)
	         })

#--------------------------------------------------------------------
## Cytoscape
#--------------------------------------------------------------------

cytoscapePing()

## INPUTS
output_name = tissue # name of graph?
cysdir = netwdir # Where to save the graph?
network_layout = 'force-directed edgeAttribute=weight'

#winfile <- gsub("/mnt/d/","D:/",cysfile)
#openSession(winfile)

module <- modules[[1]]
module_name <- names(modules)[1]
nodes <- modules[[module_name]]


# When done, save cytoscape session.
#myfile <- file.path(netsdir,paste0(part_type,".cys"))
#winfile <- gsub("/mnt/d/","D:/",myfile)
#saveSession(winfile)
