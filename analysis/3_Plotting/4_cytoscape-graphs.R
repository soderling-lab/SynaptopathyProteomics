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

##FIXME: is adjm enhanced network?
max(E(g)$weight) 

# Parse command line input:
args <- commandArgs(trailingOnly = TRUE)
msg <- c("Please specify a tissue type to be analyzed:\n",
	 "       Choose either 'Cortex' or 'Striatum'.")
if (interactive()) { 
	tissue <- "Cortex" 
} else if (!length(args == 1)) { 
	stop(msg) 
} else if (args[1] %in% c("Cortex","Striatum")) { 
	tissue <- args[1]
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
adjm_g <- graph_from_adjacency_matrix(adjm,mode="undirected",diag=FALSE,
				      weighted=TRUE)
netw_g <- graph_from_adjacency_matrix(netw,mode="undirected",diag=FALSE,
				      weighted=TRUE)
ppin_g <- graph_from_adjacency_matrix(ppin,mode="undirected",diag=FALSE)

# We have three graphs:
graphs <- list(
	       adjm=adjm_g, # The protein co-variation matrix.
	       netw=netw_g, # The enhanced network.
	       ppin=ppin_g  # The protein-protein interaction graph.
	       )

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
createCytoscapeGraph <- function(input_graph,module=NULL) {}

module_name <- names(modules)[1]
module <- modules[[1]]
input_graph <- graphs[[1]] # adjm
nodes <- module

## NOTE: Sys.sleep()'s are important to prevent R from getting
# ahead of Cytoscape!
suppressPackageStartupMessages({
	library(RCy3)
})

# Check that we are connected to Cytoscape.
response <- cytoscapePing()
if (response != "You are connected to Cytoscape!") { stop() }

# Subset graph: keep nodes in module.
idx <- match(nodes, names(V(input_graph)))
subgraph <- induced_subgraph(input_graph, vids = V(input_graph)[idx])

# Prune weak edges.
n_edges <- length(E(subgraph))

# Seq from min weight to max to generate thresholds.
max_weight <- max(E(subgraph)$weight)
min_weight <- min(E(subgraph)$weight)
cutoffs <- seq(min_weight, max_weight, by = 0.01)

# Function to check if graph is connnected at a given edge weight threshold.
is_connected <- function(graph,threshold) {
	filt_graph <- delete.edges(graph, 
			       which(E(graph)$weight <= threshold))
	return(is.connected(graph))
}

# Check if graph is connected or not at various thresholds.
checks <- sapply(cutoffs,function(threshold) is_connected(subgraph,threshold))

# Limit is maximum cutoff level at which the graph is still connected.
limit <- cutoffs[max(which(checks==TRUE))]

# Prune edges.
# NOTE: This removes all edge types.
g <- delete.edges(subgraph, which(E(subgraph)$weight <= limit))

# Write graph to file this is faster than sending to cytoscape.
myfile <- file.path(netsdir, paste0(module_name, ".gml"))
write_graph(g, myfile, format = "gml")

# Send to Cytoscape.
## FIXME: underscores from edge weight attributes are removed!
winfile <- gsub("/mnt/d/", "D:/", myfile)
cys_net <- importNetworkFromFile(winfile)
Sys.sleep(2)
unlink(myfile)

# Create a visual style.
style.name <- paste(module_name, "style", sep = "-")
# DEFAULTS:
defaults <- list(
NODE_FILL_COLOR = col2hex("gray"),
NODE_TRANSPARENCY = 200,
NODE_SIZE = 35,
NODE_SHAPE = "ellipse",
NODE_LABEL_TRANSPARENCY = 255,
NODE_LABEL_FONT_SIZE = 12,
NODE_LABEL_COLOR = col2hex("black"),
NODE_BORDER_TRANSPARENCY = 200,
NODE_BORDER_WIDTH = 4,
NODE_BORDER_PAINT = col2hex("black"),
NODE_TRANSPARENCY = 200,
EDGE_STROKE_UNSELECTED_PAINT = col2hex("black"),
EDGE_WIDTH = 2,
NETWORK_BACKGROUND_PAINT = col2hex("white")
)
# MAPPED PROPERTIES:
mappings <- list(
NODE_LABEL = mapVisualProperty("node label", "symbol", "p"),
NODE_FILL_COLOR = mapVisualProperty("node fill color", "color", "p"),
NODE_SIZE = mapVisualProperty(
"node size",
"kme",
"c",
c(min(V(g)$kme), max(V(g)$kme)),
c(25, 75)
),
EDGE_TRANSPARENCY = mapVisualProperty(
"edge transparency",
"weight",
"c",
c(min(E(g)$weight), max(E(g)$weight)),
c(155, 255)
),
EDGE_STROKE_UNSELECTED_PAINT = mapVisualProperty(
"edge stroke unselected paint",
"weight", "c",
c(min(E(g)$weight), max(E(g)$weight)),
c(col2hex("gray"), col2hex("dark red"))
)
)
# Create a visual style.
createVisualStyle(style.name, defaults = defaults, mappings = mappings)
# Apply to graph.
setVisualStyle(style.name)
Sys.sleep(3)
# Set NS nodes to gray.
anyNS <- length(names(V(g))[which(V(g)$sigProt == 0)]) > 0
if (anyNS) {
setNodePropertyBypass(
node.names = names(V(g))[which(V(g)$sigProt == 0)],
new.values = col2hex("gray"),
visual.property = "NODE_FILL_COLOR",
bypass = TRUE,
)
setNodePropertyBypass(
node.names = names(V(g))[which(V(g)$sigProt == 0)],
new.values = 200,
visual.property = "NODE_TRANSPARENCY",
bypass = TRUE,
)
}
# Collect PPI edges.
subg <- induced_subgraph(ppi_graph, vids = V(ppi_graph)[match(nodes, names(V(ppi_graph)))])
edge_list <- apply(as_edgelist(subg, names = TRUE), 1, as.list)
# If edge list is only of length 1, unnest it to avoid problems.
if (length(edge_list) == 1) {
edge_list <- unlist(edge_list, recursive = FALSE)
}
# Add PPI edges to Cytoscape graph.
if (length(edge_list) > 0) {
ppi_edges <- addCyEdges(edge_list)
# Add PPIs and set to black.
selected_edges <- selectEdges(ppi_edges, by.col = "SUID")
# Set to black with edge bend.
setEdgePropertyBypass(
edge.names = selected_edges$edges,
new.values = col2hex("black"),
visual.property = "EDGE_STROKE_UNSELECTED_PAINT",
bypass = TRUE
)
setEdgePropertyBypass(
edge.names = selected_edges$edges,
new.values = TRUE,
visual.property = "EDGE_BEND",
bypass = TRUE
)
} # Ends IF statement.
clearSelection()
Sys.sleep(2)
# Apply layout.
layoutNetwork(network_layout)
Sys.sleep(2)
fitContent()
# Save Image..
if (!is.null(image_file)) {
# If image exists, first remove it.
if (file.exists(paste(image_file, image_format, sep = "."))) {
unlink(paste(image_file, image_format, sep = "."))
}
winfile <- gsub("/mnt/d/", "D:/", image_file)
exportImage(winfile, image_format)
}

# Free up some memory.
cytoscapeFreeMemory()

# When done, save cytoscape session.
#myfile <- file.path(netsdir,paste0(part_type,".cys"))
#winfile <- gsub("/mnt/d/","D:/",myfile)
#saveSession(winfile)
