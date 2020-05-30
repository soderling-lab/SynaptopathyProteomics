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
renv::load()

# Global imports.
suppressPackageStartupMessages({
  library(RCy3)
  library(dplyr)
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
dataset <- tolower(paste0(analysis_type,"_data"))
data <- list(dataset)
eval(parse(text=paste0("data_list=",dataset)))

# Load the graph partition in root/data.
dataset <- tolower(paste0(analysis_type,"_partition"))
data(list=dataset)
eval(parse(text=paste0("partition=",dataset)))

# Load gene map.
data(gene_map)

#---------------------------------------------------------------------
## Generate PPI graph.
#---------------------------------------------------------------------

## NOTE: Coerce boolean attributes to integer to avoid warnings when
# loading into cytoscape.

# Create co-expression graph.
exp_graph <- graph_from_adjacency_matrix(adjm,mode="undirected",
				  weighted=TRUE, diag=FALSE)

# Create PPI graph.
ppi_graph <- graph_from_adjacency_matrix(adjm_ppi,mode="undirected",
				  weighted=TRUE, diag=FALSE)

# Remove NAs from PPI edges.
E(ppi_graph)$weight[which(is.na(E(ppi_graph)$weight))] <- 0

## Add attributes to igraph object.
# Add Gene symbols.
symbols <- protmap$gene[match(names(V(exp_graph)),protmap$ids)]
exp_graph <- set_vertex_attr(exp_graph,"symbol",value = symbols)

# Add sigProt vertex attribute.
anySig <- as.numeric(sigProts[[data_type]][names(V(exp_graph))])
exp_graph <- set_vertex_attr(exp_graph, "sigProt", 
			 value = anySig)

# Add any DBDprot vertex attribute.
DBDnodes <- lapply(DBDprots,function(x) names(V(exp_graph)) %in% x)
for (DBD in names(DBDnodes)){
	exp_graph <- set_vertex_attr(exp_graph, name=DBD, 
				 value = as.numeric(DBDnodes[[DBD]]))
}

# Collect PPI evidence.
myfile <- file.path(rdatdir,"3_All_PPIs.RData")
ppis <- readRDS(myfile)

# Map mouse entrez to protein ids.
ppis$ProteinA <- protmap$ids[match(ppis$osEntrezA,protmap$entrez)]
ppis$ProteinB <- protmap$ids[match(ppis$osEntrezB,protmap$entrez)]
out <- is.na(ppis$ProteinA) | is.na(ppis$ProteinB)
ppis <- ppis[!out,]

# Get the relevant columns.
ppis <- ppis %>% select(ProteinA,ProteinB,osEntrezA,osEntrezB,
			Interactor_A_Taxonomy,Interactor_B_Taxonomy,
			Source_database,Confidence_score,
			Publications,Methods)

# Save to file..
myfile <- file.path(tabsdir,paste0("3_All_PPIs.csv"))
write_excel(list(PPIs=ppis),myfile)

#---------------------------------------------------------------------
## Generate cytoscape graphs.
#---------------------------------------------------------------------

# Prompt the user to open Cytoscape if it is not open.
cytoscape_ping()

# If working with Combined data, append graphs to tissue specific 
# Cytoscape file.

if (generate_cytoscape_graphs) {

	if (data_type == "Combined") {

		cysfile <- file.path(netsdir,paste0(part_type,".cys"))

		if (file.exists(cysfile)){
			message(paste("Adding graphs to",part_type,"file!"))
			winfile <- gsub("/mnt/d/","D:/",cysfile)
			openSession(winfile)
		} else {

			message(paste("Analyze",part_type,"data first.",
			              "Combined graphs will be appended to",
			              "this file."))
	}
}

# Create graphs.

	for (i in c(1:length(modules))) {

		module_name = names(modules)[i]
		message(paste("Working on module", module_name,"..."))

		nodes = names(modules[[module_name]])
		module_kme = KME_list[[module_name]]

		network_layout = 'force-directed edgeAttribute=weight'
		image_file = file.path(dirname(figsdir),"Networks",module_name)
		image_format = "SVG"

		createCytoscapeGraph(exp_graph,ppi_graph,nodes,
			     module_kme,module_name,
			     module_colors, network_layout,
			     output_file, image_file,
			     image_format)

		# When done, save cytoscape session.
		if (i == length(modules)) {
			myfile <- file.path(netsdir,paste0(part_type,".cys"))
			winfile <- gsub("/mnt/d/","D:/",myfile)
			saveSession(winfile)
