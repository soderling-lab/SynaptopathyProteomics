#!/usr/bin/env Rscript

## Inputs:
analysis_type = "Cortex"

# Load renv.
root <- getrd()
renv::load(root)

# Imports:
library(igraph)
library(ggplot2)

# Load functions in root/R.
devtools::load_all()

# Paths to input files:
adjm_file <- c(Cortex="Cortex_Adjm.RData",Striatum="")[analysis_type]
netw_file <- c(Cortex="Cortex_NE_Adjm.RData",Striatum="")[analysis_type]

# Load the adjacency matrix and network from root/rdata.
rdatdir <- file.path(root,"rdata")
adjm <- load_data(file.path(rdatdir,adjm_file),method="RData")
netw <- load_data(file.path(rdatdir,netw_file),method="RData")

# Create a igraph graphs from the adjmatrix and network.
adjms <- lapply(list(adjm,netw),as.matrix)
names(adjms) <- c("adjm","netw")
graphs <- lapply(adjms,graph_from_adjacency_matrix,
		 weighted=TRUE,mode="directed",diag=FALSE)

# Get x-y coordinates of nodes with LGL layout.
# NOTE: This takes awhile.
# layouts: lgl, kk, fr
coords_list <- lapply(graphs,layout_with_lgl)

coords_list <- lapply(graphs,layout_with_kk)

coords_list <- lapply(graphs,layout_with_fr)

# Clean-up up coordinate matrices for ggplot.
clean_up <- function(dm){
	colnames(dm) <- c("x","y")
	rownames(dm) <- V(g)$name
	dm <- as.data.frame(dm)
	return(dm)
}

coords_list <- lapply(coords_list,clean_up)

# Plot these with ggplot.
plot <- function(coords_df) {
	p <- ggplot(coords_df,aes(x,y)) + geom_point()
	return(p)
}

# Doesn't look great...
p1 <- plot(coords_list$adjm)
p2 <- plot(coords_list$netw)
cowplot::plot_grid(p1,p2,rel_heights=1,rel_widths=c(1,1))
