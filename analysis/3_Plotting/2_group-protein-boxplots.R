#!/usr/bin/env Rscript

## Inputs:
#analysis_type = "Cortex"
analysis_type = "Striatum"

# Load renv.
root <- getrd()
renv::load(root)

# Imports:
suppressPackageStartupMessages({
	library(igraph)
	library(ggplot2)
})

# Load functions in root/R and data in root/data.
suppressWarnings({ devtools::load_all() })

# Directory for saving figures:
figsdir <- file.path(root,"figs")

# Load protein data and graph partitions:
expr_data <- paste0(tolower(analysis_type),"_data")
part_data <- paste0(tolower(analysis_type),"_partition")
data(list=c(expr_data,part_data))

# Get the data we are gonna analyze:
if (analysis_type == "Cortex") {
	data <- cortex_data$Data
	partition <- cortex_partition
} else {
	data <- striatum_data$Data
	partition <- striatum_partition
}

# NOTE: Apply is done row-wise (dim = 1), but the result is transposed.
# See://stackoverflow.com/questions/9521260/
dm_norm <- t(apply(2^data, 1, function(x) x/sum(x)))
plot <- ggplotPCAprot(dm_norm,scale=TRUE,center=TRUE)

# Save the plot.
myfile <- file.path(figsdir,paste0(analysis_type,"_Protein_PCA.pdf"))
ggsave(file=myfile,plot,width=6,height=6)
