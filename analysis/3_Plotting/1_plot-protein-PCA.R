#!/usr/bin/env Rscript

## Inputs:
tissue = "Cortex"

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Imports:
suppressPackageStartupMessages({
	library(igraph)
	library(ggplot2)
})

# Load functions in root/R and data in root/data.
suppressWarnings({ devtools::load_all() })

# Directory for saving figures:
figsdir <- file.path(root,"figs")

# Dynamically load protein data and graph partition.
expr_data <- paste0(tolower(tissue),"_data")
part_data <- paste0(tolower(tissue),"_partition")
data(list=c(expr_data,part_data))
eval(parse(text=paste0("data=",expr_data)))
eval(parse(text=paste0("partition=",part_data)))

# NOTE: Apply is done row-wise (dim = 1), but the result is transposed.
# See://stackoverflow.com/questions/9521260/
dm_norm <- t(apply(2^data, 1, function(x) x/sum(x)))
plot <- ggplotPCAprot(dm_norm,scale=TRUE,center=TRUE)

# Save the plot.
myfile <- file.path(figsdir,paste0(tissue,"_Protein_PCA.pdf"))
ggsave(file=myfile,plot,width=6,height=6)

# Add colors.

