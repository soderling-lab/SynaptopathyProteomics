#!/usr/bin/env Rscript

#' ---
#' title:
#' description: generate networks
#' authors: Tyler W. Bradshaw
#' ---

#------------------------------------------------------------------------------
## Prepare the workspace.
#------------------------------------------------------------------------------

## Parse command line input:
# Analysis (tissue) type: cortex (1) or striatum(2).
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
	message(paste0("Analyzing ", tissue, "..."))
}

# Load renv.
root <- getrd()
renv::load(root)

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(WGCNA)
  library(neten)
  library(getPPIs)
  library(data.table)
})

# Load additional functions in root/R.
suppressWarnings({ devtools::load_all() })

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#------------------------------------------------------------------------------
## Generate protein correlation matrix.
#------------------------------------------------------------------------------

# Load the normalized expression data.
# Combined and normalized data, sample level outliers removed.
data(combined_protein)

# Load sample traits.
data(samples)

# Remove QC data, log2 transform, and coerce to data.table.
idx <- match(colnames(combined_protein), rownames(samples))
out <- samples$SampleType[idx] == "QC"
data_filt <- as.data.table(log2(combined_protein[, !out]))
rownames(data_filt) <- rownames(combined_protein)

# Drop any trait rows that are not in data == remove outlier samples.
samples <- as.data.table(samples) %>% filter(SampleID %in% colnames(data_filt))

# Subset data.
sub_samples <- samples$SampleID[samples$Tissue == tissue]
data <- as.matrix(data_filt %>% select(all_of(sub_samples)))

# Fix rownames.
rownames(data) <- rownames(data_filt)

# Create signed adjacency (correlation) matrices.
adjm <- silently({ WGCNA::bicor(t(data)) })

# Perform network enhancement.
message("\nPerforming network enhancement, this will take several minutes...")
netw <- neten(adjm)

#---------------------------------------------------------------------
## Generate PPI graph.
#---------------------------------------------------------------------

# Load protein identifier map.
prot_map <- readRDS(file.path(rdatdir,"Protein_ID_Map.RData"))

# Load mouse interactome from getPPIs.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
orgs <- c(10090, 9606, 10116)
idx <- musInteractome$Interactor_A_Taxonomy %in% orgs
ppis <- subset(musInteractome, idx)

# Save PPIs data frame--this contains PPI evidence information.
#myfile <- file.path(rdatdir,"PPI_Data.csv")
#fwrite(ppis,myfile)

# Get entrez IDs for all proteins in the network.
entrez <- prot_map$entrez[match(rownames(data),prot_map$ids)]

# Build a ppi graph.
message("\nBuilding PPI graph...")
g <- buildNetwork(ppis, entrez, taxid = 10090)

# Get ppi adjacency matrix.
PPI_adjm <- as.matrix(as_adjacency_matrix(g))

# Map entrez back to protein ids.
ids <- prot_map$ids[match(colnames(PPI_adjm),prot_map$entrez)]
colnames(PPI_adjm) <- rownames(PPI_adjm) <- ids

# One protein is mapped to the same entrez id.
missing <- prot_map$id[prot_map$ids %notin% colnames(PPI_adjm)]
names(missing) <- prot_map$entrez[match(missing,prot_map$id)]
duplicated <- prot_map$id[which(prot_map$entrez ==  names(missing))]

# Add missing column and row. 
# Just copy column and row for other protein.
missing_col <- PPI_adjm[,duplicated[duplicated %notin% missing]]
PPI_adjm <- cbind(PPI_adjm,missing_col)
colnames(PPI_adjm)[ncol(PPI_adjm)] <- missing
missing_row <- PPI_adjm[duplicated[duplicated %notin% missing],]
PPI_adjm <- rbind(PPI_adjm,missing_row)
rownames(PPI_adjm)[nrow(PPI_adjm)] <- missing

# Evaluate scale free fit.
dc <- apply(PPI_adjm,2,sum) # Node degree.
fit <- WGCNA::scaleFreeFitIndex(dc,nBreaks=10,removeFirst=FALSE)
r <- fit$Rsquared.SFT
message(paste("\nScale free fit of PPI graph:",round(r,3)))

#--------------------------------------------------------------------
## Save everything to file.
#--------------------------------------------------------------------

# Cortex and Striatum data will be saved as an R object in root/data.

# Save list containing cortex data, adjm, and network:
data_list <- list(Data = data,
		  Adjm = adjm,
		  Netw = netw,
		  Meta = samples %>% filter(Tissue == tissue),
		  Description = c("Final normalized protein data.",
		  		  "Bicor correlation matrix.",
				  "Enhanced correlation matrix.",
				  "Sample meta data."))
myfile <- file.path(datadir, paste0(tissue,"_data.rda"))
object <- tolower(paste0(tissue,"_data"))
eval(parse(text=paste(object,"= data_list")))
eval(parse(text=paste0("save(",object, ", file = myfile, version = 2)")))

# Save correlation matrices as csv.
myfile <- file.path(rdatdir,paste0(tissue,"_Adjm.csv"))
adjm <- as.data.table(adjm,keep.rownames="Accession")
fwrite(adjm,myfile)

# Write enhanced coorelation matrices (networks) as csv.
myfile <- file.path(rdatdir,paste0(tissue,"_NE_Adjm.csv"))
netw <- as.data.table(netw,keep.rownames="Accession")
fwrite(netw,myfile)

# Done!
end <- Sys.time()
message(paste("\nCompleted analysis at:", end))
message(paste("Elapsed time:",
	      round(difftime(end,start,units="secs"),2),"seconds."))
