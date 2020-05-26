#!/usr/bin/env Rscript

#' ---
#' title:
#' description: generate networks
#' authors: Tyler W. Bradshaw
#' ---

#------------------------------------------------------------------------------
## Generate protein correlation matrix.
#------------------------------------------------------------------------------

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
TBmiscr::load_all()

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Load the normalized expression data.
# Combined and normalized data, sample level outliers removed.
myfile <- file.path(datadir, "combined_protein.rda")
load(myfile)

# Load sample traits.
myfile <- file.path(datadir, "samples.rda")
load(myfile)
traits <- samples

# Remove QC data, log2 transform, and coerce to data.table.
idx <- match(colnames(combined_protein), rownames(traits))
out <- traits$SampleType[idx] == "QC"
data <- as.data.table(log2(combined_protein[, !out]))
rownames(data) <- rownames(combined_protein)

# Drop any trait rows that are not in data == remove outlier samples.
traits <- as.data.table(traits) %>% filter(SampleID %in% colnames(data))

# Cortex and Striatum samples.
samples <- list(
  "Cortex" = traits$SampleID[traits$Tissue == "Cortex"],
  "Striatum" = traits$SampleID[traits$Tissue == "Striatum"]
)

# Subset data.
data_list <- lapply(samples, function(x) as.matrix(data %>% select(all_of(x))))

# Fix rownames.
data_list <- lapply(data_list, function(x) {
  rownames(x) <- rownames(data)
  return(x)
})


# Create signed adjacency (correlation) matrices.
adjm_list <- lapply(data_list, function(x) {
  silently({
    WGCNA::bicor(t(x))
  })
})

# Perform network enhancement.
message("\nPerforming network enhancement, this will take several minutes...")
netw_list <- lapply(adjm_list, neten)

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

# Save PPIs data frame--this contains evidence information.
myfile <- file.path(rdatdir,"PPI_Data.csv")
fwrite(ppis,myfile)

# Get entrez IDs for all proteins in co-expression networks.
entrez <- prot_map$entrez

# Build a ppi graph.
message("\nBuilding PPI graph...")
g <- buildNetwork(ppis, entrez, taxid = 10090)

# Get ppi adjacency matrix.
PPIadjm <- as.matrix(as_adjacency_matrix(g))

# Map entrez back to protein ids.
ids <- prot_map$ids[match(colnames(PPIadjm),prot_map$entrez)]
colnames(PPIadjm) <- rownames(PPIadjm) <- ids

# One protein is mapped to the same entrez id.
missing <- prot_map$id[prot_map$ids %notin% colnames(PPIadjm)]
names(missing) <- prot_map$entrez[match(missing,prot_map$id)]
duplicated <- prot_map$id[which(prot_map$entrez ==  names(missing))]

# Add missing column and row. 
# Just copy column and row for other protein.
missing_col <- PPIadjm[,duplicated[duplicated %notin% missing]]
PPIadjm <- cbind(PPIadjm,missing_col)
colnames(PPIadjm)[ncol(PPIadjm)] <- missing
missing_row <- PPIadjm[duplicated[duplicated %notin% missing],]
PPIadjm <- rbind(PPIadjm,missing_row)
rownames(PPIadjm)[nrow(PPIadjm)] <- missing

# Evaluate scale free fit.
dc <- apply(PPIadjm,2,sum) # Node degree.
fit <- WGCNA::scaleFreeFitIndex(dc,nBreaks=10,removeFirst=FALSE)
r <- fit$Rsquared.SFT
message(paste("\nScale free fit of PPI graph:",round(r,3)))

#--------------------------------------------------------------------
## Save everything to file.
#--------------------------------------------------------------------

# Cortex and Striatum data will be saved as an R object in root/data.

# Save list containing cortex data, adjm, and network:
cortex_data <- list(Data = data_list$Cortex,
		    Adjm = adjm_list$Cortex,
		    Netw = netw_list$Cortex,
		    Description = c("Final normalized protein data.",
				    "Bicor correlation matrix.",
				    "Enhanced correlation matrix."))
myfile <- file.path(datadir, "cortex_data.rda")
save(cortex_data, file = myfile, version = 2)

# Save list containing striatum data, adjm, and network:
striatum_data <- list(Data = data_list$Striatum,
		      Adjm = adjm_list$Striatum,
		      Netw = netw_list$Striatum,
		      Description = c("Final normalized protein data.",
		  		      "Bicor correlation matrix.",
				      "Enhanced correlation matrix."))
myfile <- file.path(datadir, "striatum_data.rda")
save(striatum_data, file = myfile, version = 2)

# Save correlation matrices as csv.
myfile <- file.path(rdatdir,"Cortex_Adjm.csv")
adjm <- as.data.table(adjm_list$Cortex,keep.rownames="Accession")
fwrite(adjm,myfile)

myfile <- file.path(rdatdir,"Striatum_Adjm.csv")
adjm <- as.data.table(adjm_list$Striatum,keep.rownames="Accession")
fwrite(adjm,myfile)

# Write enhanced coorelation matrices (networks) as csv.
myfile <- file.path(rdatdir,"Cortex_NE_Adjm.csv")
netw <- as.data.table(netw_list$Cortex,keep.rownames="Accession")
fwrite(netw,myfile)

myfile <- file.path(rdatdir,"Striatum_NE_Adjm.csv")
netw <- as.data.table(adjm_list$Striatum,keep.rownames="Accession")
fwrite(netw,myfile)

message("Done!\n")
