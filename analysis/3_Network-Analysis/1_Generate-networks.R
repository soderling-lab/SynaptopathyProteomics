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

# Functions.
TBmiscr::load_all()

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Load the normalized expression data.
# Combined and normalized data, sample level outliers removed.
myfile <- file.path(datadir, "combined.rda")
load(myfile)
cleanDat <- combined

# Load sample traits.
myfile <- file.path(datadir, "samples.rda")
load(myfile)
traits <- samples

# Remove QC data, log2 transform, and coerce to data.table.
idx <- match(colnames(cleanDat), rownames(traits))
out <- traits$SampleType[idx] == "QC"
data <- as.data.table(log2(cleanDat[, !out]))
rownames(data) <- rownames(cleanDat)

# Drop any trait rows that are not in data == remove outlier samples.
traits <- as.data.table(traits) %>% filter(SampleID %in% colnames(data))

# Cortex and Striatum samples.
samples <- list(
  "Cortex" = traits$SampleID[traits$Tissue == "Cortex"],
  "Striatum" = traits$SampleID[traits$Tissue == "Striatum"]
)

# Subset data.
subDat <- lapply(samples, function(x) as.matrix(data %>% select(all_of(x))))

# Fix rownames.
subDat <- lapply(subDat, function(x) {
  rownames(x) <- rownames(data)
  return(x)
})

# Save expression data to file.
myfiles <- file.path(datadir, paste0(tolower(names(subDat)), ".rda"))
invisible(mapply(function(x, y) save(x, file=y,version=2), subDat, myfiles))

# Create signed adjacency (correlation) matrices.
adjm <- lapply(subDat, function(x) {
  silently({
    WGCNA::bicor(t(x))
  })
})

# Perform network enhancement.
message("\nPerforming network enhancement, this will take several minutes...")
adjm_ne <- lapply(adjm,neten)
names(adjm_ne) <- paste(names(adjm_ne),"NE",sep="_")

# Combine with other networks.
adjm <- c(adjm,adjm_ne)

# Coerce adjm to data.tables.
adjm <- lapply(adjm, as.data.table)

# Fix names of adjmatrices.
adjm <- lapply(adjm, function(x) {
  rownames(x) <- colnames(x) <- rownames(data)
  return(x)
})

# Write correlation matrices to .csv.
myfiles <- file.path(datadir, paste0(tolower(names(adjm)), "_adjm.csv"))
invisible(mapply(function(x, y) as.data.table(x,keep.rownames="Accession") %>% fwrite(y, row.names = TRUE), adjm, myfiles))

# Save correlation matrices as RData.
myfiles <- file.path(rdatdir, paste0(tolower(names(adjm)), "_adjm.rda"))
invisible(mapply(function(x, y) save(x,file=y,version=2), adjm, myfiles))

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

# Write to file.
myfile <- file.path(rdatdir,"PPI_Adjm.csv")
fwrite(as.data.table(PPIadjm),myfile,row.names=TRUE)

# Save as Rdata.
myfile <- file.path(rdatdir,"PPI_Adjm.RData")
saveRDS(PPIadjm,myfile)

message("Done!\n")
