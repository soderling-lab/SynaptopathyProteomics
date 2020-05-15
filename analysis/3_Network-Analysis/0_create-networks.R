#!/usr/bin/env Rscript 

#' ---
#' title: 
#' description: generate co-expresion adjacency matrices.
#' authors: Tyler W. Bradshaw
#' ---

## Inputs:
# Input data should be in root/rdata/:
input_data = "Combined_norm_protein.csv"
input_data = "Combined_tidy_protein.csv"

## Other parameters:

## Output for downstream analysis:
# Stored in root/rdata/
# * Cortex_Adjm.csv
# * Striatum_Adjm.csv
# * Cortex_NE_Adjm.csv
# * Striatum_NE_Adjm.csv

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------

# Load renv.
renv::load(getrd())

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(WGCNA)
  library(neten)
  library(data.table)
})

# Directories.
root <- getrd()
rdatdir <- file.path(root, "rdata")

# Additional functions.
TBmiscr::load_all()

# Load final normalized data.
myfile <- file.path(rdatdir,input_data)
data <- fread(myfile)

# Log transform data.
data$Abundance <- log2(data$Abundance)

# Split into cortex and striatum datasets.
data_list <- data %>% group_by(Tissue) %>% group_split()
names(data_list) <- sapply(data_list,function(x) unique(x$Tissue))

# Function to coerce tidy data to data matrices.
tidyProt_to_matrix <- function(x) {
	dm <- as.data.table(x) %>% 
		dcast(Accession ~ Sample, value.var="Abundance") %>%
		as.matrix(rownames="Accession") %>% log2()
	return(dm)
}

# Coerce to data matrix.
dm_list <- lapply(data_list,tidyProt_to_matrix)

# Create signed adjacency (correlation) matrices.
adjm_list <- lapply(dm_list, function(x) { WGCNA::bicor(t(x)) })

# Perform network enhancement.
message("Performing network enhancement, this will take several minutes...")
ne_list <- lapply(adjm_list, neten)

# Write correlation matrices to file.
myfiles <- file.path(rdatdir, paste0(names(adjm_list),"_Adjm.csv"))
namen <- names(myfiles) <- names(adjm_list)
for (name in namen){
	as.data.table(adjm_list[[name]],keep.rownames="Accession") %>% 
		fwrite(myfiles[name])
}

# Write enhanced networks to file.
myfiles <- file.path(rdatdir, paste0(names(adjm_list),"_NE_Adjm.csv"))
namen <- names(myfiles) <- names(adjm_list)
for (name in namen){
	as.data.table(ne_list[[name]],keep.rownames="Accession") %>% 
		fwrite(myfiles[name])
}

quit()

#--------------------------------------------------------------------
## Identify subset of highly reproducible proteins
#--------------------------------------------------------------------

# Split the data into a list of proteins.
prot_list <- data %>% group_by(Accession) %>% group_split()
names(prot_list) <- sapply(prot_list,function(x) unique(x$Accession))

# Define a function that checks reproducibility of a protein.
check_reproducibility <- function(prot_list,protein,genotype,tissue.type,
				  fun="bicor",threshold=0.6) {
	# Which proteins are highly reproducible between replicates.
	# For a given protein, cast the data into a matrix: Fraction ~ Replicate.
	dt <- prot_list[[protein]] %>% 
		filter(Treatment == genotype,Tissue==tissue.type) %>%
		as.data.table() 
	# Treat all mutants the same:
	dt$Treatment <- gsub("HET","KO",dt$Treatment) 
	dt$Replicate <- as.numeric(as.factor(dt$Channel))
        dt_temp <- dcast(dt, Genotype + Tissue ~ Treatment + Replicate, 
		      value.var="Abundance") 
	value_cols <- grep("_",colnames(dt_temp))
	dm <- dt_temp %>% select(all_of(value_cols)) %>% 
		as.matrix(rownames.value=paste(dt_temp$Genotype,
					       dt_temp$Tissue,sep="."))
	# missing values can arise in dm if outlier sample was removed.
	# Calculate the pairwise correlation between samples.
	opts <- "pairwise.complete.obs"
	if (fun == "bicor") { cormat <- bicor(dm,use=opts) }
	if (fun == "pearson") { cormat <- cor(dm,method="pearson",use=opts) }
	if (fun == "spearman") { cormat <- cor(dm,method="spearman",use=opts) }
	# Ignore self- and duplicate- comparisons.
	diag(cormat) <- NA
	cormat[cormat==0] <- NA
	cormat[lower.tri(cormat)] <- NA
	# Melt into a vector of correlation values.
	x <- reshape2::melt(cormat,na.rm=TRUE,value.name="cor")[["cor"]]
	# Check if all values are > threshold.
	check <- sum(x > threshold)
	return(check)
}

# Check protein reproducibility.
check <- list()
for (prot in names(prot_list)) {
	check[[prot]] <- check_reproducibility(prot_list, prot,
				       genotype="WT",
				       tissue.type="Cortex",
				       fun="bicor",
				       threshold=0.6)
}

checks <- unlist(check)
message(paste("Numer of reproducible proteins:",sum(checks)))

#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
## Generate PPI graph.
#---------------------------------------------------------------------

# Load protein identifier map.
prot_map <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))


# Load mouse interactome.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
orgs <- c(10090, 9606, 10116)
idx <- musInteractome$Interactor_A_Taxonomy %in% orgs
ppis <- subset(musInteractome, idx)

# Save PPIs data frame--this contains evidence information.
myfile <- file.path(rdatdir, "3_All_PPIs.RData")
saveRDS(ppis, myfile)

# Get entrez IDs for all proteins in co-expression networks.
entrez <- prot_map$entrez

# Build a ppi graph.
message("Building PPI graph...")
g <- buildNetwork(ppis, entrez, taxid = 10090)

# Get ppi adjacency matrix.
PPIadjm <- as.matrix(as_adjacency_matrix(g))

# Map entrez back to protein ids.
ids <- prot_map$ids[match(colnames(PPIadjm), prot_map$entrez)]
colnames(PPIadjm) <- rownames(PPIadjm) <- ids

# One protein is mapped to the same entrez id.
missing <- prot_map$id[prot_map$ids %notin% colnames(PPIadjm)]
names(missing) <- prot_map$entrez[match(missing, prot_map$id)]
duplicated <- prot_map$id[which(prot_map$entrez == names(missing))]

# Add missing column and row.
# Just copy column and row for other protein.
missing_col <- PPIadjm[, duplicated[duplicated %notin% missing]]
PPIadjm <- cbind(PPIadjm, missing_col)
colnames(PPIadjm)[ncol(PPIadjm)] <- missing
missing_row <- PPIadjm[duplicated[duplicated %notin% missing], ]
PPIadjm <- rbind(PPIadjm, missing_row)
rownames(PPIadjm)[nrow(PPIadjm)] <- missing

# Evaluate scale free fit.
dc <- apply(PPIadjm, 2, sum) # Node degree.
fit <- WGCNA::scaleFreeFitIndex(dc, nBreaks = 10, removeFirst = FALSE)
r <- fit$Rsquared.SFT
message(paste("Scale free fit of PPI graph:", round(r, 3)))

# Write to file.
myfile <- file.path(rdatdir, "3_PPI_Adjm.csv")
fwrite(as.data.table(PPIadjm), myfile, row.names = TRUE)

# Save as Rdata.
myfile <- file.path(rdatdir, "3_PPI_Adjm.RData")
saveRDS(PPIadjm, myfile)

message("Done!\n")
