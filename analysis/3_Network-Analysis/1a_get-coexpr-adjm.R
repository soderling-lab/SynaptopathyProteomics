#!/usr/bin/env Rscript 

#' ---
#' title: 1a_calc-adjacency.R
#' description: generate co-expresion adjacency matrices.
#' authors: Tyler W. Bradshaw
#' ---

#------------------------------------------------------------------------------
## Generate protein correlation matrix.
#------------------------------------------------------------------------------

# Load renv.
renv::load(getrd())

# Imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(WGCNA)
  library(neten)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatadir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
funcdir <- file.path(root, "R")

# Additional functions.
TBmiscr::load_all()

# Load the normalized expression data.
# Combined and normalized data, sample level outliers removed.
myfile <- file.path(rdatadir, "2_Combined_cleanDat.RData")
cleanDat <- readRDS(myfile)

# Load sample traits.
myfile <- file.path(rdatadir, "2_Combined_traits.RData")
traits <- readRDS(myfile)
traits$Model <- gsub(" ","",traits$Model) # remove spaces

# Remove QC data, log2 transform, and coerce to data.table.
out <- traits$SampleType[match(colnames(cleanDat), rownames(traits))] == "QC"
data <- as.data.table(log2(cleanDat[, !out]))
rownames(data) <- rownames(cleanDat)

# Drop any trait rows that are not in data == remove outlier samples.
traits <- as.data.table(traits) %>% filter(SampleID %in% colnames(data))

# WT and KO samples.
samples <- list(
  "WT" = traits$SampleID[traits$SampleType == "WT"],
  "KO" = traits$SampleID[traits$SampleType == "KO" | traits$SampleType == "HET"],
  "Cortex" = traits$SampleID[traits$Tissue == "Cortex"],
  "Striatum" = traits$SampleID[traits$Tissue == "Striatum"]
)

# Subset WT and KO data.
subDat <- lapply(samples, function(x) as.matrix(data %>% select(all_of(x))))

# Fix rownames.
subDat <- lapply(subDat, function(x) {
  rownames(x) <- rownames(data)
  return(x)
})

# Add combined data.
subDat[["Combined"]] <- data

# Save expression data to file.
myfiles <- file.path(rdatadir, paste0("3_", names(subDat), "_cleanDat.RData"))
invisible(mapply(function(x, y) saveRDS(x, y), subDat, myfiles))

# Create signed adjacency (correlation) matrices.
adjm <- lapply(subDat, function(x) {
  silently({
    WGCNA::bicor(t(x))
  })
})

# Perform network enhancement.
message("Performing network enhancement, this will take several minutes...")
adjm_ne <- lapply(adjm, neten)
names(adjm_ne) <- paste(names(adjm_ne), "NE", sep = "_")

# Combine with other networks.
adjm <- c(adjm, adjm_ne)

# Coerce adjm to data.tables.
adjm <- lapply(adjm, as.data.table)

# Fix names of adjmatrices.
adjm <- lapply(adjm, function(x) {
  rownames(x) <- colnames(x) <- rownames(data)
  return(x)
})

# Write correlation matrices to .csv.
myfiles <- file.path(rdatadir, paste0("3_", names(adjm), "_Adjm.csv"))
invisible(mapply(function(x, y) fwrite(x, y, row.names = TRUE), adjm, myfiles))

# Save correlation matrices RData.
myfiles <- file.path(rdatadir, paste0("3_", names(adjm), "_Adjm.RData"))
invisible(mapply(function(x, y) saveRDS(x, y), adjm, myfiles))

quit()

#--------------------------------------------------------------------
## Identify subset of highly reproducible proteins
#--------------------------------------------------------------------

# Tidy up the data.
dt <- reshape2::melt(cleanDat,value.name="Abundance")
colnames(dt)[c(1,2)] <- c("Accession", "Sample")
dt$Sample <- as.character(dt$Sample)

# Annotate with additional meta data.
traits$Sample <- as.character(traits$SampleID)
traits$Channel <- sapply(strsplit(traits$SampleID,"\\."),"[",2)
data <- left_join(dt,traits, by=c("Sample"))
data <- data %>% select(Accession,Tissue,Model,SampleType, Channel, Abundance)

# Remove QC data.
data <- data %>% filter(SampleType != "QC")

# Split the data into a list--each element is the data for a protein.
data_list <- data %>% group_by(Accession) %>% group_split()
names(data_list) <- sapply(data_list,function(x) unique(x$Accession)) # Prots

# Define a function that checks reproducibility of a protein.
check_reproducibility <- function(data_list,protein,tissue.type,
				  fun="bicor",threshold=0.6) {
	# Which proteins are highly reproducible between replicates.
	# For a given protein, cast the data into a matrix: Fraction ~ Replicate.
	df <- data_list[[protein]] %>% 
		filter(SampleType != "WT",Tissue==tissue.type) 
	df <- as.data.table(df) 
	# Treat all mutants the same:
	df$SampleType <- gsub("HET","KO",df$SampleType) 
	df$Replicate <- as.numeric(as.factor(df$Channel))
        temp <- dcast(df, Model + Tissue ~ SampleType + Replicate, 
		      value.var="Abundance") 
	cols <- grep("_",colnames(temp))
	dm <- temp %>% select(all_of(cols)) %>% 
		as.matrix(rownames.value=paste(temp$Model,temp$Tissue,sep="."))
	dm <- log2(dm)
	# missing values can arise in dm if outlier sample was removed.
	# Calculate the pairwise correlation between samples.
	opts <- "pairwise.complete.obs"
	if (fun == "bicor") { cormat <- bicor(dm,use=opts) }
	if (fun == "pearson") { cormat <- cor(dm,method="pearson",use=opts) }
	if (fun == "spearman") { cormat <- cor(dm,method="spearman",use=opts) }
	# Ignore self- and duplicate- comparisons.
	diag(cormat) <- NA
	cormat[lower.tri(cormat)] <- NA
	# Melt into a vector of correlation values.
	x <- reshape2::melt(cormat,na.rm=TRUE,value.name="cor")[["cor"]]
	# Check if all values are > threshold.
	check <- all(x > threshold)
	return(check)
}

# Check protein reproducibility.
check <- sapply(names(data_list),
		function(prot) {
			check_reproducibility(data_list,prot,
					      tissue.type="Cortex",
					      fun="bicor",threshold=0.6)
		})
message(paste("Numer of reproducible proteins:",sum(check)))

check <- sapply(names(data_list),
		function(prot) {
			check_reproducibility(data_list,prot,
					      tissue.type="Striatum",
					      fun="bicor",threshold=0.6)
		})
message(paste("Numer of reproducible proteins:",sum(check)))

ids <- names(which(check))


# Status.
n_markers <- formatC(sum(check),big.mark=",")
message(paste("Identified",n_markers,"potential protein markers."))

# List of potential markers:
potential_markers <- list("id" = names(which(ids)),
			  "entrez" = names(which(entrez)))

