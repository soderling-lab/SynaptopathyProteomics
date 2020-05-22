#!/usr/bin/env Rscript

#' ---
#' title: 1a_calc-adjacency.R
#' description: generate co-expresion adjacency matrices.
#' authors: Tyler W. Bradshaw
#' ---

#------------------------------------------------------------------------------
## Generate protein correlation matrix.
#------------------------------------------------------------------------------

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

# Functions.
myfun <- list.files(funcdir,full.names = TRUE)
invisible(sapply(myfun, source))

# Load the normalized expression data.
# Combined and normalized data, sample level outliers removed.
myfile <- file.path(rdatadir, "2_Combined_cleanDat.RData")
cleanDat <- readRDS(myfile)

# Load sample traits.
myfile <- file.path(rdatadir, "2_Combined_traits.RData")
traits <- readRDS(myfile)

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
myfiles <- file.path(rdatadir, paste0("3_", names(adjm), "_Adjm.csv"))
invisible(mapply(function(x, y) fwrite(x, y, row.names = TRUE), adjm, myfiles))

# Save correlation matrices RData.
myfiles <- file.path(rdatadir, paste0("3_", names(adjm), "_Adjm.RData"))
invisible(mapply(function(x, y) saveRDS(x, y), adjm, myfiles))
#!/usr/bin/env Rscript

## User parameters.
calc_GO_Sem_Sim = FALSE

# Imports.
suppressPackageStartupMessages({
	library(GOSemSim) # Other packages for GOSemSim:
	library(org.Mm.eg.db)
	library(AnnotationDbi)
	library(stats4) 
	library(BiocGenerics) 
	library(parallel)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatdir <- file.path(root,"rdata")
funcdir <- file.path(root,"R")

# Functions.
myfun <- list.files(funcdir,full.names=TRUE)
invisible(sapply(myfun,source))

# Load Cortex partitions.
myfile <- list.files(rdatdir,patter="10773682",full.names=TRUE)
cox_parts <- readRDS(myfile)

# Load protein identifier map.
myfile <- file.path(rdatdir,"2_Protein_ID_Map.RData")
prot_map <- readRDS(myfile)

# Get Entrez IDs.
entrez <- prot_map$entrez[match(names(cox_parts[[1]]),prot_map$ids)]

# Calculate semantic similarity between all pairwise comparisons of genes.
if (calc_GO_Sem_Sim) {
	message(paste("Calculating GO semantic similarity, 
		      this will take several hoursle."))
	# Build list of GO ontologies.
	msGO <- lapply(c("BP","CC","MF"), function(x) {
			       godata('org.Mm.eg.db', ont=x) })
	names(msGO) <- c("BP","CC","MF")
	# Calculate GOSemSim for all ontologies.
	gosemsim <- lapply(msGO,function(x) { 
				   mgeneSim(genes=entrez, 
					    semData=x, 
					    measure="Wang", # Graph-based
					    verbose=TRUE) })
	# Save to file.
	myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_Adjms.RData")
	saveRDS(gosemsim,myfile)
} else {
	# Load from file.
	message("Loading GO semantic similarity data from file!")
	myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_Adjms.RData")
	gosemsim <- readRDS(myfile)
}

# Create an empty adjm. We need to insure that GO semantic similarity matrices
# and the protein co-expression matrix are the same dimensions.
adjm <- matrix(nrow=length(entrez),ncol=length(entrez))
colnames(adjm) <- rownames(adjm) <- entrez

# Loop to enforce consistent dimensions...
output <- list()
for (i in seq_along(gosemsim)){
	dm <- gosemsim[[i]]
	missing <- colnames(adjm)[colnames(adjm) %notin% colnames(dm)]
	missing_rows <- matrix(NA,nrow=length(missing),ncol=ncol(dm))
	rownames(missing_rows) <- missing
	dm1 <- rbind(dm,missing_rows)
	missing_cols <- matrix(NA,nrow=nrow(dm1),ncol=length(missing))
	colnames(missing_cols) <- missing
	dm2 <- cbind(dm1,missing_cols)
	idx <- idy <- match(colnames(adjm),colnames(dm2))
	dm3 <- dm2[idx,idy]
	# Check.
	check <- all(colnames(dm3) == colnames(adjm)) & 
		all(rownames(dm3) == rownames(adjm))
	if (!check) { print("Problem, matrix indices don't match!") }
	output[[i]] <- dm3
}
names(output) <- names(gosemsim)

# Combined BP + MF + CC results by calculating RMS.
# This is fast, but any missing value (NA) will cause combined, RMS 
# value to be NA.
#x1 <- output[[1]]
#x2 <- output[[2]]
#x3 <- output[[3]]
#rms <- sqrt((x1^2 + x2^2 + x3^2)/3) 
#n_missing <- sum(is.na(rms))
#message(paste("Percent missing:",round(100*(n_missing/length(rms)),2)))
# Approx. 25% missing values!

# Lets try another way...
# Melt and merge GO semantic similarity data into single df.
# Use cbind, because its fast.
melted_gosemsim <- lapply(output,reshape2::melt) 
gosemsim_df <- do.call(cbind,melted_gosemsim)
gosemsim_df <- gosemsim_df[,-c(grep("Var*",colnames(gosemsim_df))[-c(1,2)])]
colnames(gosemsim_df) <- c("EntrezA","EntrezB","BP.GOsim","CC.GOsim","MF.GOsim")
gosemsim_df <- tibble::add_column(gosemsim_df,
				  Index=c(1:nrow(gosemsim_df)),.after=2)

# Loop with progress bar to calculate RMS.
# Ignore NA. Slow bc the data is huge!
# Could use apply, but pbar is helpful.
message(paste("Calculating combined GO semantic similarity score..."))
pbar <- txtProgressBar(min=1,max=nrow(gosemsim_df),style=3)
rms <- vector(mode="numeric",length=nrow(gosemsim_df))
for (i in 1:nrow(gosemsim_df)){
	setTxtProgressBar(pbar,i)
	rms[i] <- sqrt(mean(as.numeric(gosemsim_df[i,c(4:6)])^2,na.rm=TRUE))
	if (i==nrow(gosemsim_df)) { close(pbar); message("\n") }
}

# Cast into dm. 
# Should be in correct order as we are just casting the previously
# melted matrix.
gosemsim_df$GOsim.RMS <- rms
rms_adjm <- matrix(rms,nrow=nrow(adjm),ncol=ncol(adjm))
colnames(rms_adjm) <- rownames(rms_adjm) <- colnames(adjm)

# Replace missing values with 0. 
idx <- is.na(rms_adjm)
rms_adjm[idx] <- 0

# Map names back to protein ids.
ids <- prot_map$ids[match(colnames(rms_adjm),prot_map$entrez)]
colnames(rms_adjm) <- rownames(rms_adjm) <- ids

# Write to file for LA clusting!
myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_RMS_Adjm.csv")
data.table::fwrite(rms_adjm,myfile,row.names=TRUE)
#!/usr/bin/env Rscript

#' ---
#' title: 1c_get-ppi-adjm.R
#' description: Generate PPI adjacency matrix.
#' authors: Tyler W Bradshaw
#' ---

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
	library(getPPIs)
	library(ggplot2)
	library(data.table)
})

# Directories.
here <- getwd()
subdir <- dirname(here)
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables",subdir)

# Functions.
suppressWarnings({ devtools::load_all() })

# Load protein identifier map.
prot_map <- readRDS(file.path(rdatdir,"2_Protein_ID_Map.RData"))

#---------------------------------------------------------------------
## Generate PPI graph.
#---------------------------------------------------------------------

# Load mouse interactome.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
orgs <- c(10090, 9606, 10116)
idx <- musInteractome$Interactor_A_Taxonomy %in% orgs
ppis <- subset(musInteractome, idx)

# Save PPIs data frame--this contains evidence information.
myfile <- file.path(rdatdir,"3_All_PPIs.RData")
saveRDS(ppis,myfile)

# Get entrez IDs for all proteins in co-expression networks.
entrez <- prot_map$entrez

# Build a ppi graph.
message("Building PPI graph...")
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
message(paste("Scale free fit of PPI graph:",round(r,3)))

# Write to file.
myfile <- file.path(rdatdir,"3_PPI_Adjm.csv")
fwrite(as.data.table(PPIadjm),myfile,row.names=TRUE)

# Save as Rdata.
myfile <- file.path(rdatdir,"3_PPI_Adjm.RData")
saveRDS(PPIadjm,myfile)

message("Done!\n")
