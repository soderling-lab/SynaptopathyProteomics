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
