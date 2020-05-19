#!/usr/bin/env Rscript

#' ---
#' title: TAMPOR.R
#' description: TAMPOR normalization of preprocessed TMT data.
#' authors: Tyler W Bradshaw, Eric B Dammer (TAMPOR).
#' ---

## Parse command line arguments:
if (interactive()) {
	# If interactive, then define analysis_type:
	analysis_type = "Striatum" # Tissue type for analysis.
} else if (!interactive()) {
	## If not interactive, check that only 1 arg is passed.
	args <- commandArgs(trailingOnly=TRUE)
	if (length(args) == 1) { 
		analysis_type = commandArgs(trailingOnly=TRUE)[1]
	} else { 
		stop("Specify either 'Cortex' or 'Striatum'.",call.=FALSE) 
	}
}

# Input data should be in root/rdata:
input_data = list("Cortex" = "Cortex_filt_protein.csv",
		  "Striatum"= "Striatum_filt_protein.csv")

## Other parameters:
output_name = analysis_type
alpha_threshold = 0.1 # FDR threshold for signficance.
n_threads = parallel::detectCores() - 1 # Numb of cores for parallel processing.

## Main Outputs:
# Stored in root/tables/

## Output for downstream analysis:
# Stored in root/rdata/

## Order of data processing operations:

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. 
# Load custom functions and prepare the project directory for saving 
# output files.

# Load renv.
root <- getrd()
renv::load(root)

# Load required packages.
suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(tibble)
  library(data.table)
})

# Load additional functions:
TBmiscr::load_all()

# Set any other directories.
rdatdir <- file.path(root, "rdata")

#---------------------------------------------------------------------
## Reformat the data for TAMPOR
#---------------------------------------------------------------------

# Collect cortex and striatum data in a list and then bind together.
data_list <- lapply(input_data,function(x) fread(file.path(rdatdir,x)))
combined_dt <- bind_rows(data_list) %>% as.data.table()
data_list[["Combined"]] <- combined_dt

# Which subset are we analyzing?
message(paste0("\nAnalyzing ",analysis_type,"..."))
combined_dt <- data_list[[analysis_type]] %>% 
	filter(Treatment != "QC") %>% as.data.table()

# Add batch.channel annotation.
batch <- interaction(combined_dt$Tissue,combined_dt$Genotype) %>% as.numeric()
batch <- paste0("b",batch)
combined_dt$batch <- batch
combined_dt$ID <- interaction(batch,combined_dt$Channel) %>% as.character()

# Cast into a data.matrix.
dm <- combined_dt %>% dcast(Accession ~ ID, value.var="Intensity") %>% 
	as.matrix(rownames="Accession")

# Remove rows with any missing values.
missing_vals <- apply(dm,1,function(x) any(is.na(x)))
if (sum(missing_vals) > 0) {
	message(paste("\nRemoving", sum(missing_vals),
		      "rows that contain missing values."))
	dm <- dm[!missing_vals,]
}

# Combined data for TAMPOR:
combined_dm <- dm

# Traits must include batch and ID columns.
traits <- combined_dt %>% select(batch,ID) %>% unique()
traits <- as.data.frame(traits)
rownames(traits) <- traits[,"ID"]

# Global internal standards (GIS) for normalization is all WT samples.
# WT Cortex and WT Striatum samples will be scaled by TAMPOR to be 
controls <- unique(filter(combined_dt,Treatment=="WT")[["ID"]])

#---------------------------------------------------------------------
## TAMPOR Normalization.
#---------------------------------------------------------------------
# Perform TAMPOR normalization.

data_tampor <- TAMPOR(
  dat = combined_dm,
  traits = traits,
  batchPrefixInSampleNames = TRUE,
  samplesToIgnore = "None",
  GISchannels = controls,
  parallelThreads = n_threads
)

# Collect normalized, relative abundance data post-TAMPOR.
tampor_protein <- data_tampor$cleanRelAbun
dt <- reshape2::melt(tampor_protein,value.name="Abundance")
colnames(dt)[c(1,2)] <- c("Accession","ID")
dt$ID <- as.character(dt$ID)
dt$Accession <- as.character(dt$Accession)
norm_protein <- left_join(combined_dt,dt,by=c("Accession","ID"))

#---------------------------------------------------------------------
## Identify Sample outliers.
#---------------------------------------------------------------------

# Calculate Oldham's normalized sample connectivity (zK) in order to 
# identify outlier samples.
# This approach was adapted from Oldham et al., 2012 (pmid: 22691535).
zK <- sampleConnectivity(norm_protein %>% filter(Treatment != "QC"))

# Outliers:
sample_connectivity_threshold = 2.5
outlier_samples <- c(names(zK)[zK < -sample_connectivity_threshold],
		     names(zK)[zK > sample_connectivity_threshold])

# Remove sample outliers.
norm_protein <- norm_protein %>% filter(Sample %notin% outlier_samples)

# Status:
if (length(outlier_samples) == 0) {
	message(paste("\nFinal number of samples:",
		      length(unique(norm_protein$Sample))))
} else {
	message(paste0("\nRemoved outlier sample(s):\n",
	       paste(outlier_samples,collapse="\n")))
	message(paste("Final number of samples:",
		      length(unique(norm_protein$Sample))))
}

#--------------------------------------------------------------------
## Identify subset of highly reproducible proteins
#--------------------------------------------------------------------
# NOTE: scaling WT's to be equal scrubs WT variability and consequently,
# there are few to zero 'highly reproducible' proteins.

# A 'highly reproducible' protein is a protein whose expression profile
# is highly reproducible-- that is, any given WT sample is highly coorelated
# with all other WT replicates.

# Check reproducibility of WT protein expression after IRS normalization,
# This can also be done after fitting the glm, but the number of 
# reproducible proteins may be inflated.
message(paste("\nChecking reproducibility of WT protein expression..."))

# Split the data into a list of proteins.
prot_list <- norm_protein %>% group_by(Accession) %>% group_split()
names(prot_list) <- sapply(prot_list,function(x) unique(x$Accession))

# Check protein reproducibility.
check <- list()
for (prot in names(prot_list)) {
	check[[prot]] <- check_reproducibility(prot_list, prot,
				       treatment.subset="WT",
				       fun="pearson",
				       threshold=0.8)
}

# Status.
checks <- unlist(check)
n <- max(checks)
reproducible_prots <- names(which(checks==n))
message(paste("Number of highly reproducible proteins (potential markers):",
	      length(reproducible_prots)))

#---------------------------------------------------------------------
## Evaluate protein differential abundance.
#---------------------------------------------------------------------
# Use EdgeR glm to evaluate differential protein abundance.
# Compare WT v Mutant.

# Status.
message(paste("\nAnalyzing protein differential abundance",
	      "with EdgeR GLM..."))

# Asses changes in protein abundance using a glm to account for
# differences in genetic background (genotype).
# NOTE: This script combines WTs!

## Pick a model
# model = "~0 + group"
comparisons = c("Genotype","Treatment")
model = "~background + group"
glm_results <- glmDA(norm_protein, comparisons, model, value.var="Abundance")

# Extract contrasts of interest.
idx <- c((length(glm_results)-3):length(glm_results))
glm_results <- glm_results[idx]
names(glm_results) <-  c("Shank3","Syngap1","Ube3a","Shank2")

#  Annotate with gene IDs.
f <- function(x) {
	Symbol <- norm_protein$Symbol[match(x$Accession,norm_protein$Accession)]
	Entrez <- norm_protein$Entrez[match(x$Accession,norm_protein$Accession)]
	x <-  tibble::add_column(x,Symbol,.after="Accession")
	x <- tibble::add_column(x,Entrez,.after="Symbol")
	return(x)
}
glm_results <- lapply(glm_results,f)

# Summary of DA proteins.
message(paste0("Summary of differentially abundant proteins at FDR <",
	      alpha_threshold,":"))
n <- sapply(glm_results,function(x) sum(x$FDR < alpha_threshold))
knitr::kable(as.data.table(n,keep.rownames="Genotype"))

# Check top gene in each condition.
message("Top differentially abundant proteins:")
knitr::kable(sapply(glm_results,function(x) head(x$Symbol,3)))

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:

## Save key results.
#message("\nSaving data for downstream analysis.")

# [output_name]_tidy_peptide.csv -- raw peptide data.
#myfile <- file.path(rdatdir,paste(output_name,"tidy_peptide.csv",sep="_"))
#fwrite(tidy_peptide,myfile)

message("\nDone!")
