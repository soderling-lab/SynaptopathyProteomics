#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

analysis_type = "Combined"

## Optional parameters:
n_threads = parallel::detectCores()

## Input data in root/rdata/:
input_data = list(Cortex = "Cortex_norm_protein.csv",
		  Striatum = "Striatum_norm_protein.csv")

## Output for downstream analysis, stored in root/rdata/
output_name = analysis_type # Prefix for naming output files.

# Stored in root/rdata/
# * [output_name]_tidy_peptide.csv  - tidy, raw peptide data.
# * [output_name]_norm_protein      - final preprocessed data.

## Order of data processing operations:
# * Load the data from PD.
# * Initial sample loading normalization -- done within an experiment.
# * Impute missing peptide values with KNN algorithm (k=10).
# * Examine QC peptide reproducibility -- remove peptide outliers.
# * Summarize to protein level by summing all peptides for a protein.
# * Protein-level sample loading normalization -- done across all experiments.
# * Intra-batch ComBat.
# * IRS normalization -- equalizes proteins quantified by different peptides.
# * Protein level filtering -- remove proteins identified by a single peptide;
#      remove proteins with too many missing values; 
#      remove proteins that are not reproducible 
#      (exhibit high inter-experimental variablility across the 3x 
#      biological replicates.
# * Perform TAMPOR normalization in order to scales WT's to be equal.
# * Identify and remove sample outliers.
# * Fit a glm to assess protein differential abundance. 

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. 

# Load renv -- use load NOT activate!
root <- getrd() # See .Rprofile alias or TBmiscr::getrd().
renv::load(root)

# Load required packages and functions.
suppressPackageStartupMessages({
	library(dplyr)
	library(edgeR)
	library(getPPIs)
	library(data.table)
})

# Load additional functions.
TBmiscr::load_all()

# Project directories:
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

#---------------------------------------------------------------------
## Perform TAMPOR normalization to scale WT samples.
#---------------------------------------------------------------------

# Load the TMT data.
data_list <- lapply(input_data,function(x) fread(file.path(rdatdir,x)))
prot_dt <- bind_rows(data_list) %>% 
	filter(Treatment != "QC") %>% as.data.table()

# Define experimental batches:
batch <- interaction(prot_dt$Tissue,prot_dt$Genotype) %>% 
	as.numeric()
batch <- paste0("b",batch)

# Annotate data with batch and batch.channel (ID).
prot_dt$batch <- batch
prot_dt$ID <- interaction(batch,prot_dt$Channel) %>% 
	as.character()

# Cast into a data.matrix for TAMPOR.
dm <- prot_dt %>% dcast(Accession ~ ID, value.var="Intensity") %>% 
	as.matrix(rownames="Accession")

# Check: there should be no rows with missing values at this point.
missing_vals <- apply(dm,1,function(x) any(is.na(x)))
if (sum(missing_vals) > 0) {
	message(paste("\nRemoving", sum(missing_vals),
		      "rows that contain missing values."))
	dm <- dm[!missing_vals,]
}

# Define traits df for TAMPOR.
# Traits must include batch and ID columns.
traits <- prot_dt %>% dplyr::select(batch,ID) %>% unique()
traits <- as.data.frame(traits)
rownames(traits) <- traits[,"ID"] # rownames must be batch.channel.

# Define controls or global internal standards (GIS) for 
# TAMPOR normalization. Use all WT samples. 
controls <- unique(filter(prot_dt,Treatment=="WT")[["ID"]])

# Perform the normalization.
data_tampor <- TAMPOR(dat = dm, traits = traits, batchPrefixInSampleNames = TRUE,
		      samplesToIgnore = "none", GISchannels = controls, 
		      parallelThreads = n_threads)

# Collect normalized, relative abundance data post-TAMPOR.
# NOTE: Abundance is the normalized data.
suppressWarnings({ # Suppress warnings about coercing factors to characters.
tampor_protein <- data_tampor$cleanRelAbun %>% # Extract normalized data from list.
	reshape2::melt() %>% # Melt into a tidy df.
	setNames(c("Accession","ID","Abundance")) %>% # Set column names.
	as.data.table() %>% # Coerce to data.table.
	left_join(prot_dt,by=c("Accession","ID")) %>% # Combine with meta data from input.
	as.data.table()
})

# Clean-up -- remove unncessary cols.
tampor_protein$batch <- NULL
tampor_protein$ID <- NULL

#--------------------------------------------------------------------
## Identify and remove sample outliers.
#--------------------------------------------------------------------

# Calculate Oldham's normalized sample connectivity (zK) in order to 
# identify outlier samples.
# This approach was adapted from Oldham et al., 2012 (pmid: 22691535).

n_samples <- length(unique(tampor_protein$Sample))
message(paste("Initial number of samples:",n_samples))

# Loop to iteratively remove outliers.
i = 0
sample_connectivity_threshold = 2.5
outlier_samples <- NULL
while (i < 5) {
	data_in <- tampor_protein %>% 
		filter(Sample %notin% outlier_samples)
	zk <- sampleConnectivity(data_in,value.var="Abundance",log=TRUE)
	outlier_samples <- unique({
		c(outlier_samples,
		  names(zk)[zk < -sample_connectivity_threshold],
		  names(zk)[zk > sample_connectivity_threshold])
	})
	i = i + 1
}

# Remove sample outliers.
final_protein <- tampor_protein %>% filter(Sample %notin% outlier_samples)
n_samples <- length(unique(final_protein$Sample))

# Status:
if (length(outlier_samples) == 0) {
	message(paste("\nFinal number of samples:",
		      length(unique(final_protein$Sample))))
} else {
	message(paste0("\nRemoved outlier sample(s):\n",
	       paste(outlier_samples,collapse="\n")))
	message(paste("Final number of samples:",
		      length(unique(final_protein$Sample))))
}

#---------------------------------------------------------------------
## Evaluate protein differential abundance.
#---------------------------------------------------------------------
# Use EdgeR glm to evaluate differential protein abundance.
# Compare WT v Mutant.

# Status.
message(paste("\nAnalyzing protein differential abundance",
	      "with EdgeR GLM..."))

# Perform analysis WITH grouping of WT samples and fit a glm
# to account for differences in genetic background..

## FIXME: extract fitted data from glm object.
comparisons = c("Tissue","Genotype","Treatment")
treatment.ignore = "QC"
value.var = "Abundance"

glm_data <- glmDA(final_protein, treatment.ignore, 
		     value.var, comparisons,
		     combine.WT = FALSE)
glm_results <- glm_data$results

# Fit data is normally distributed.
fit_dt <- glm_data$fit.vals %>% filter(Treatment != "QC") %>% 
	as.data.table()

# Obs data is log-normal.
obs_dt <- glm_data$obs.vals %>% filter(Treatment != "QC") %>% 
	as.data.table()

# Combine.
cols <- intersect(colnames(fit_dt),colnames(obs_dt))
glm_protein <- left_join(fit_dt,obs_dt,by=cols)

# Annotate glm_results (stats list) with gene identifiers.
annotate_genes <- function(x) {
	idx <- match(x$Accession,gene_map$Accession)
	Symbol <- gene_map$Symbol[idx]
	Entrez <- gene_map$Entrez[idx]
	x <-  tibble::add_column(x,Symbol,.after="Accession")
	x <- tibble::add_column(x,Entrez,.after="Symbol")
	return(x)
}
glm_results <- lapply(glm_results,annotate_genes)

# Summary of DA proteins.
message(paste0("\nSummary of differentially abundant proteins at FDR <",
	      alpha_threshold,":"))
n <- sapply(glm_results,function(x) sum(x$FDR < alpha_threshold))
knitr::kable(as.data.table(n,keep.rownames="Genotype"))

# Check top gene in each condition.
message("\nTop differentially abundant proteins:")
knitr::kable(sapply(glm_results,function(x) head(x$Symbol,3)))

# Collect stats in a single df.
glm_stats <- bind_rows(glm_results,.id="Genotype")

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
prot_list <- tampor_protein %>% group_by(Accession) %>% group_split()
names(prot_list) <- sapply(prot_list,function(x) unique(x$Accession))

#ggplotPCA(glm_protein,value.var="Obs.Intensity",log=TRUE)

# Check protein reproducibility.
check <- list()
for (prot in names(prot_list)) {
	check[[prot]] <- check_reproducibility(prot_list, prot,
				       treatment.subset="WT",
				       value.var="Abundance",
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
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:

## Save key results.
message("\nSaving data for downstream analysis.")

# [output_name]_tidy_peptide.csv -- raw peptide data.
myfile <- file.path(rdatdir,paste(output_name,"tidy_peptide.csv",sep="_"))
fwrite(tidy_peptide,myfile)

# [output_name]_norm_protein.csv -- 
myfile <- file.path(rdatdir,paste(output_name,"norm_protein.csv",sep="_"))
fwrite(glm_protein,myfile)

# [output_name]_glm_stats.csv
myfile <- file.path(rdatdir,paste(output_name,"glm_stats.csv",sep="_"))
write_excel(glm_results,myfile)

# [output_name]_GLM_Results.xlsx
myfile <- file.path(tabsdir,paste(output_name,"GLM_Results.xlsx",sep="_"))
write_excel(glm_results,myfile)

message("\nDone!")
