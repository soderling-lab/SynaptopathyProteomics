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
input_data = list("Cortex" = "Cortex_norm_protein.csv",
		  "Striatum"= "Striatum_norm_protein.csv")

## Other parameters:
output_name = analysis_type
alpha_threshold = 0.05 # FDR threshold for signficance.
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

quit()

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

# Define a function that checks reproducibility of a protein.
check_reproducibility <- function(prot_list,protein,treatment.subset="WT",
				  fun="bicor",threshold=0.8) {
	# Which proteins are highly reproducible between replicates.
	# For a given protein, cast the data into a matrix: Fraction ~ Replicate.
	dt <- prot_list[[protein]] %>% 
		filter(Treatment == treatment.subset) %>%
		as.data.table() 
        dt_temp <- dcast(dt, Genotype ~ Treatment + Channel, 
		      value.var="Intensity") 
	value_cols <- grep("_",colnames(dt_temp))
	dm <- dt_temp %>% dplyr::select(all_of(value_cols)) %>% 
		as.matrix(rownames.value=dt_temp$Genotype)
	# missing values can arise in dm if outlier sample was removed.
	# Calculate the pairwise correlation between samples.
	opts <- "pairwise.complete.obs"
	if (fun == "bicor") { cormat <- WGCNA::bicor(log2(dm),
						     use=opts) }
	if (fun == "pearson") { cormat <- cor(log2(dm),method="pearson",
					      use=opts) }
	if (fun == "spearman") { cormat <- cor(log2(dm),method="spearman",
					       use=opts) }
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

# model = "~0 + group"
# model = "~background + group"

comparisons = c("Genotype","Treatment")
model = "~background + group"
data_glm <- glmDA(norm_protein, comparisons, model,value.var="Abundance")

# Extract statistical results from glm data.
glm_results <- data_glm$results

# Extract data from glm object.
glm_protein <- data_glm$data

#--------------------------------------------------------------------
## Tidy-up glm statistical results.
#--------------------------------------------------------------------

# We are interested in the four mutant mouse model contrasts.
glm_results <- glm_results[c((length(glm_results)-3):length(glm_results))]
names(glm_results) <- c("Shank3","Syngap1","Ube3a","Shank2")
glm_results <- glm_results[c("Shank2","Shank3","Syngap1","Ube3a")] # Sort.
glm_results <- lapply(glm_results,as.data.table)

# Summary of DA proteins.
message(paste0("Summary of differentially abundant proteins at FDR <",
	      alpha_threshold,":"))
tab <- sapply(glm_results,function(x) sum(as.numeric(x$FDR) < alpha_threshold))
knitr::kable(t(tab))

# Check top gene in each condition.
message("Top differentially abundant proteins:")
glm_summary <- sapply(glm_results,function(x) {
			      x %>% dplyr::select(Symbol) %>% head(3)
		  })
dt <- as.data.table(do.call(cbind,glm_summary))
colnames(dt) <- gsub(".Symbol"," Top3:",colnames(dt))
knitr::kable(dt)

# Merge glm_results by shared column names:
cols <- Reduce(intersect, lapply(glm_results,colnames))
glm_stats <- lapply(glm_results, function(x) {
			    as.data.table(x) %>% 
				    dplyr::select(all_of(cols)) }) %>%
				    bind_rows(.id="Genotype")

# Annotate with tissue type.
glm_stats <- tibble::add_column(glm_stats, Tissue = analysis_type,
				.after="Genotype")

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:
# * [output_name]_tidy_peptide.csv  - tidy, raw peptide data.
# * [output_name]_norm_protein      - final preprocessed data.
# * [output_name]_glm_stats.csv     - glm statistical results.
# * [output_name]_glm_protein.csv   - glm fitted protein values.
# * [output_name]_GLM_Results.xlsx  - glm statistical results.
# * [output_name]_repro_prots.RData - potential marker proteins.

## Save key results.
message("\nSaving data for downstream analysis.")

# [output_name]_tidy_peptide.csv -- raw peptide data.
myfile <- file.path(rdatdir,paste(output_name,"tidy_peptide.csv",sep="_"))
fwrite(tidy_peptide,myfile)

# [output_name]_norm_protein.csv -- final, normalized and regressed data.
myfile <- file.path(rdatdir,paste(output_name,"norm_protein.csv",sep="_"))
fwrite(norm_protein,myfile)

# [output_name]_glm_stats.csv -- tidy statistical results.
myfile <- file.path(rdatdir,paste(output_name,"glm_stats.csv",sep="_"))
fwrite(glm_stats,myfile)

# [output_name]_glm_protein.csv -- tidy glm fitted protein values.
myfile <- file.path(rdatdir,paste(output_name,"glm_protein.csv",sep="_"))
fwrite(glm_protein, myfile)

# [output_name]_GLM_Results.xlsx -- statistical results.
myfile <- file.path(tabsdir,paste(output_name,"GLM_Results.xlsx",sep="_"))
write_excel(glm_results,myfile)

# [output_name]_repro_prots.RData -- subset of highly reproducible prots.
myfile <- file.path(rdatdir,paste0(output_name,"_repro_prots.RData"))
saveRDS(reproducible_prots,myfile)

message("\nDone!")
