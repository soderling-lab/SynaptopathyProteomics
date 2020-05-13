#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Inputs:
# Input data should be in root/data/:
# 1. TMT-samples.csv - sample meta data.
input_samples = "4227_TMT_Cortex_Combined_traits.csv"
#input_samples = "4227_TMT_Striatum_Combined_traits.csv"

# 2. TMT-raw-peptide.csv - raw peptide data from PD.
input_data = "4227_TMT_Cortex_Combined_PD_Peptide_Intensity.csv"
#input_data = "4227_TMT_Striatum_Combined_PD_Peptide_Intensity.csv"

## Other parameters:
output_name = "Cortex" # Prefix for naming output files.
#output_name = "Striatum"
sample_connectivity_threshold = 2.5 # Threshold for detecting sample level outliers.
alpha = 0.1 # FDR threshold for differential abundance.

## Main Outputs:
# Stored in root/tables/
# 0. [output_name]_TMT_Results.xlsx

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. tidy_peptide.csv - tidy, raw peptide data.

## Order of data processing operations:
# * Load the data from PD.
# * Initial sample loading normalization -- done within an experiment.
# * Impute missing peptide values with KNN algorithm (k=10).
# * Examine QC peptide reproducibility -- remove peptide outliers.
# * Summarize to protein level by summing all peptides for a protein.
# * Protein-level sample loading normalization -- done across all experiments.
# * IRS normalization -- equalizes proteins quantified by different peptides.
# * Protein level filtering -- remove proteins identified by a single peptide;
#   remove proteins with too many missing values; remove proteins that are not 
#   reproducible (exhibit high inter-experimental variablility across the 3x 
#   biological replicates.
# * Reformat data for final TAMPOR normalization.

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. 

# Load renv -- use load NOT activate!
renv::load(getrd()) # NOTE: getrd() is a function in my .Rprofile.

# Load required packages and functions.
suppressPackageStartupMessages({
	library(dplyr)
	library(WGCNA)
	library(TBmiscr)
	library(getPPIs)
	library(ggplot2)
	library(data.table)
})

# Load additional functions.
TBmiscr::load_all()

# Project directories:
rootdir <- getrd()
funcdir <- file.path(rootdir, "R")
datadir <- file.path(rootdir, "data")
rdatdir <- file.path(rootdir, "rdata")
downdir <- file.path(rootdir, "downloads")
figsdir <- file.path(rootdir, "figs")
tabsdir <- file.path(rootdir, "tables")

#---------------------------------------------------------------------
## Load the raw data and sample info.
#---------------------------------------------------------------------

# Load the TMT data.
myfile <- file.path(datadir,input_data)
peptides <- fread(myfile)

# Load sample information.
myfile <- file.path(datadir,input_samples)
samples <- fread(myfile)

#---------------------------------------------------------------------
## Map all Uniprot accession numbers to stable entrez IDs.
#---------------------------------------------------------------------

# 0. Remove non-mouse and other spurious proteins.
# 1. Map Uniprot IDs to Entrez using MGI batch query function.
# 2. Map Entrez IDs to gene symbols using getPPIs::getIDs().
message("\nCreating gene identifier map.")

# First, remove any non-mouse proteins from the data.
peptides <- peptides %>% filter(grepl("OS=Mus musculus",Description))

# Remove these Immunoglobin proteins:
ig_prots <- c("P01631","P01646","P01665","P01680","P01746","P01750",
	      "P01786","P01864","P01878","P03975","P06330","P03987")
peptides <- peptides %>% filter(Accession %notin% ig_prots)

# Collect all uniprot IDs.
uniprot <- unique(peptides$Accession)

# Map Uniprot IDs to entrez using MGI database.
# This takes a couple minutes because the function currently
# downloads the MGI data each time the function is called.
entrez <- getPPIs::mgi_batch_query(uniprot,quiet=FALSE,download=FALSE)
names(entrez) <- uniprot

# Check: we have successfully mapped all uniprot ids.
check <- sum(is.na(entrez)) == 0
if (!check) { stop("Unable to map all Uniprot IDs to stable gene identifiers!") }

# Map entrez ids to gene symbols.
symbols <- getPPIs::getIDs(entrez,from="entrez",to="symbol",species="mouse")

# Check: we have successfully mapped all gene symbols.
check <- sum(is.na(symbols)) == 0
if (!check) { stop("Unable to map all Entrez IDs to gene Symbols!") }

# Create mapping data.table.
# We will save this in rdata/.
gene_map <- data.table(uniprot = names(entrez),
                       entrez = entrez,
	               symbol = symbols)
gene_map$id <- paste(gene_map$symbol,gene_map$uniprot,sep="|")

#---------------------------------------------------------------------
## Tidy-up the input data from Proteome Discover.
#---------------------------------------------------------------------

# Convert PD df into tidy df. 
# Samples should contain the following columns:
# Treatment, Channel, Sample, Experiment
message("\nLoading raw data from Proteome Discover...")
cols <- colnames(peptides)[!grepl("Abundance",colnames(peptides))]
tidy_peptide <- tidyProt(peptides,intensity.cols=cols)

# Annotate tidy data with additional meta data from samples.
tidy_peptide <- left_join(tidy_peptide,samples,by="Sample")

#---------------------------------------------------------------------
## Perform sample loading normalization.
#---------------------------------------------------------------------

# Perform sample normalization. Normalization is down for each 
# experiment independently (group by Experiment:Sample).
message("\nPerforming sample loading normalization.")
sl_peptide <- normSL(tidy_peptide, groupBy=c("Experiment","Sample"))

#---------------------------------------------------------------------
## Impute missing peptide values.
#---------------------------------------------------------------------

# Impute missing peptide values with KNN algorithm for MNAR data.
# * Missing QC values will not be imputed.
# * Peptides (rows) with more than 50% missingness will not be imputed.
# Values in these rows are masked (replaced with NA).
message("\nImputing a small number of missing peptide values.")
imputed_peptide <- imputeKNNpep(sl_peptide, groupBy="Genotype",
				samples_to_ignore="QC",quiet=FALSE) 

#---------------------------------------------------------------------
## Examine reproducibility of QC measurements.
#---------------------------------------------------------------------

# Examine reproducibility of QC measurements.
# This strategy was adapted from Ping et al., 2018 (pmid: 29533394).
# For each experiment, calculate the ratio of QC measurements.
# Bin these ratios based on the average Intensity of QC peptides into
# 5 bins. For each bin, remove measurements that are outside 
# +/- 4x SD from the mean of all ratios within that bin.
# This threshold can be set with nSD.
message("\nRemoving peptides with irreproducible QC measurements.")
filt_peptide <- filtQC(imputed_peptide, grouping.col='Treatment',
		       controls='QC',quiet=FALSE)

#---------------------------------------------------------------------
## Summarize to protein level.
#---------------------------------------------------------------------

# Protein intensity is just sum of all peptides for that protein.
message("\nSummarizing proteins as the sum of their peptides.")
proteins <- summarize_prot(filt_peptide)

# Perform SL normalization across all experiments (group by Sample).
message("\nPerforming sample loading normalization between experiments.")
sl_protein <- normSL(proteins, groupBy="Sample")

#---------------------------------------------------------------------
## Insure there are no QC outlier samples.
#---------------------------------------------------------------------

# Calculate Oldham's normalized sample connectivity (zK) in order to 
# identify outlier samples.
# This approach was adapted from Oldham et al., 2012 (pmid: 22691535).
zK <- sampleConnectivity(sl_protein %>% filter(Treatment == "QC"))

outlier_samples <- c(names(zK)[zK < -sample_connectivity_threshold],
		     names(zK)[zK > sample_connectivity_threshold])

# There are no sample outliers.
check <- length(outlier_samples) == 0
if (!check) { stop("Why are there outlier samples?") }

#---------------------------------------------------------------------
## Perform IRS Normalization.
#---------------------------------------------------------------------

# Equalize QC measurements between experiments. Adjusts protein 
# measurements of biological replicates simultaneously.
# Accounts for protein quantification by different peptides in 
# each experiment.
message(paste("\nStandardizing protein measurements between",
	      "experiments by IRS normalization."))
irs_protein <- normIRS(sl_protein,controls="QC",robust=TRUE)

#---------------------------------------------------------------------
## Protein level filtering.
#---------------------------------------------------------------------

# Remove proteins:
# Remove proteins that were identified by a single peptide.
# Remove proteins with too many missing values.
# Remove proteins with any missing QC values.
# FIXME: Remove proteins that have outlier measurements.
message(paste("\nFiltering proteins..."))
filt_protein <- filtProt(irs_protein,
			 controls="QC",nbins=5,nSD=4,summary=TRUE)

# There should be no missing values at this point.
check <- sum(is.na(filt_protein$Intensity)) == 0
if (!check) { stop("Why are there missing values!") }

#---------------------------------------------------------------------
## Identify Sample outliers.
#---------------------------------------------------------------------

# Calculate Oldham's normalized sample connectivity (zK) in order to 
# identify outlier samples.
# This approach was adapted from Oldham et al., 2012 (pmid: 22691535).
zK <- sampleConnectivity(filt_protein)
outlier_samples <- c(names(zK)[zK < -sample_connectivity_threshold],
		     names(zK)[zK > sample_connectivity_threshold])

# Status:
message(paste0("Outlier samples:\n",
	       paste(outlier_samples,collapse="\n")))

# Remove sample outliers.
filt_protein <- filt_protein %>% filter(Sample %notin% outlier_samples)

#---------------------------------------------------------------------
# Reformat data for TAMPOR normalization script.
#---------------------------------------------------------------------

# Reformat final normalized data for TAMPOR Normalization.
dm_tampor <- reformat_TAMPOR(filt_protein,samples)

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. [output_name]_tidy_peptide.csv -- raw peptide data.
# 2. [output_name]_cleanDat.RData" -- data for TAMPOR.

## Save key results.
message("\nSaving data for downstream analysis")

# 0. gene_map.RData   - gene identifier map.
myfile <- file.path(rdatdir,"gene_map.RData")
saveRDS(gene_map,myfile)

# 1. [output_name]_tidy_peptide.csv -- raw peptide data.
myfile <- file.path(rdatdir,paste(output_name,"tidy_peptide.csv",sep="_"))
fwrite(tidy_peptide,myfile)

# 2. [output_name]_cleanDat.RData" -- data for TAMPOR.
myfile <- file.path(rdatdir,paste(output_name,"cleanDat.RData",sep="_"))
fwrite(dm_tampor,myfile)
