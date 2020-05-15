#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

# Function to check command line input.
parse_args <- function(arg_opts = c("Cortex","Striatum"), nargs=1, err=NULL){
	args <- commandArgs(trailingOnly=TRUE)
	if (interactive()) { 
		set.seed(as.numeric(Sys.time()))
		n <- sample(c(1,2),1)
		warning(paste0("Running R interactively. ",
			       "Analyzing ",arg_opts[n]," data."))
		return(arg_opts[n])
	} else if (length(args) == nargs & args[nargs] %in% arg_opts) {
		return(args[nargs])
	} else {
		stop(err, call.=FALSE)
	}
}

## Parse input arguments.
msg <- "Please specify either 'Cortex' or 'Striatum' for the analysis."
analysis_type <- parse_args(err=msg)

## Inputs:
# Input data should be in root/data/:
# 1. TMT-samples.csv - sample meta data.
sample_files = list("Cortex" = "4227_TMT_Cortex_Combined_traits.csv",
		     "Striatum" = "4227_TMT_Striatum_Combined_traits.csv")
input_samples <- sample_files[[analysis_type]]

# 2. TMT-raw-peptide.csv - raw peptide data from PD.
data_files = list("Cortex" = "4227_TMT_Cortex_Combined_PD_Peptide_Intensity.csv",
		  "Striatum" = "4227_TMT_Striatum_Combined_PD_Peptide_Intensity.csv")
input_data <- data_files[[analysis_type]]

## Other parameters:
output_name = analysis_type # Prefix for naming output files.
sample_connectivity_threshold = 2.5 # Sample level outlier threshold. 

## Output for downstream analysis:
# Stored in root/rdata/
# 0. [output_name]_gene_map.RData     - gene identifier map.
# 1. [output_name]_tidy_peptide.csv   - tidy, raw peptide data.
# 2. [output_name]_preprocessed.RData - preprocessed data for TAMPOR.

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
# Alternatively, use: 
# renv::load("/mnt/d/projects/SynaptopathyProteomics") # Path to project root.

# Load required packages and functions.
suppressPackageStartupMessages({
	library(dplyr)
	library(TBmiscr)
	library(getPPIs)
	library(data.table)
})

# Load additional functions.
TBmiscr::load_all()

# Project directories:
rootdir <- getrd()
funcdir <- file.path(rootdir, "R")
datadir <- file.path(rootdir, "data")
figsdir <- file.path(rootdir, "figs")
rdatdir <- file.path(rootdir, "rdata")
tabsdir <- file.path(rootdir, "tables")
downdir <- file.path(rootdir, "downloads")

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

# Map any missing ids by hand.
not_mapped <- is.na(entrez)
if (sum(not_mapped) > 0) {
	message("Mapping missing IDs by hand.")
	mapped_by_hand <- c("P10853" = 319180)
	entrez[names(mapped_by_hand)] <- mapped_by_hand
}

# Check: no unmapped ids?
check <- sum(is.na(entrez)) == 0
if (!check) { stop( "Unable to map all uniprot IDs to stable entrez IDs.") }

# Check: we have successfully mapped all uniprot ids.
check <- sum(is.na(entrez)) == 0
if (!check) { stop("Unable to map all Uniprot IDs to stable gene identifiers!") }

# Map entrez ids to gene symbols.
symbols <- getPPIs::getIDs(entrez,from="entrez",to="symbol",species="mouse")

# Check: we have successfully mapped all gene symbols.
check <- sum(is.na(symbols)) == 0
if (!check) { stop("Unable to map all Entrez IDs to gene Symbols!") }

# Create mapping data.table.
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
# This strategy was adapted from Ping et al., 2018 (pmid:29533394).
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
## Intra-batch Protein-lavel ComBat.
#---------------------------------------------------------------------

# Before tackling inter-experimental batch effects, address the 
# intra-experimental batch effect.

# Each experimental cohort of 8 was prepared in two batches. 
# This was necessary because the ultra-centrifuge rotor used to 
# prepare purified synaptosomes holds a maximum of 6 samples. 
# This intra-batch batch effect was recorded for 6/8 batches. 
# Here I will utilize the sva::ComBat() function in order to account 
# for this batch effect before correcting for inter-batch batch
# effects with IRS normalization. 

# NOTE: In the absence of evidence of a batch effect 
# (not annotated or cor(PCA,batch)<0.1), then ComBat is not applied.
# If ComBat is not applied, then the data is returned un-regressed.

# NOTE: The values of QC samples are not adjusted by ComBat.
# QC samples were prepared from a seperate batch of mice and
# represent a single batch.

message("\nPerforming ComBat to remove intra-batch batch effect.")

# Perform ComBat for each dataset.
data_combat <- intrabatch_combat(sl_protein, samples,
				 group="Shank2", batch="PrepDate")
data_combat <- intrabatch_combat(data_combat, samples,
				 group="Shank3", batch="PrepDate")
data_combat <- intrabatch_combat(data_combat, samples,
				 group="Syngap1", batch="PrepDate")
data_combat <- intrabatch_combat(data_combat, samples,
				 group="Ube3a", batch="PrepDate")
norm_protein <- data_combat

#---------------------------------------------------------------------
## Insure there are no QC outlier samples.
#---------------------------------------------------------------------

# Remove QC outliers before performing IRS normalization.

# Calculate Oldham's normalized sample connectivity (zK) in order to 
# identify outlier samples.
# This approach was adapted from Oldham et al., 2012 (pmid: 22691535).
zK <- sampleConnectivity(norm_protein %>% filter(Treatment == "QC"))

outlier_samples <- c(names(zK)[zK < -sample_connectivity_threshold],
		     names(zK)[zK > sample_connectivity_threshold])

# There are no QC sample outliers.
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
irs_protein <- normIRS(norm_protein,controls="QC",robust=TRUE)

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

# Remove sample outliers.
filt_protein <- filt_protein %>% filter(Sample %notin% outlier_samples)

# Status:
# Good: by removing outlier protein measurments, we now have one less
# sample level outlier.
message(paste0("\nRemoving outlier sample(s):\n",
	       paste(outlier_samples,collapse="\n")))
message(paste("Final number of samples:",
	      length(unique(filt_protein$Sample))))

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. [output_name]_tidy_peptide.csv -- raw peptide data.
# 2. [output_name]_cleanDat.RData" -- data for TAMPOR.

## Save key results.
message("\nSaving data for downstream analysis...")

# 0. gene_map.RData   - gene identifier map.
myfile <- file.path(rdatdir,paste(output_name,"gene_map.RData",sep="_"))
saveRDS(gene_map,myfile)

# 1. [output_name]_tidy_peptide.csv -- raw peptide data.
myfile <- file.path(rdatdir,paste(output_name,"tidy_peptide.csv",sep="_"))
fwrite(tidy_peptide,myfile)

# 2. [output_name]_preprocessed.RData -- data for TAMPOR.
myfile <- file.path(rdatdir,paste(output_name,"preprocessed.RData",sep="_"))
saveRDS(filt_protein,myfile)

message("\nDone!")
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

# 2. TMT-raw-peptide.csv - raw peptide data from PD.
input_data = "4227_TMT_Cortex_Combined_PD_Peptide_Intensity.csv"

## Other parameters:
output_name = "Cortex" # Prefix for naming output files.
sample_connectivity_threshold = 2.5 # Threshold for detecting sample level outliers.
alpha = 0.1 # FDR threshold for differential abundance.

## Main Outputs:
# Stored in root/tables/
# 0. [output_name]_TMT_Results.xlsx

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. tidy_peptide.csv - raw peptide data.
# 2. tidy_protein.csv - final normalized protein data.
# 3. norm_protein.csv - final normalized protein data matrix.
# 4. adjm.csv         - adjacency matrix. 
# 5. ne_adjm.csv      - enhanced adjacency matrix. 
# 6. glm_stats.RData  - statistical results from edgeR glm

## Order of data processing operations:
# * Load the data from PD.
# * Initial sample loading normalization -- done within an experiment.
# * Impute missing peptide values with KNN algorithm (k=10).
# * Examine QC peptide reproducibility -- remove peptide outliers.
# * Summarize to protein level by summing all peptides for a protein.
# * Protein-level sample loading normalization -- done across all experiments.
# * IRS normalization -- equalizes protein measurements made from different
#   peptides.
# * Protein level filtering -- remove proteins identified by a single peptide;
#   remove proteins with too many missing values; remove proteins that are not 
#   reproducible (exhibit high inter-experimental variablility across the 3x 
#   biological replicates.
# * Regression of covariates --
# * Asses differential abundance with a general linear model 
#   ([Abundance] ~ 0 + groups) as implemented by the Edge R package 
#   and its functions gmQLFfit() and glmQLFTest().

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
## Final Normalization: Regression of covariates.
#---------------------------------------------------------------------

# Batch and sex are used as covariates.
eblm_protein <- eBLM_regression(filt_protein, 
				traits=samples, ignore="QC")

#---------------------------------------------------------------------
# Reformat data for TAMPOR normalization script.
#---------------------------------------------------------------------

# Reformat final normalized data for TAMPOR Normalization.
reformat_TAMPOR <- function(tp,samples){}

tp <- filt_protein
tp_in <- tp <- as.data.table(tp)
dt <- tp %>% dcast(Accession ~ Sample,value.var="Intensity")
dm <- as.matrix(dt,rownames="Accession")

# Row names are Symbol|Accession
idx <- match(rownames(dm),gene_map$uniprot)
rownames(dm) <- paste(gene_map$symbol[idx],rownames(dm),sep="|")

# Column names are batch.channel
idx <- match(colnames(dm), samples$Sample)
batch <- paste0("b",as.numeric(as.factor(samples$Genotype)))
channel <- samples$Channel
colnames(dm) <- paste(batch,channel,sep=".")[idx]

# Reorder based on batch.channel.
dm <- dm[,order(colnames(dm))]

## Save data as RData...
# Save raw data.
myfile <- file.path(Rdatadir, paste0(outputMatName, "_raw_peptide.RData"))
saveRDS(raw_peptide, myfile)

# Save cleanDat as RData.
myfile <- file.path(Rdatadir, paste0(outputMatName, "_cleanDat.RData"))
saveRDS(cleanDat, myfile)

quit()

#---------------------------------------------------------------------
## Protein differential abundance.
#---------------------------------------------------------------------

# Use EdgeR glm to evaluate differential protein abundance.
# Compare WT v Mutant within a fraction.
message(paste("\nEvaluating differential abundance between WT and",
	      "Swip mutant samples."))
comparisons <- "Genotype.Treatment"
results <- glmDA(eblm_protein,comparisons,samples,
		 samples_to_ignore="QC")

# Collect the results.
tidy_protein <- results$data
glm_results <- results$results

# Annotate glm_results with gene symbols and entrez ids.
annotate_genes <- function(dt,gene_map){
	idx <- match(dt$Accession,gene_map$uniprot)
	Entrez <- gene_map$entrez[idx]
	Symbol <- gene_map$symbol[idx]
	dt <- tibble::add_column(dt,Symbol,.after="Accession")
	dt <- tibble::add_column(dt,Entrez,.after="Symbol")
	return(dt)
}
glm_results <- lapply(glm_results,function(x) annotate_genes(x,gene_map))

# Summary of DA proteins:
message(paste0("Summary of differentially abundant proteins ",
	      "in each subceulluar fraction (FDR < ",alpha,"):"))
results$summary

#---------------------------------------------------------------------
## Examine PCA of final normalized data.
#---------------------------------------------------------------------

# Generate the plot.
# FIXME: customize colors!
plot <- ggplotPCA(tidy_protein)

# Save
myfile <- file.path(figsdir,"PCA.pdf")
ggsave(myfile,plot,height=5,width=5)

#---------------------------------------------------------------------
## Create protein co-expression matrix.
#---------------------------------------------------------------------

# Create data matrix with final normalized data.
norm_protein <- tidy_protein %>%
	dcast(Accession ~ Sample,value.var="Intensity") %>% 
	as.matrix(rownames=TRUE) %>% log2

# Remove QC samples when calculating co-expression matrices.
samples_to_ignore <- "SPQC"
idy <- grepl(samples_to_ignore,colnames(norm_protein))
subdat <- norm_protein[,!idy]

# Get list of data: all, just wt, just mutant.
samples_to_keep <- c(All="",WT="Control",KO="Mutant")
data_list <- sapply(samples_to_keep, function(x) { 
			    subdat[,grep(x,colnames(subdat))]
		 })

# Calculate bicor correlation matrix for all subsets of the data.
# Silence suppress unwanted output generated by WGNA::bicor().
message("\nCreating protein co-expression network.")
silence({ adjm_list <- lapply(data_list,function(x) bicor(t(x))) })

# Perform network enhancment for all matrices.
message(paste("\nPerforming network enhancement to denoise network."))
ne_list <- lapply(adjm_list, function(x) neten::neten(x))

#---------------------------------------------------------------------
## Save the TMT data and stats as a single excel workbook
#---------------------------------------------------------------------

## Save TMT data as a single excel document...
## Create an excel workbook with the following sheets:
# * Samples
# * Input - raw peptides
# * Output - normalized protein
# * Statistical results:
#     - Seperate sheet for each fraction/comparison.
#     - Include contrast specific data.

# Sort the statistical results.
glm_results <- glm_results[c("F4","F5","F6","F7","F8","F9","F10")]

# Loop to add normalized protein data to glm statistical results.
message("\nSaving TMT data and statistical results.")
for (i in 1:length(glm_results)){
	df <- glm_results[[i]]
	namen <- names(glm_results)[i]
	subsamples <- samples$Sample[grepl(namen,samples$Fraction)]
	idy <- match(subsamples,colnames(norm_protein))
	dm <- norm_protein[,idy]
	# Sort the data by Exp.Sample
	ids <- sapply(sapply(strsplit(colnames(dm),"\\."),"[",c(6,7),
			     simplify=FALSE), paste,collapse=".")
	names(ids) <- colnames(dm)
	ids <- ids[order(ids)]
	dm <- dm[,names(ids)]
	# Coerce dm to dt and merge with stats df.
	dt <- as.data.table(dm,keep.rownames="Accession")
	dt_out <- left_join(df,dt,by="Accession") %>% as.data.table
	glm_results[[i]] <- dt_out
} # Ends loop.

# Map uniprot to entrez IDs.
uniprot <- glm_results[[1]][["Accession"]]
idx <- match(uniprot,gene_map$uniprot)
entrez <- gene_map$entrez[idx]
names(entrez) <- uniprot

# Map missing ids by hand.
is_missing <- is.na(entrez)
missing_entrez <- entrez[is_missing]
mapped_by_hand <- c("P00848"=17705,
		    "P62806"=326619,
		    "Q64525"=319189,
		    "P06330"=111507, 
		    "Q5PR69"=320827, 
		    "Q80WG5"=241296,
		    "P03975"=15598, 
		    "Q80TK0"=100503185)
entrez[names(mapped_by_hand)] <- mapped_by_hand

# Check that we have successfully mapped all uniprot ids to entrez.
check <- sum(is.na(entrez))
if (check != 0) { stop("There are unmapped uniprot IDs.") }

# Map entrez to gene symbols.
idx <- match(entrez,gene_map$entrez)
symbols <- gene_map$symbol[idx]
names(symbols) <- entrez

# Loop to add Entrez IDs and Gene Symbols to data.
# Also fix column titles at the same time. 
for (i in 1:length(glm_results)){
	df <- glm_results[[i]]
	colnames(df) <- gsub("\\."," ",colnames(df))
	idy <- grepl("Abundance",colnames(df))
	colnames(df)[idy] <- paste("Log2",colnames(df)[idy])
	df <- tibble::add_column(df,"Entrez"=entrez[df$Accession],
				 .after="Accession")
	df <- tibble::add_column(df,"Gene"=symbols[as.character(df$Entrez)],
				 .after="Entrez")
	glm_results[[i]] <- df
}

## Save as excel workboook.
names(glm_results) <- paste(names(glm_results),"Results")
results <- c(list("Samples" = samples,
		"Normalized Protein" = tidy_protein),
	     glm_results)
myfile <- file.path(tabsdir,"Swip_TMT_Results.xlsx")
write_excel(results,myfile,rowNames=FALSE)

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. tidy_peptide.csv - raw peptide data.
# 2. tidy_protein.csv - final normalized protein data.
# 3. norm_protein.csv - final normalized protein data matrix.
# 4. adjm.csv         - adjacency matrix. 
# 5. ne_adjm.csv      - enhanced adjacency matrix. 
# 6. glm_stats.RData  - statistical results from edgeR glm

## Save key results.
message("\nSaving data for downstream analysis")

# 0. Save gene map.
myfile <- file.path(rdatdir,"gene_map.RData")
saveRDS(gene_map,myfile)

# 1. Save tidy_peptide.
myfile <- file.path(rdatdir,"tidy_peptide.csv")
fwrite(tidy_peptide,myfile)

# 2. Save tidy_protein.
myfile <- file.path(rdatdir,"tidy_protein.csv")
fwrite(tidy_protein,myfile)

# 3. Save normalized protein data matrix.
myfile <- file.path(rdatdir,"norm_protein.csv")
norm_protein %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(myfile)

# Save WT subset of data matrix.
myfile <- file.path(rdatdir,"norm_wt_protein.csv")
data_list[["WT"]] %>% as.data.table(keep.rownames="Accession") %>%
	       fwrite(myfile)
myfile <- file.path(rdatdir,"norm_ko_protein.csv")

# Save KO subset of data matrix.
data_list[["KO"]] %>% as.data.table(keep.rownames="Accession") %>%
	       fwrite(myfile)

# 4. Save adjms.
myfile <- file.path(rdatdir,"adjm.csv")
adjm_list[["All"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

myfile <- file.path(rdatdir,"wt_adjm.csv")
adjm_list[["WT"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

myfile <- file.path(rdatdir,"ko_adjm.csv")
adjm_list[["KO"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

# 5. Save enhanced adjms.
myfile <- file.path(rdatdir,"ne_adjm.csv")
ne_list[["All"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

myfile <- file.path(rdatdir,"wt_ne_adjm.csv")
ne_list[["WT"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

myfile <- file.path(rdatdir,"ko_ne_adjm.csv")
ne_list[["KO"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)
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

## Main Outputs:
# Stored in root/tables/
# 0. [output_name]_TMT_Results.xlsx

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. tidy_peptide.csv - tidy, raw peptide data.
# 2. [output_name]_cleanDat.RData - preprocessed data for TAMPOR.

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
# Alternatively, use: 
# renv::load("/mnt/d/projects/SynaptopathyProteomics") # Path to project root.

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

# Remove QC outliers before performing IRS normalization.

# Calculate Oldham's normalized sample connectivity (zK) in order to 
# identify outlier samples.
# This approach was adapted from Oldham et al., 2012 (pmid: 22691535).
zK <- sampleConnectivity(sl_protein %>% filter(Treatment == "QC"))

outlier_samples <- c(names(zK)[zK < -sample_connectivity_threshold],
		     names(zK)[zK > sample_connectivity_threshold])

# There are no QC sample outliers.
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
# Good: by removing outlier protein measurments, we now have one less
# sample level outlier.
message(paste0("\nRemoving outlier sample(s):\n",
	       paste(outlier_samples,collapse="\n")))

# Remove sample outliers.
filt_protein <- filt_protein %>% filter(Sample %notin% outlier_samples)

#---------------------------------------------------------------------
# Reformat data for TAMPOR normalization script.
#---------------------------------------------------------------------

# Reformat final normalized data for TAMPOR Normalization.
cleanDat <- reformat_TAMPOR(filt_protein,samples)

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. [output_name]_tidy_peptide.csv -- raw peptide data.
# 2. [output_name]_cleanDat.RData" -- data for TAMPOR.

## Save key results.
message("\nSaving data for downstream analysis...")

# 0. gene_map.RData   - gene identifier map.
myfile <- file.path(rdatdir,"gene_map.RData")
saveRDS(gene_map,myfile)

# 1. [output_name]_tidy_peptide.csv -- raw peptide data.
myfile <- file.path(rdatdir,paste(output_name,"tidy_peptide.csv",sep="_"))
fwrite(tidy_peptide,myfile)

# 2. [output_name]_cleanDat.RData -- data for TAMPOR.
myfile <- file.path(rdatdir,paste(output_name,"cleanDat.RData",sep="_"))
saveRDS(cleanDat,myfile)

message("\nDone!")
* Project '/mnt/d/projects/SynaptopathyProteomics' loaded. [renv 0.10.0]
Loading functions in SynaptopathyProteomics/R

Creating gene identifier map.

Loading raw data from Proteome Discover...

Performing sample loading normalization.

Imputing a small number of missing peptide values.
Imputing missing values in group: Shank2
... Percent missing values: 0.014% (n=71).
Imputing missing values in group: Shank3
... Percent missing values: 0.02% (n=103).
Imputing missing values in group: Syngap1
... Percent missing values: 0.039% (n=202).
Imputing missing values in group: Ube3a
... Percent missing values: 0.009% (n=48).

Removing peptides with irreproducible QC measurements.
Total number of QC outlier peptides identified: 611

Summarizing proteins as the sum of their peptides.

Performing sample loading normalization between experiments.

Performing ComBat to remove intra-batch batch effect.
Warning message:
For group: Shank2 - No detectable batch effect! Not applying ComBat. 
Initial coorelation between batch and samples: 0.828
Performing ComBat for group: Shank3
Final coorelation between batch and samples: 0.102 

Initial coorelation between batch and samples: 0.543
Performing ComBat for group: Syngap1
Final coorelation between batch and samples: 0.114 

Warning message:
For group: Ube3a - No detectable batch effect! Not applying ComBat. 

Standardizing protein measurements between experiments by IRS normalization.

Filtering proteins...
Removed proteins identified by one peptide.....: 95
Removed proteins with too many missing values..: 52
Removed proteins with any missing QC values....: 430
Final number of reproducibly quantified proteins: 2,860

Removing outlier sample(s):
Abundance: F4: 127C, Sample, KO, 41335, Shank3, Cortex
Final number of samples: 43

Saving data for downstream analysis...

Done!
#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

# Function to check command line input.
parse_args <- function(arg_opts = c("Cortex","Striatum"), nargs=1, err=NULL){
	args <- commandArgs(trailingOnly=TRUE)
	if (interactive()) { 
		set.seed(as.numeric(Sys.time()))
		n <- sample(c(1,2),1)
		warning(paste0("Running R interactively. ",
			       "Analyzing ",arg_opts[n]," data."))
		return(arg_opts[n])
	} else if (length(args) == nargs & args[nargs] %in% arg_opts) {
		return(args[nargs])
	} else {
		stop(err, call.=FALSE)
	}
}

## Parse input arguments.
msg <- "Please specify either 'Cortex' or 'Striatum' for the analysis."
analysis_type <- parse_args(err=msg)

## Inputs:
# Input data should be in root/data/:
# 1. TMT-samples.csv - sample meta data.
sample_files = list("Cortex" = "4227_TMT_Cortex_Combined_traits.csv",
		     "Striatum" = "4227_TMT_Striatum_Combined_traits.csv")
input_samples <- sample_files[[analysis_type]]

# 2. TMT-raw-peptide.csv - raw peptide data from PD.
data_files = list("Cortex" = "4227_TMT_Cortex_Combined_PD_Peptide_Intensity.csv",
		  "Striatum" = "4227_TMT_Striatum_Combined_PD_Peptide_Intensity.csv")
input_data <- data_files[[analysis_type]]

## Other parameters:
output_name = analysis_type # Prefix for naming output files.
sample_connectivity_threshold = 2.5 # Sample level outlier threshold. 

## Output for downstream analysis:
# Stored in root/rdata/
# 0. [output_name]_gene_map.RData     - gene identifier map.
# 1. [output_name]_tidy_peptide.csv   - tidy, raw peptide data.
# 2. [output_name]_preprocessed.RData - preprocessed data for TAMPOR.

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
# Alternatively, use: 
# renv::load("/mnt/d/projects/SynaptopathyProteomics") # Path to project root.

# Load required packages and functions.
suppressPackageStartupMessages({
	library(dplyr)
	library(TBmiscr)
	library(getPPIs)
	library(data.table)
})

# Load additional functions.
TBmiscr::load_all()

# Project directories:
rootdir <- getrd()
funcdir <- file.path(rootdir, "R")
datadir <- file.path(rootdir, "data")
figsdir <- file.path(rootdir, "figs")
rdatdir <- file.path(rootdir, "rdata")
tabsdir <- file.path(rootdir, "tables")
downdir <- file.path(rootdir, "downloads")

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

# Map any missing ids by hand.
not_mapped <- is.na(entrez)
if (sum(not_mapped) > 0) {
	message("Mapping missing IDs by hand.")
	mapped_by_hand <- c("P10853" = 319180)
	entrez[names(mapped_by_hand)] <- mapped_by_hand
}

# Check: no unmapped ids?
check <- sum(is.na(entrez)) == 0
if (!check) { stop( "Unable to map all uniprot IDs to stable entrez IDs.") }

# Check: we have successfully mapped all uniprot ids.
check <- sum(is.na(entrez)) == 0
if (!check) { stop("Unable to map all Uniprot IDs to stable gene identifiers!") }

# Map entrez ids to gene symbols.
symbols <- getPPIs::getIDs(entrez,from="entrez",to="symbol",species="mouse")

# Check: we have successfully mapped all gene symbols.
check <- sum(is.na(symbols)) == 0
if (!check) { stop("Unable to map all Entrez IDs to gene Symbols!") }

# Create mapping data.table.
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
# This strategy was adapted from Ping et al., 2018 (pmid:29533394).
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
## Intra-batch Protein-lavel ComBat.
#---------------------------------------------------------------------

# Before tackling inter-experimental batch effects, address the 
# intra-experimental batch effect.

# Each experimental cohort of 8 was prepared in two batches. 
# This was necessary because the ultra-centrifuge rotor used to 
# prepare purified synaptosomes holds a maximum of 6 samples. 
# This intra-batch batch effect was recorded for 6/8 batches. 
# Here I will utilize the sva::ComBat() function in order to account 
# for this batch effect before correcting for inter-batch batch
# effects with IRS normalization. 

# NOTE: In the absence of evidence of a batch effect 
# (not annotated or cor(PCA,batch)<0.1), then ComBat is not applied.
# If ComBat is not applied, then the data is returned un-regressed.

# NOTE: The values of QC samples are not adjusted by ComBat.
# QC samples were prepared from a seperate batch of mice and
# represent a single batch.

message("\nPerforming ComBat to remove intra-batch batch effect.")

# Perform ComBat for each dataset.
data_combat <- intrabatch_combat(sl_protein, samples,
				 group="Shank2", batch="PrepDate")
data_combat <- intrabatch_combat(data_combat, samples,
				 group="Shank3", batch="PrepDate")
data_combat <- intrabatch_combat(data_combat, samples,
				 group="Syngap1", batch="PrepDate")
data_combat <- intrabatch_combat(data_combat, samples,
				 group="Ube3a", batch="PrepDate")
norm_protein <- data_combat

#---------------------------------------------------------------------
## Insure there are no QC outlier samples.
#---------------------------------------------------------------------

# Remove QC outliers before performing IRS normalization.

# Calculate Oldham's normalized sample connectivity (zK) in order to 
# identify outlier samples.
# This approach was adapted from Oldham et al., 2012 (pmid: 22691535).
zK <- sampleConnectivity(norm_protein %>% filter(Treatment == "QC"))

outlier_samples <- c(names(zK)[zK < -sample_connectivity_threshold],
		     names(zK)[zK > sample_connectivity_threshold])

# There are no QC sample outliers.
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
irs_protein <- normIRS(norm_protein,controls="QC",robust=TRUE)

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

# Remove sample outliers.
filt_protein <- filt_protein %>% filter(Sample %notin% outlier_samples)

# Status:
# Good: by removing outlier protein measurments, we now have one less
# sample level outlier.
message(paste0("\nRemoving outlier sample(s):\n",
	       paste(outlier_samples,collapse="\n")))
message(paste("Final number of samples:",
	      length(unique(filt_protein$Sample))))

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. [output_name]_tidy_peptide.csv -- raw peptide data.
# 2. [output_name]_cleanDat.RData" -- data for TAMPOR.

## Save key results.
message("\nSaving data for downstream analysis...")

# 0. gene_map.RData   - gene identifier map.
myfile <- file.path(rdatdir,paste(output_name,"gene_map.RData",sep="_"))
saveRDS(gene_map,myfile)

# 1. [output_name]_tidy_peptide.csv -- raw peptide data.
myfile <- file.path(rdatdir,paste(output_name,"tidy_peptide.csv",sep="_"))
fwrite(tidy_peptide,myfile)

# 2. [output_name]_preprocessed.RData -- data for TAMPOR.
myfile <- file.path(rdatdir,paste(output_name,"preprocessed.RData",sep="_"))
saveRDS(filt_protein,myfile)

message("\nDone!")
* Project '/mnt/d/projects/SynaptopathyProteomics' loaded. [renv 0.10.0]
Loading functions in SynaptopathyProteomics/R

Creating gene identifier map.
Warning message:
In getPPIs::mgi_batch_query(uniprot, quiet = FALSE, download = FALSE) :
  Unable to map 1 gene identifiers to entrez.
Mapping missing IDs by hand.

Loading raw data from Proteome Discover...

Performing sample loading normalization.

Imputing a small number of missing peptide values.
Imputing missing values in group: Shank2
... Percent missing values: 0.029% (n=148).
Imputing missing values in group: Shank3
... Percent missing values: 0.014% (n=73).
Imputing missing values in group: Syngap1
... Percent missing values: 0.017% (n=86).
Imputing missing values in group: Ube3a
... Percent missing values: 0.02% (n=104).

Removing peptides with irreproducible QC measurements.
Total number of QC outlier peptides identified: 423

Summarizing proteins as the sum of their peptides.

Performing sample loading normalization between experiments.

Performing ComBat to remove intra-batch batch effect.
Initial coorelation between batch and samples: 0.588
Performing ComBat for group: Shank2
Final coorelation between batch and samples: 0.04 

Warning message:
For group: Shank3 - Cannot perform ComBat with one batch! 
Warning message:
For group: Syngap1 - Cannot perform ComBat with one batch! 
Warning message:
For group: Ube3a - No detectable batch effect! Not applying ComBat. 

Standardizing protein measurements between experiments by IRS normalization.

Filtering proteins...
Removed proteins identified by one peptide.....: 124
Removed proteins with too many missing values..: 51
Removed proteins with any missing QC values....: 605
Final number of reproducibly quantified proteins: 2,858

Removing outlier sample(s):
Abundance: F2: 128C, Sample, KO, 38894, Ube3a, Striatum
Abundance: F2: 131N, Sample, KO, 38898, Ube3a, Striatum
Final number of samples: 42

Saving data for downstream analysis...

Done!
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

# 2. TMT-raw-peptide.csv - raw peptide data from PD.
input_data = "4227_TMT_Cortex_Combined_PD_Peptide_Intensity.csv"

## Other parameters:
output_name = "Cortex" # Prefix for naming output files.
sample_connectivity_threshold = 2.5 # Threshold for detecting sample level outliers.
alpha = 0.1 # FDR threshold for differential abundance.

## Main Outputs:
# Stored in root/tables/
# 0. [output_name]_TMT_Results.xlsx

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. tidy_peptide.csv - raw peptide data.
# 2. tidy_protein.csv - final normalized protein data.
# 3. norm_protein.csv - final normalized protein data matrix.
# 4. adjm.csv         - adjacency matrix. 
# 5. ne_adjm.csv      - enhanced adjacency matrix. 
# 6. glm_stats.RData  - statistical results from edgeR glm

## Order of data processing operations:
# * Load the data from PD.
# * Initial sample loading normalization -- done within an experiment.
# * Impute missing peptide values with KNN algorithm (k=10).
# * Examine QC peptide reproducibility -- remove peptide outliers.
# * Summarize to protein level by summing all peptides for a protein.
# * Protein-level sample loading normalization -- done across all experiments.
# * IRS normalization -- equalizes protein measurements made from different
#   peptides.
# * Protein level filtering -- remove proteins identified by a single peptide;
#   remove proteins with too many missing values; remove proteins that are not 
#   reproducible (exhibit high inter-experimental variablility across the 3x 
#   biological replicates.
# * Regression of covariates --
# * Asses differential abundance with a general linear model 
#   ([Abundance] ~ 0 + groups) as implemented by the Edge R package 
#   and its functions gmQLFfit() and glmQLFTest().

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
## Final Normalization: Regression of covariates.
#---------------------------------------------------------------------

# Batch and sex are used as covariates.
eblm_protein <- eBLM_regression(filt_protein, 
				traits=samples, ignore="QC")

#---------------------------------------------------------------------
# Reformat data for TAMPOR normalization script.
#---------------------------------------------------------------------

# Reformat final normalized data for TAMPOR Normalization.
reformat_TAMPOR <- function(tp,samples){}

tp <- filt_protein
tp_in <- tp <- as.data.table(tp)
dt <- tp %>% dcast(Accession ~ Sample,value.var="Intensity")
dm <- as.matrix(dt,rownames="Accession")

# Row names are Symbol|Accession
idx <- match(rownames(dm),gene_map$uniprot)
rownames(dm) <- paste(gene_map$symbol[idx],rownames(dm),sep="|")

# Column names are batch.channel
idx <- match(colnames(dm), samples$Sample)
batch <- paste0("b",as.numeric(as.factor(samples$Genotype)))
channel <- samples$Channel
colnames(dm) <- paste(batch,channel,sep=".")[idx]

# Reorder based on batch.channel.
dm <- dm[,order(colnames(dm))]

## Save data as RData...
# Save raw data.
myfile <- file.path(Rdatadir, paste0(outputMatName, "_raw_peptide.RData"))
saveRDS(raw_peptide, myfile)

# Save cleanDat as RData.
myfile <- file.path(Rdatadir, paste0(outputMatName, "_cleanDat.RData"))
saveRDS(cleanDat, myfile)

quit()

#---------------------------------------------------------------------
## Protein differential abundance.
#---------------------------------------------------------------------

# Use EdgeR glm to evaluate differential protein abundance.
# Compare WT v Mutant within a fraction.
message(paste("\nEvaluating differential abundance between WT and",
	      "Swip mutant samples."))
comparisons <- "Genotype.Treatment"
results <- glmDA(eblm_protein,comparisons,samples,
		 samples_to_ignore="QC")

# Collect the results.
tidy_protein <- results$data
glm_results <- results$results

# Annotate glm_results with gene symbols and entrez ids.
annotate_genes <- function(dt,gene_map){
	idx <- match(dt$Accession,gene_map$uniprot)
	Entrez <- gene_map$entrez[idx]
	Symbol <- gene_map$symbol[idx]
	dt <- tibble::add_column(dt,Symbol,.after="Accession")
	dt <- tibble::add_column(dt,Entrez,.after="Symbol")
	return(dt)
}
glm_results <- lapply(glm_results,function(x) annotate_genes(x,gene_map))

# Summary of DA proteins:
message(paste0("Summary of differentially abundant proteins ",
	      "in each subceulluar fraction (FDR < ",alpha,"):"))
results$summary

#---------------------------------------------------------------------
## Examine PCA of final normalized data.
#---------------------------------------------------------------------

# Generate the plot.
# FIXME: customize colors!
plot <- ggplotPCA(tidy_protein)

# Save
myfile <- file.path(figsdir,"PCA.pdf")
ggsave(myfile,plot,height=5,width=5)

#---------------------------------------------------------------------
## Create protein co-expression matrix.
#---------------------------------------------------------------------

# Create data matrix with final normalized data.
norm_protein <- tidy_protein %>%
	dcast(Accession ~ Sample,value.var="Intensity") %>% 
	as.matrix(rownames=TRUE) %>% log2

# Remove QC samples when calculating co-expression matrices.
samples_to_ignore <- "SPQC"
idy <- grepl(samples_to_ignore,colnames(norm_protein))
subdat <- norm_protein[,!idy]

# Get list of data: all, just wt, just mutant.
samples_to_keep <- c(All="",WT="Control",KO="Mutant")
data_list <- sapply(samples_to_keep, function(x) { 
			    subdat[,grep(x,colnames(subdat))]
		 })

# Calculate bicor correlation matrix for all subsets of the data.
# Silence suppress unwanted output generated by WGNA::bicor().
message("\nCreating protein co-expression network.")
silence({ adjm_list <- lapply(data_list,function(x) bicor(t(x))) })

# Perform network enhancment for all matrices.
message(paste("\nPerforming network enhancement to denoise network."))
ne_list <- lapply(adjm_list, function(x) neten::neten(x))

#---------------------------------------------------------------------
## Save the TMT data and stats as a single excel workbook
#---------------------------------------------------------------------

## Save TMT data as a single excel document...
## Create an excel workbook with the following sheets:
# * Samples
# * Input - raw peptides
# * Output - normalized protein
# * Statistical results:
#     - Seperate sheet for each fraction/comparison.
#     - Include contrast specific data.

# Sort the statistical results.
glm_results <- glm_results[c("F4","F5","F6","F7","F8","F9","F10")]

# Loop to add normalized protein data to glm statistical results.
message("\nSaving TMT data and statistical results.")
for (i in 1:length(glm_results)){
	df <- glm_results[[i]]
	namen <- names(glm_results)[i]
	subsamples <- samples$Sample[grepl(namen,samples$Fraction)]
	idy <- match(subsamples,colnames(norm_protein))
	dm <- norm_protein[,idy]
	# Sort the data by Exp.Sample
	ids <- sapply(sapply(strsplit(colnames(dm),"\\."),"[",c(6,7),
			     simplify=FALSE), paste,collapse=".")
	names(ids) <- colnames(dm)
	ids <- ids[order(ids)]
	dm <- dm[,names(ids)]
	# Coerce dm to dt and merge with stats df.
	dt <- as.data.table(dm,keep.rownames="Accession")
	dt_out <- left_join(df,dt,by="Accession") %>% as.data.table
	glm_results[[i]] <- dt_out
} # Ends loop.

# Map uniprot to entrez IDs.
uniprot <- glm_results[[1]][["Accession"]]
idx <- match(uniprot,gene_map$uniprot)
entrez <- gene_map$entrez[idx]
names(entrez) <- uniprot

# Map missing ids by hand.
is_missing <- is.na(entrez)
missing_entrez <- entrez[is_missing]
mapped_by_hand <- c("P00848"=17705,
		    "P62806"=326619,
		    "Q64525"=319189,
		    "P06330"=111507, 
		    "Q5PR69"=320827, 
		    "Q80WG5"=241296,
		    "P03975"=15598, 
		    "Q80TK0"=100503185)
entrez[names(mapped_by_hand)] <- mapped_by_hand

# Check that we have successfully mapped all uniprot ids to entrez.
check <- sum(is.na(entrez))
if (check != 0) { stop("There are unmapped uniprot IDs.") }

# Map entrez to gene symbols.
idx <- match(entrez,gene_map$entrez)
symbols <- gene_map$symbol[idx]
names(symbols) <- entrez

# Loop to add Entrez IDs and Gene Symbols to data.
# Also fix column titles at the same time. 
for (i in 1:length(glm_results)){
	df <- glm_results[[i]]
	colnames(df) <- gsub("\\."," ",colnames(df))
	idy <- grepl("Abundance",colnames(df))
	colnames(df)[idy] <- paste("Log2",colnames(df)[idy])
	df <- tibble::add_column(df,"Entrez"=entrez[df$Accession],
				 .after="Accession")
	df <- tibble::add_column(df,"Gene"=symbols[as.character(df$Entrez)],
				 .after="Entrez")
	glm_results[[i]] <- df
}

## Save as excel workboook.
names(glm_results) <- paste(names(glm_results),"Results")
results <- c(list("Samples" = samples,
		"Normalized Protein" = tidy_protein),
	     glm_results)
myfile <- file.path(tabsdir,"Swip_TMT_Results.xlsx")
write_excel(results,myfile,rowNames=FALSE)

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. tidy_peptide.csv - raw peptide data.
# 2. tidy_protein.csv - final normalized protein data.
# 3. norm_protein.csv - final normalized protein data matrix.
# 4. adjm.csv         - adjacency matrix. 
# 5. ne_adjm.csv      - enhanced adjacency matrix. 
# 6. glm_stats.RData  - statistical results from edgeR glm

## Save key results.
message("\nSaving data for downstream analysis")

# 0. Save gene map.
myfile <- file.path(rdatdir,"gene_map.RData")
saveRDS(gene_map,myfile)

# 1. Save tidy_peptide.
myfile <- file.path(rdatdir,"tidy_peptide.csv")
fwrite(tidy_peptide,myfile)

# 2. Save tidy_protein.
myfile <- file.path(rdatdir,"tidy_protein.csv")
fwrite(tidy_protein,myfile)

# 3. Save normalized protein data matrix.
myfile <- file.path(rdatdir,"norm_protein.csv")
norm_protein %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(myfile)

# Save WT subset of data matrix.
myfile <- file.path(rdatdir,"norm_wt_protein.csv")
data_list[["WT"]] %>% as.data.table(keep.rownames="Accession") %>%
	       fwrite(myfile)
myfile <- file.path(rdatdir,"norm_ko_protein.csv")

# Save KO subset of data matrix.
data_list[["KO"]] %>% as.data.table(keep.rownames="Accession") %>%
	       fwrite(myfile)

# 4. Save adjms.
myfile <- file.path(rdatdir,"adjm.csv")
adjm_list[["All"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

myfile <- file.path(rdatdir,"wt_adjm.csv")
adjm_list[["WT"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

myfile <- file.path(rdatdir,"ko_adjm.csv")
adjm_list[["KO"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

# 5. Save enhanced adjms.
myfile <- file.path(rdatdir,"ne_adjm.csv")
ne_list[["All"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

myfile <- file.path(rdatdir,"wt_ne_adjm.csv")
ne_list[["WT"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

myfile <- file.path(rdatdir,"ko_ne_adjm.csv")
ne_list[["KO"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)
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

## Main Outputs:
# Stored in root/tables/
# 0. [output_name]_TMT_Results.xlsx

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. tidy_peptide.csv - tidy, raw peptide data.
# 2. [output_name]_cleanDat.RData - preprocessed data for TAMPOR.

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
# Alternatively, use: 
# renv::load("/mnt/d/projects/SynaptopathyProteomics") # Path to project root.

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

# Remove QC outliers before performing IRS normalization.

# Calculate Oldham's normalized sample connectivity (zK) in order to 
# identify outlier samples.
# This approach was adapted from Oldham et al., 2012 (pmid: 22691535).
zK <- sampleConnectivity(sl_protein %>% filter(Treatment == "QC"))

outlier_samples <- c(names(zK)[zK < -sample_connectivity_threshold],
		     names(zK)[zK > sample_connectivity_threshold])

# There are no QC sample outliers.
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
# Good: by removing outlier protein measurments, we now have one less
# sample level outlier.
message(paste0("\nRemoving outlier sample(s):\n",
	       paste(outlier_samples,collapse="\n")))

# Remove sample outliers.
filt_protein <- filt_protein %>% filter(Sample %notin% outlier_samples)

#---------------------------------------------------------------------
# Reformat data for TAMPOR normalization script.
#---------------------------------------------------------------------

# Reformat final normalized data for TAMPOR Normalization.
cleanDat <- reformat_TAMPOR(filt_protein,samples)

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. [output_name]_tidy_peptide.csv -- raw peptide data.
# 2. [output_name]_cleanDat.RData" -- data for TAMPOR.

## Save key results.
message("\nSaving data for downstream analysis...")

# 0. gene_map.RData   - gene identifier map.
myfile <- file.path(rdatdir,"gene_map.RData")
saveRDS(gene_map,myfile)

# 1. [output_name]_tidy_peptide.csv -- raw peptide data.
myfile <- file.path(rdatdir,paste(output_name,"tidy_peptide.csv",sep="_"))
fwrite(tidy_peptide,myfile)

# 2. [output_name]_cleanDat.RData -- data for TAMPOR.
myfile <- file.path(rdatdir,paste(output_name,"cleanDat.RData",sep="_"))
saveRDS(cleanDat,myfile)

message("\nDone!")
* Project '/mnt/d/projects/SynaptopathyProteomics' loaded. [renv 0.10.0]
Loading functions in SynaptopathyProteomics/R

Creating gene identifier map.

Loading raw data from Proteome Discover...

Performing sample loading normalization.

Imputing a small number of missing peptide values.
Imputing missing values in group: Shank2
... Percent missing values: 0.014% (n=71).
Imputing missing values in group: Shank3
... Percent missing values: 0.02% (n=103).
Imputing missing values in group: Syngap1
... Percent missing values: 0.039% (n=202).
Imputing missing values in group: Ube3a
... Percent missing values: 0.009% (n=48).

Removing peptides with irreproducible QC measurements.
Total number of QC outlier peptides identified: 611

Summarizing proteins as the sum of their peptides.

Performing sample loading normalization between experiments.

Performing ComBat to remove intra-batch batch effect.
Warning message:
For group: Shank2 - No detectable batch effect! Not applying ComBat. 
Initial coorelation between batch and samples: 0.828
Performing ComBat for group: Shank3
Final coorelation between batch and samples: 0.102 

Initial coorelation between batch and samples: 0.543
Performing ComBat for group: Syngap1
Final coorelation between batch and samples: 0.114 

Warning message:
For group: Ube3a - No detectable batch effect! Not applying ComBat. 

Standardizing protein measurements between experiments by IRS normalization.

Filtering proteins...
Removed proteins identified by one peptide.....: 95
Removed proteins with too many missing values..: 52
Removed proteins with any missing QC values....: 430
Final number of reproducibly quantified proteins: 2,860

Removing outlier sample(s):
Abundance: F4: 127C, Sample, KO, 41335, Shank3, Cortex
Final number of samples: 43

Saving data for downstream analysis...

Done!
#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

# Function to check command line input.
parse_args <- function(arg_opts = c("Cortex","Striatum"), nargs=1, err=NULL){
	args <- commandArgs(trailingOnly=TRUE)
	if (interactive()) { 
		set.seed(as.numeric(Sys.time()))
		n <- sample(c(1,2),1)
		warning(paste0("Running R interactively. ",
			       "Analyzing ",arg_opts[n]," data."))
		return(arg_opts[n])
	} else if (length(args) == nargs & args[nargs] %in% arg_opts) {
		return(args[nargs])
	} else {
		stop(err, call.=FALSE)
	}
}

## Parse input arguments.
msg <- "Please specify either 'Cortex' or 'Striatum' for the analysis."
analysis_type <- parse_args(err=msg)

## Inputs:
# Input data should be in root/data/:
# 1. TMT-samples.csv - sample meta data.
sample_files = list("Cortex" = "4227_TMT_Cortex_Combined_traits.csv",
		     "Striatum" = "4227_TMT_Striatum_Combined_traits.csv")
input_samples <- sample_files[[analysis_type]]

# 2. TMT-raw-peptide.csv - raw peptide data from PD.
data_files = list("Cortex" = "4227_TMT_Cortex_Combined_PD_Peptide_Intensity.csv",
		  "Striatum" = "4227_TMT_Striatum_Combined_PD_Peptide_Intensity.csv")
input_data <- data_files[[analysis_type]]

## Other parameters:
output_name = analysis_type # Prefix for naming output files.
sample_connectivity_threshold = 2.5 # Sample level outlier threshold. 

## Output for downstream analysis:
# Stored in root/rdata/
# 0. [output_name]_gene_map.RData     - gene identifier map.
# 1. [output_name]_tidy_peptide.csv   - tidy, raw peptide data.
# 2. [output_name]_preprocessed.RData - preprocessed data for TAMPOR.

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
# Alternatively, use: 
# renv::load("/mnt/d/projects/SynaptopathyProteomics") # Path to project root.

# Load required packages and functions.
suppressPackageStartupMessages({
	library(dplyr)
	library(TBmiscr)
	library(getPPIs)
	library(data.table)
})

# Load additional functions.
TBmiscr::load_all()

# Project directories:
rootdir <- getrd()
funcdir <- file.path(rootdir, "R")
datadir <- file.path(rootdir, "data")
figsdir <- file.path(rootdir, "figs")
rdatdir <- file.path(rootdir, "rdata")
tabsdir <- file.path(rootdir, "tables")
downdir <- file.path(rootdir, "downloads")

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

# Map any missing ids by hand.
not_mapped <- is.na(entrez)
if (sum(not_mapped) > 0) {
	message("Mapping missing IDs by hand.")
	mapped_by_hand <- c("P10853" = 319180)
	entrez[names(mapped_by_hand)] <- mapped_by_hand
}

# Check: no unmapped ids?
check <- sum(is.na(entrez)) == 0
if (!check) { stop( "Unable to map all uniprot IDs to stable entrez IDs.") }

# Check: we have successfully mapped all uniprot ids.
check <- sum(is.na(entrez)) == 0
if (!check) { stop("Unable to map all Uniprot IDs to stable gene identifiers!") }

# Map entrez ids to gene symbols.
symbols <- getPPIs::getIDs(entrez,from="entrez",to="symbol",species="mouse")

# Check: we have successfully mapped all gene symbols.
check <- sum(is.na(symbols)) == 0
if (!check) { stop("Unable to map all Entrez IDs to gene Symbols!") }

# Create mapping data.table.
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
# This strategy was adapted from Ping et al., 2018 (pmid:29533394).
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
## Intra-batch Protein-lavel ComBat.
#---------------------------------------------------------------------

# Before tackling inter-experimental batch effects, address the 
# intra-experimental batch effect.

# Each experimental cohort of 8 was prepared in two batches. 
# This was necessary because the ultra-centrifuge rotor used to 
# prepare purified synaptosomes holds a maximum of 6 samples. 
# This intra-batch batch effect was recorded for 6/8 batches. 
# Here I will utilize the sva::ComBat() function in order to account 
# for this batch effect before correcting for inter-batch batch
# effects with IRS normalization. 

# NOTE: In the absence of evidence of a batch effect 
# (not annotated or cor(PCA,batch)<0.1), then ComBat is not applied.
# If ComBat is not applied, then the data is returned un-regressed.

# NOTE: The values of QC samples are not adjusted by ComBat.
# QC samples were prepared from a seperate batch of mice and
# represent a single batch.

message("\nPerforming ComBat to remove intra-batch batch effect.")

# Perform ComBat for each dataset.
data_combat <- intrabatch_combat(sl_protein, samples,
				 group="Shank2", batch="PrepDate")
data_combat <- intrabatch_combat(data_combat, samples,
				 group="Shank3", batch="PrepDate")
data_combat <- intrabatch_combat(data_combat, samples,
				 group="Syngap1", batch="PrepDate")
data_combat <- intrabatch_combat(data_combat, samples,
				 group="Ube3a", batch="PrepDate")
norm_protein <- data_combat

#---------------------------------------------------------------------
## Insure there are no QC outlier samples.
#---------------------------------------------------------------------

# Remove QC outliers before performing IRS normalization.

# Calculate Oldham's normalized sample connectivity (zK) in order to 
# identify outlier samples.
# This approach was adapted from Oldham et al., 2012 (pmid: 22691535).
zK <- sampleConnectivity(norm_protein %>% filter(Treatment == "QC"))

outlier_samples <- c(names(zK)[zK < -sample_connectivity_threshold],
		     names(zK)[zK > sample_connectivity_threshold])

# There are no QC sample outliers.
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
irs_protein <- normIRS(norm_protein,controls="QC",robust=TRUE)

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

# Remove sample outliers.
filt_protein <- filt_protein %>% filter(Sample %notin% outlier_samples)

# Status:
# Good: by removing outlier protein measurments, we now have one less
# sample level outlier.
message(paste0("\nRemoving outlier sample(s):\n",
	       paste(outlier_samples,collapse="\n")))
message(paste("Final number of samples:",
	      length(unique(filt_protein$Sample))))

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. [output_name]_tidy_peptide.csv -- raw peptide data.
# 2. [output_name]_cleanDat.RData" -- data for TAMPOR.

## Save key results.
message("\nSaving data for downstream analysis...")

# 0. gene_map.RData   - gene identifier map.
myfile <- file.path(rdatdir,paste(output_name,"gene_map.RData",sep="_"))
saveRDS(gene_map,myfile)

# 1. [output_name]_tidy_peptide.csv -- raw peptide data.
myfile <- file.path(rdatdir,paste(output_name,"tidy_peptide.csv",sep="_"))
fwrite(tidy_peptide,myfile)

# 2. [output_name]_preprocessed.RData -- data for TAMPOR.
myfile <- file.path(rdatdir,paste(output_name,"preprocessed.RData",sep="_"))
saveRDS(filt_protein,myfile)

message("\nDone!")
* Project '/mnt/d/projects/SynaptopathyProteomics' loaded. [renv 0.10.0]
Loading functions in SynaptopathyProteomics/R

Creating gene identifier map.
Warning message:
In getPPIs::mgi_batch_query(uniprot, quiet = FALSE, download = FALSE) :
  Unable to map 1 gene identifiers to entrez.
Mapping missing IDs by hand.

Loading raw data from Proteome Discover...

Performing sample loading normalization.

Imputing a small number of missing peptide values.
Imputing missing values in group: Shank2
... Percent missing values: 0.029% (n=148).
Imputing missing values in group: Shank3
... Percent missing values: 0.014% (n=73).
Imputing missing values in group: Syngap1
... Percent missing values: 0.017% (n=86).
Imputing missing values in group: Ube3a
... Percent missing values: 0.02% (n=104).

Removing peptides with irreproducible QC measurements.
Total number of QC outlier peptides identified: 423

Summarizing proteins as the sum of their peptides.

Performing sample loading normalization between experiments.

Performing ComBat to remove intra-batch batch effect.
Initial coorelation between batch and samples: 0.588
Performing ComBat for group: Shank2
Final coorelation between batch and samples: 0.04 

Warning message:
For group: Shank3 - Cannot perform ComBat with one batch! 
Warning message:
For group: Syngap1 - Cannot perform ComBat with one batch! 
Warning message:
For group: Ube3a - No detectable batch effect! Not applying ComBat. 

Standardizing protein measurements between experiments by IRS normalization.

Filtering proteins...
Removed proteins identified by one peptide.....: 124
Removed proteins with too many missing values..: 51
Removed proteins with any missing QC values....: 605
Final number of reproducibly quantified proteins: 2,858

Removing outlier sample(s):
Abundance: F2: 128C, Sample, KO, 38894, Ube3a, Striatum
Abundance: F2: 131N, Sample, KO, 38898, Ube3a, Striatum
Final number of samples: 42

Saving data for downstream analysis...

Done!
message(paste0("Summary of differentially abundant proteins ",
	      "in each subceulluar fraction (FDR < ",alpha,"):"))
results$summary

#---------------------------------------------------------------------
## Examine PCA of final normalized data.
#---------------------------------------------------------------------

# Generate the plot.
# FIXME: customize colors!
plot <- ggplotPCA(tidy_protein)

# Save
myfile <- file.path(figsdir,"PCA.pdf")
ggsave(myfile,plot,height=5,width=5)

#---------------------------------------------------------------------
## Create protein co-expression matrix.
#---------------------------------------------------------------------

# Create data matrix with final normalized data.
norm_protein <- tidy_protein %>%
	dcast(Accession ~ Sample,value.var="Intensity") %>% 
	as.matrix(rownames=TRUE) %>% log2

# Remove QC samples when calculating co-expression matrices.
samples_to_ignore <- "SPQC"
idy <- grepl(samples_to_ignore,colnames(norm_protein))
subdat <- norm_protein[,!idy]

# Get list of data: all, just wt, just mutant.
samples_to_keep <- c(All="",WT="Control",KO="Mutant")
data_list <- sapply(samples_to_keep, function(x) { 
			    subdat[,grep(x,colnames(subdat))]
		 })

# Calculate bicor correlation matrix for all subsets of the data.
# Silence suppress unwanted output generated by WGNA::bicor().
message("\nCreating protein co-expression network.")
silence({ adjm_list <- lapply(data_list,function(x) bicor(t(x))) })

# Perform network enhancment for all matrices.
message(paste("\nPerforming network enhancement to denoise network."))
ne_list <- lapply(adjm_list, function(x) neten::neten(x))

#---------------------------------------------------------------------
## Save the TMT data and stats as a single excel workbook
#---------------------------------------------------------------------

## Save TMT data as a single excel document...
## Create an excel workbook with the following sheets:
# * Samples
# * Input - raw peptides
# * Output - normalized protein
# * Statistical results:
#     - Seperate sheet for each fraction/comparison.
#     - Include contrast specific data.

# Sort the statistical results.
glm_results <- glm_results[c("F4","F5","F6","F7","F8","F9","F10")]

# Loop to add normalized protein data to glm statistical results.
message("\nSaving TMT data and statistical results.")
for (i in 1:length(glm_results)){
	df <- glm_results[[i]]
	namen <- names(glm_results)[i]
	subsamples <- samples$Sample[grepl(namen,samples$Fraction)]
	idy <- match(subsamples,colnames(norm_protein))
	dm <- norm_protein[,idy]
	# Sort the data by Exp.Sample
	ids <- sapply(sapply(strsplit(colnames(dm),"\\."),"[",c(6,7),
			     simplify=FALSE), paste,collapse=".")
	names(ids) <- colnames(dm)
	ids <- ids[order(ids)]
	dm <- dm[,names(ids)]
	# Coerce dm to dt and merge with stats df.
	dt <- as.data.table(dm,keep.rownames="Accession")
	dt_out <- left_join(df,dt,by="Accession") %>% as.data.table
	glm_results[[i]] <- dt_out
} # Ends loop.

# Map uniprot to entrez IDs.
uniprot <- glm_results[[1]][["Accession"]]
idx <- match(uniprot,gene_map$uniprot)
entrez <- gene_map$entrez[idx]
names(entrez) <- uniprot

# Map missing ids by hand.
is_missing <- is.na(entrez)
missing_entrez <- entrez[is_missing]
mapped_by_hand <- c("P00848"=17705,
		    "P62806"=326619,
		    "Q64525"=319189,
		    "P06330"=111507, 
		    "Q5PR69"=320827, 
		    "Q80WG5"=241296,
		    "P03975"=15598, 
		    "Q80TK0"=100503185)
entrez[names(mapped_by_hand)] <- mapped_by_hand

# Check that we have successfully mapped all uniprot ids to entrez.
check <- sum(is.na(entrez))
if (check != 0) { stop("There are unmapped uniprot IDs.") }

# Map entrez to gene symbols.
idx <- match(entrez,gene_map$entrez)
symbols <- gene_map$symbol[idx]
names(symbols) <- entrez

# Loop to add Entrez IDs and Gene Symbols to data.
# Also fix column titles at the same time. 
for (i in 1:length(glm_results)){
	df <- glm_results[[i]]
	colnames(df) <- gsub("\\."," ",colnames(df))
	idy <- grepl("Abundance",colnames(df))
	colnames(df)[idy] <- paste("Log2",colnames(df)[idy])
	df <- tibble::add_column(df,"Entrez"=entrez[df$Accession],
				 .after="Accession")
	df <- tibble::add_column(df,"Gene"=symbols[as.character(df$Entrez)],
				 .after="Entrez")
	glm_results[[i]] <- df
}

## Save as excel workboook.
names(glm_results) <- paste(names(glm_results),"Results")
results <- c(list("Samples" = samples,
		"Normalized Protein" = tidy_protein),
	     glm_results)
myfile <- file.path(tabsdir,"Swip_TMT_Results.xlsx")
write_excel(results,myfile,rowNames=FALSE)

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:
# Stored in root/rdata/
# 0. gene_map.RData   - gene identifier map.
# 1. tidy_peptide.csv - raw peptide data.
# 2. tidy_protein.csv - final normalized protein data.
# 3. norm_protein.csv - final normalized protein data matrix.
# 4. adjm.csv         - adjacency matrix. 
# 5. ne_adjm.csv      - enhanced adjacency matrix. 
# 6. glm_stats.RData  - statistical results from edgeR glm

## Save key results.
message("\nSaving data for downstream analysis")

# 0. Save gene map.
myfile <- file.path(rdatdir,"gene_map.RData")
saveRDS(gene_map,myfile)

# 1. Save tidy_peptide.
myfile <- file.path(rdatdir,"tidy_peptide.csv")
fwrite(tidy_peptide,myfile)

# 2. Save tidy_protein.
myfile <- file.path(rdatdir,"tidy_protein.csv")
fwrite(tidy_protein,myfile)

# 3. Save normalized protein data matrix.
myfile <- file.path(rdatdir,"norm_protein.csv")
norm_protein %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(myfile)

# Save WT subset of data matrix.
myfile <- file.path(rdatdir,"norm_wt_protein.csv")
data_list[["WT"]] %>% as.data.table(keep.rownames="Accession") %>%
	       fwrite(myfile)
myfile <- file.path(rdatdir,"norm_ko_protein.csv")

# Save KO subset of data matrix.
data_list[["KO"]] %>% as.data.table(keep.rownames="Accession") %>%
	       fwrite(myfile)

# 4. Save adjms.
myfile <- file.path(rdatdir,"adjm.csv")
adjm_list[["All"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

myfile <- file.path(rdatdir,"wt_adjm.csv")
adjm_list[["WT"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

myfile <- file.path(rdatdir,"ko_adjm.csv")
adjm_list[["KO"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

# 5. Save enhanced adjms.
myfile <- file.path(rdatdir,"ne_adjm.csv")
ne_list[["All"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

myfile <- file.path(rdatdir,"wt_ne_adjm.csv")
ne_list[["WT"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

myfile <- file.path(rdatdir,"ko_ne_adjm.csv")
ne_list[["KO"]] %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)
