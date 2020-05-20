#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## If running interactively, define the tissue type to be analyzed:
if (interactive()) {
	analysis_type = "Cortex"
## Otherwise, parse command line arguments:
} else if (!interactive()) {
	args <- commandArgs(trailingOnly=TRUE)
	if (length(args) == 1) { 
		# Check that only 1 arg is passed.
		analysis_type = commandArgs(trailingOnly=TRUE)[1]
	} else { 
		stop("Specify either 'Cortex' or 'Striatum'.",call.=FALSE) 
	}
}

## Optional parameters:
alpha_threshold = 0.1
sample_connectivity_threshold = 2.5
n_threads = parallel::detectCores()

## Processing switches
group.WT = FALSE
scale.WT = FALSE
rand.batch = FALSE
remove.protein.outliers = FALSE

## Input data in root/data/:
# 1. TMT-samples.csv - sample meta data.
input_samples = list(Cortex = "Cortex_Samples.csv",
		     Striatum = "Striatum_Samples.csv")[[analysis_type]]

# 2. TMT-raw-peptide.csv - raw peptide data from PD.
input_data = list(Cortex = "Cortex_Peptides.csv",
		  Striatum = "Striatum_Peptides.csv")[[analysis_type]]

## Output for downstream analysis, stored in root/rdata/
if (group.WT) {
	output_name = paste(analysis_type,"grouped_WT",sep="_")
} else {
	output_name = analysis_type # Prefix for naming output files.
}

if (scale.WT) {
	output_name = paste(output_name,"scaled_WT",sep="_")
} else {
	# pass
	output_name = output_name
}

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
tabsdir <- file.path(root, "tables")

#---------------------------------------------------------------------
## Load the raw data and sample info.
#---------------------------------------------------------------------

message(paste0("\nAnalyzing ",analysis_type,"..."))

# Load the TMT data.
myfile <- file.path(datadir,input_data)
peptides <- fread(myfile)

# Load sample information.
myfile <- file.path(datadir,input_samples)
samples <- fread(myfile)

#---------------------------------------------------------------------
## Map all Uniprot accession numbers to stable entrez IDs.
#---------------------------------------------------------------------
# Create a gene identifier map:
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

# Collect all uniprot accession IDs.
all_accession <- unique(peptides$Accession)

# Map Uniprot IDs to entrez using MGI database.
# This can take a couple minutes if download == TRUE, as 
# the function currently downloads the MGI data each time it 
# is function is called.
entrez <- getPPIs::mgi_batch_query(all_accession,quiet=FALSE,download=FALSE)
names(entrez) <- all_accession

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
gene_map <- data.table(Accession = names(entrez),
                       Entrez = entrez,
	               Symbol = symbols)
gene_map$ID <- paste(gene_map$Symbol,gene_map$Accession,sep="|")

#---------------------------------------------------------------------
## Tidy-up the input data from Proteome Discover.
#---------------------------------------------------------------------
# Convert PD data.frame into tidy data.table.
# NOTE: Samples should contain the following columns:
# Treatment, Channel, Sample, Genotype

message("\nLoading raw peptide data from Proteome Discover.")

cols <- colnames(peptides)[!grepl("Abundance",colnames(peptides))]
tidy_peptide <- tidyProt(peptides,intensity.cols=cols)

# Annotate tidy data with additional meta data from samples.
tidy_peptide <- left_join(tidy_peptide,samples,by="Sample")

#---------------------------------------------------------------------
## Perform sample loading normalization.
#---------------------------------------------------------------------
# Perform sample loading normalization. Normalization is done for each 
# experiment independently (group by Genotype:Sample).

message("\nPerforming sample loading normalization.")

sl_peptide <- normSL(tidy_peptide, groupBy=c("Genotype","Sample"))

#---------------------------------------------------------------------
## Impute missing peptide values.
#---------------------------------------------------------------------
# Impute missing peptide values with the KNN algorithm for data that
# is not missing at random (MNAR). Missing values are inferred to be
# missing because they are low abundant species. This can be confirmed
# by plotting the density distribution of peptides with and without 
# missing values.
# * NOTE: Missing QC values will not be imputed.
# * NOTE: Peptides (rows) with more than 50% missingness will not be imputed.
# Values in these rows are masked (replaced with NA). 

message("\nImputing a small number of missing peptide values...")

imputed_peptide <- imputeKNNpep(sl_peptide, groupBy="Genotype",
				ignore.samples="QC",quiet=FALSE) 

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

message("\nRemoving peptides with irreproducible QC measurements...")

filt_peptide <- filtQC(imputed_peptide, grouping.col='Treatment',
		       controls='QC',quiet=FALSE)

#---------------------------------------------------------------------
## Summarize to protein level.
#---------------------------------------------------------------------
# Protein intensity is just sum of all peptides for that protein.

# Sum to protein level.
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

# NOTE: The values of QC samples are not adjusted by ComBat
# because QC samples were prepared as a seperate batch and
# represent a single batch.

if (rand.batch){
	idx <- grepl("Syngap1",samples$Sample) & !grepl("QC",samples$Treatment)
	set.seed(as.numeric(Sys.time()))
	samples$PrepDate[idx] <- sample(c(rep(0,4),rep(1,4)))
	x = samples$PrepDate[idx]
	names(x) <- samples$Sample[idx]
	print(x)
}

# Perform ComBat for each dataset.
message("\nPerforming ComBat to remove intra-batch batch effect...")
all_combat <- function(sl_protein,samples) {
	# intrabatch_combat returns the input data, adjusted for batch effect.
	data_combat <- intrabatch_combat(sl_protein, samples, 
					 group="Shank2", batch="PrepDate")
	data_combat <- intrabatch_combat(data_combat, samples, 
					 group="Shank3", batch="PrepDate")
	data_combat <- intrabatch_combat(data_combat, samples, 
					 group="Syngap1", batch="PrepDate")
	data_combat <- intrabatch_combat(data_combat, samples, 
					 group="Ube3a", batch="PrepDate")
	return(data_combat)
}
combat_protein <- all_combat(sl_protein,samples)

#---------------------------------------------------------------------
## Insure there are no QC outlier samples.
#---------------------------------------------------------------------
# Remove QC outliers before performing IRS normalization.
# This approach was adapted from Oldham et al., 2012 (pmid: 22691535).

# Calculate Oldham's normalized sample connectivity (zK) in order to 
# identify outlier samples.
sample_connectivity_threshold = 2.5
zK <- sampleConnectivity(combat_protein %>% filter(Treatment == "QC"))
qc_outliers <- c(names(zK)[zK < -sample_connectivity_threshold],
		     names(zK)[zK > sample_connectivity_threshold])

# Check, there should be no QC sample outliers.
check <- length(qc_outliers) == 0
if (!check) { stop("Why are there outlier samples?") }

#---------------------------------------------------------------------
## Perform IRS Normalization.
#---------------------------------------------------------------------
# Equalize QC measurements between experiments. IRS normalization 
# accounts for protein quantification by different peptides in 
# each experiment.

message(paste("\nStandardizing protein measurements between",
	      "experiments by IRS normalization."))
irs_protein <- normIRS(combat_protein,controls="QC",robust=TRUE)

#---------------------------------------------------------------------
## Protein level filtering.
#---------------------------------------------------------------------
# Remove proteins that are:
# * Identified by a single peptide.
# * Have too many missing values.
# * Have any missing QC values.
# * Have outlier outlier measurments.
message(paste("\nFiltering proteins..."))

# Ignore proteins encoded by manipulated genes.
symbols_ignore = c("Shank2","Shank3","Syngap1","Ube3a")
uniprot_ignore <- gene_map$Accession[match(symbols_ignore,gene_map$Symbol)]
names(uniprot_ignore) <- symbols_ignore

# Filter proteins.
filt_protein <- filtProt(irs_protein,
			 controls="QC",
			 remove.protein.outliers=remove.protein.outliers,
			 ignore = uniprot_ignore,
			 nbins=5,nSD=4,
			 summary=TRUE)

# There should be no missing values at this point.
check <- sum(is.na(filt_protein$Intensity)) == 0
if (!check) { stop("Why are there missing values!") }

#---------------------------------------------------------------------
## Scale WT's to be ~equal.
#---------------------------------------------------------------------

# Scale WT?
if (scale.WT) {
	scale_protein <- normIRS(filt_protein,controls="WT",robust=TRUE)
} else {
	scale_protein <- filt_protein
}

#---------------------------------------------------------------------
## Check for sample outliers -- remove if found.
#---------------------------------------------------------------------

# Calculate Oldham's normalized sample connectivity (zK) in order to 
# identify outlier samples.
# This approach was adapted from Oldham et al., 2012 (pmid: 22691535).

n_samples <- length(unique(scale_protein$Sample))
message(paste("Initial number of samples:",n_samples))

# Loop to iteratively remove outliers.
i = 0
outlier_samples <- NULL
while (i < 5) {
	zk <- sampleConnectivity({
		scale_protein %>% 
			filter(Treatment != "QC") %>%
			filter(Sample %notin% outlier_samples)
	})
	outlier_samples <- unique({
		c(outlier_samples,
		  names(zk)[zk < -sample_connectivity_threshold],
		  names(zk)[zk > sample_connectivity_threshold])
	})
	i = i + 1
}

# Remove sample outliers.
final_protein <- scale_protein %>% filter(Sample %notin% outlier_samples)
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
## Annotate final protein data with sample information.
#---------------------------------------------------------------------

# Annotate protein data with sample meta data.
# Shared column names:
cols <- intersect(colnames(final_protein),colnames(samples))
final_protein <- left_join(final_protein,samples,by=cols) %>% 
	as.data.table()

# Add entrez ids and gene symbols to data.
idx <- match(final_protein$Accession,gene_map$Accession)
Symbol <- gene_map$Symbol[idx]
Entrez <- gene_map$Entrez[idx]
final_protein <- tibble::add_column(final_protein,Symbol,.after="Accession")
final_protein <- tibble::add_column(final_protein,Entrez,.after="Symbol")

#---------------------------------------------------------------------
## DO TAMPOR HERE.
#---------------------------------------------------------------------

# Define experimental batches:
batch <- interaction(final_protein$Tissue,final_protein$Genotype) %>% 
	as.numeric()
batch <- paste0("b",batch)

# Annotate data with batch and batch.channel (ID).
final_protein$batch <- batch
final_protein$ID <- interaction(batch,final_protein$Channel) %>% 
	as.character()

# Cast into a data.matrix for TAMPOR.
dm <- final_protein %>% dcast(Accession ~ ID, value.var="Intensity") %>% 
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
traits <- final_protein %>% dplyr::select(batch,ID) %>% unique()
traits <- as.data.frame(traits)
rownames(traits) <- traits[,"ID"] # rownames must be batch.channel.

# Define controls or global internal standards (GIS) for 
# TAMPOR normalization. Use all WT samples. 
controls <- unique(filter(final_protein,Treatment=="WT")[["ID"]])

# Define samples to be ignored -- QC samples.
qc_samples <- unique(filter(final_protein,Treatment=="QC")[["ID"]])

# Perform the normalization.
data_tampor <- TAMPOR(dat = dm, traits = traits, batchPrefixInSampleNames = TRUE,
		      samplesToIgnore = qc_samples, GISchannels = controls, parallelThreads = n_threads)

# Collect normalized, relative abundance data post-TAMPOR.
tampor_protein <- data_tampor$cleanRelAbun %>% 
	reshape2::melt(value.name="Abundance")
colnames(tampor_protein)[c(1,2)] <- c("Accession","ID")

# Clean-up and merge with input data.
tampor_protein$ID <- as.character(tampor_protein$ID)
tampor_protein$Accession <- as.character(tampor_protein$Accession)
norm_protein <- left_join(final_protein,tampor_protein,by=c("Accession","ID"))
norm_protein$batch <- NULL
norm_protein$ID <- NULL

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

group.WT = TRUE

if (group.WT) {
	# Perform analysis without WT grouping. All contrasts are performed
	# within a genotype.
	comparisons = c("Genotype","Treatment")
	glm_results <- glmDA(norm_protein, treatment.ignore = "QC", 
			     comparisons, value.var="Abundance",
			     combine.WT = TRUE)
} else {
	# Perform analysis WITH WT grouping.
	# Assess changes in protein abundance using a glm to account for
	# differences in genetic background (genotype).
	comparisons = c("Genotype","Treatment")
	glm_results <- glmDA(norm_protein, treatment.ignore = "QC", 
			     comparisons, value.var="Abundance",
			     combine.WT = FALSE)
}

#  Annotate glm_results with gene identifiers.
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

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------

## Output for downstream analysis:

## Save key results.
message("\nSaving data for downstream analysis.")

# [output_name]_tidy_peptide.csv -- raw peptide data.
myfile <- file.path(rdatdir,paste(output_name,"final_protein.csv",sep="_"))
fwrite(final_protein,myfile)

# [output_name]_GLM_Results.xlsx
myfile <- file.path(tabsdir,paste(output_name,"GLM_Results.xlsx",sep="_"))
write_excel(glm_results,myfile)

message("\nDone!")
