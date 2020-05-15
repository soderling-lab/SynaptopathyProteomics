#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Parse input arguments.
msg <- "Please specify either 'Cortex' or 'Striatum' for the analysis."
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) { stop(msg) } else { analysis_type <- args[1] }

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
all_combat <- function(sl_protein,samples) {
	data_combat <- intrabatch_combat(sl_protein, samples, group="Shank2", batch="PrepDate")
	data_combat <- intrabatch_combat(data_combat, samples, group="Shank3", batch="PrepDate")
	data_combat <- intrabatch_combat(data_combat, samples, group="Syngap1", batch="PrepDate")
	data_combat <- intrabatch_combat(data_combat, samples, group="Ube3a", batch="PrepDate")
	return(data_combat)
}

norm_protein <- all_combat(sl_protein,samples)

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
# NOTE: Don't remove proteins with outlier measurements.

symbols_ignore = c("Shank2","Shank3","Syngap1","Ube3a")
uniprot_ignore <- gene_map$uniprot[match(symbols_ignore,gene_map$symbol)]
names(uniprot_ignore) <- symbols_ignore

message(paste("\nFiltering proteins..."))
filt_protein <- filtProt(irs_protein,
			 controls="QC",
			 remove.protein.outliers=TRUE,
			 ignore = uniprot_ignore,
			 nbins=5,nSD=4,
			 summary=TRUE)

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

quit()

#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
## Final Normalization: Regression of covariates.
#---------------------------------------------------------------------

# Batch and sex are used as covariates.
eblm_protein <- eBLM_regression(filt_protein, 
				traits=samples, ignore="QC")

# Less sigprots after regression....
# More sigprots after TAMPOR...

#---------------------------------------------------------------------
## Protein differential abundance.
#---------------------------------------------------------------------

# Use EdgeR glm to evaluate differential protein abundance.
# Compare WT v Mutant within a fraction.
message(paste("\nEvaluating differential abundance between WT and",
	      "Swip mutant samples."))

alpha <- 0.1
comparisons <- "Genotype.Treatment"
data_in <- filt_protein
results <- glmDA(data_in,comparisons,samples,
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

write_excel(glm_results,"foo.xlsx")

#--------------------------------------------------------------------
## Identify subset of highly reproducible proteins
#--------------------------------------------------------------------

sub_samples <- samples %>% filter(Sample %in% filt_protein$Sample)

data_in <- left_join(filt_protein,sub_samples,by=c("Sample","Experiment","Channel","Treatment"))

# Split the data into a list of proteins.
prot_list <- data_in %>% group_by(Accession) %>% group_split()
names(prot_list) <- sapply(prot_list,function(x) unique(x$Accession))

# Define a function that checks reproducibility of a protein.
check_reproducibility <- function(prot_list,protein,genotype,
				  fun="bicor",threshold=0.6) {
	# Which proteins are highly reproducible between replicates.
	# For a given protein, cast the data into a matrix: Fraction ~ Replicate.
	dt <- prot_list[[protein]] %>% 
		filter(Treatment == genotype) %>%
		as.data.table() 
	# Treat all mutants the same:
	dt$Treatment <- gsub("HET","KO",dt$Treatment) 
	dt$Replicate <- as.numeric(as.factor(dt$Channel))
        dt_temp <- dcast(dt, Genotype + Tissue ~ Treatment + Replicate, 
		      value.var="Intensity") 
	value_cols <- grep("_",colnames(dt_temp))
	dm <- dt_temp %>% dplyr::select(all_of(value_cols)) %>% 
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
				       fun="bicor",
				       threshold=0.8)
}

checks <- unlist(check)
sum(checks==6)

message(paste("Numer of reproducible proteins:",sum(checks)))
