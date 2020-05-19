#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
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

## Other optional parameters:
alpha_threshold = 0.1 # FDR significance threshold.
sample_connectivity_threshold = 2.5 # Sample level outlier threshold. 

## Input data in root/data/:
# 1. TMT-samples.csv - sample meta data.
input_samples = list(Cortex = "Cortex_Samples.csv",
		     Striatum = "Striatum_Samples.csv")[[analysis_type]]

# 2. TMT-raw-peptide.csv - raw peptide data from PD.
input_data = list(Cortex = "Cortex_Peptides.csv",
		  Striatum = "Striatum_Peptides.csv")[[analysis_type]]

## Output for downstream analysis, stored in root/rdata/
output_name = analysis_type # Prefix for naming output files.
# * [output_name]_tidy_peptide.csv  - tidy, raw peptide data.
# * [output_name]_norm_protein      - final preprocessed data.
# * [output_name]_glm_stats.csv     - glm statistical results.
# * [output_name]_glm_protein.csv   - glm fitted protein values.
# * [output_name]_GLM_Results.xlsx  - glm statistical results.
# * [output_name]_repro_prots.RData - potential marker proteins.

## Final tables, stored in root/tables/
# * [output_name]_GLM_Results.xlsx - glm results results.

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

# Perform ComBat for each dataset.
message("\nPerforming ComBat to remove intra-batch batch effect...")
all_combat <- function(sl_protein,samples) {
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
zK <- sampleConnectivity(combat_protein %>% filter(Treatment == "QC"))

qc_outliers <- c(names(zK)[zK < -sample_connectivity_threshold],
		     names(zK)[zK > sample_connectivity_threshold])

# There are no QC sample outliers.
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
uniprot_ignore <- gene_map$uniprot[match(symbols_ignore,gene_map$symbol)]
names(uniprot_ignore) <- symbols_ignore

# Filter proteins.
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
zK <- sampleConnectivity(filt_protein %>% filter(Treatment != "QC"))
if (length(zK) != 32) { stop("Whoops, something went wrong.") }

# Outliers:
outlier_samples <- c(names(zK)[zK < -sample_connectivity_threshold],
		     names(zK)[zK > sample_connectivity_threshold])

# Remove sample outliers.
norm_protein <- filt_protein %>% filter(Sample %notin% outlier_samples)

# Annotate protein data with sample meta data.
# Shared column names:
cols <- intersect(colnames(norm_protein),colnames(samples))
norm_protein <- left_join(norm_protein,samples,by=cols) %>% 
	as.data.table()

# Add entrez ids and gene symbols to data.
idx <- match(norm_protein$Accession,gene_map$uniprot)
Symbol <- gene_map$symbol[idx]
Entrez <- gene_map$entrez[idx]
norm_protein <- tibble::add_column(norm_protein,Symbol,.after="Accession")
norm_protein <- tibble::add_column(norm_protein,Entrez,.after="Symbol")

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

# Drop QC samples before EdgeR analysis.
data_in <- dt_eblm %>% filter(Treatment != "QC")

# Asses changes in protein abundance using a glm to account for
# differences in genetic background (genotype).
data_glm <- glmDA(data_in,comparisons=c("Genotype","Treatment"),
		  model="~background + group", samples, gene_map)

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
