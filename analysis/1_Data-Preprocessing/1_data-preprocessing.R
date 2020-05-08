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

## Parameters:
sample_connectivity_threshold = 2.5 # Threshold for sample level outliers
alpha = 0.1 # FDR threshold for differential abundance

## Main Outputs:
# Stored in root/tables/
# 0. Swip_TMT_Results.xlsx

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
# * Examine QC peptide reproducibility -- remove outliers.
# * Summarize to protein level by summing all peptides for a protein.
# * Protein-level sample loading normalization -- done across all experiments.
# * IRS normalization -- equalizes protein measurements made from different
#   peptides.
# * Sample pool normalization -- assumes that mean of QC technical replicates is
#   equal to mean of biological replicates..
# * Impute missing protein values with KNN (k=10).
# * Protein level filtering -- remove proteins identified by a single peptide;
#   remove proteins with too many missing values; remove proteins that are not 
#   reproducible (exhibit high inter-experimental variablility across the 3x 
#   biological replicates.
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
# downlaods the MGI data each time the function is called.
entrez <- mgi_batch_query(uniprot,quiet=FALSE) # Warning  not mapped.
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
cols=colnames(peptides)[!grepl("Abundance",colnames(peptides))]
tidy_peptide <- tidyProt(peptides,intensity.cols=cols)

# Annotate tidy data with additional meta data from samples.
tidy_peptide <- left_join(tidy_peptide,samples,by="Sample")

#---------------------------------------------------------------------
## Perform sample loading normalization.
#---------------------------------------------------------------------

# Perform sample normalization. Normalization is down for each 
# experiment independently (group by Experiment:Sample).
message("\nPerforming sample loading normalization.")
sl_peptide <- normSL(tidy_peptide, groupBy=c("Model","Sample"))

#---------------------------------------------------------------------
## Impute missing peptide values.
#---------------------------------------------------------------------

# Impute missing peptide values with KNN algorithm for MNAR data.
# * Missing QC values will not be imputed.
# * Peptides (rows) with more than 50% missingness will not be imputed.
# Values in these rows are masked (replaced with NA).
message("\nImputing a small number of missing peptide values.")
imputed_peptide <- imputeKNNpep(sl_peptide, groupBy="Model",
				samples_to_ignore="QC") 

#---------------------------------------------------------------------
## Examine reproducibility of QC measurements.
#---------------------------------------------------------------------

# Examine reproducibility of QC measurements.
# This strategy was adapted from Ping et al., 2018 (pmid: 29533394).
# For each experiment, calculate the ratio of QC measurements.
# Bin these ratios based on the average Intensity of QC peptides into
# 5 bins. For each bin, remove measuremnts that are outside 
# +/- 4x SD from the mean of all ratios within that bin.
message("\nRemoving peptides with irreproducible QC measurements.")
filt_peptide <- filtQC(imputed_peptide,controls="SPQC")

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
## Perform IRS Normalization.
#---------------------------------------------------------------------

# Equalize QC measurements between experiments. Adjusts protein 
# measurements of biological replicates simultaneously.
# Accounts for protein quantification by different peptides in 
# each experiment.
message(paste("\nStandardizing protein measurements between",
	      "experiments by IRS normalization."))
irs_protein <- normIRS(sl_protein,controls="SPQC",robust=TRUE)

#---------------------------------------------------------------------
## Perform Sample Pool Normalization.
#---------------------------------------------------------------------

# QC samples were generated by pooling all biological replicates.
# Normalize mean of all QC samples to be equal to the mean of all 
# biological replicates.
message(paste("\nAccounting for experimental batch effect by",
	      "performing sample pool normalization."))
spn_protein <- normSP(irs_protein,pool=c("Control","Mutant"))

#---------------------------------------------------------------------
## Protein level imputing.
#---------------------------------------------------------------------

# Impute missing values.
# Proteins with more than 50% missing values are ignored.
message(paste("\nImputing missing protein values."))
imputed_protein <- imputeKNNprot(spn_protein,ignore="SPQC")

#---------------------------------------------------------------------
## Protein level filtering.
#---------------------------------------------------------------------

# Remove proteins:
# Remove proteins that were identified by a single peptide.
# Remove proteins with too many missing values.
# Remove proteins with any missing QC values.
# Remove proteins that have outlier measurements.
message(paste("\nFiltering proteins..."))
filt_protein <- filtProt(imputed_protein,
			 controls="SPQC",nbins=5,nSD=4,summary=TRUE)

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

# There are no sample outliers.
check <- length(outlier_samples) == 0
if (!check) { stop("Why are there outlier samples?") }

#---------------------------------------------------------------------
## Protein differential abundance.
#---------------------------------------------------------------------

# Use EdgeR glm to evaluate differential protein abundance.
# Compare WT v Mutant within a fraction.
message(paste("\nEvaluating differential abundance between WT and",
	      "Swip mutant samples."))
comparisons <- "Genotype.Fraction"
results <- glmDA(filt_protein,comparisons,samples,
		 samples_to_ignore="SPQC",summary=TRUE)

# Collect the results.
tidy_protein <- results$data
glm_results <- results$results

# Summary of DA proteins:
message(paste0("Summary of differentially abundant proteins ",
	      "in each subceulluar fraction (FDR < ",alpha,"):"))
results$summary

# Proteins that are commonly dysregulated:
combined_results <- rbindlist(glm_results,idcol="Fraction")
idx <- match(combined_results$Accession,gene_map$uniprot)
combined_results$Gene <- gene_map$symbol[idx]
df <- combined_results %>% group_by(Accession) %>% 
	summarize(Gene = unique(Gene),
		  nSig = sum(FDR<0.1),
		  nFractions = length(FDR)) %>% 
	filter(nSig==nFractions)
message(paste("Proteins that are differentially abundant in all fractions:\n",
	      paste(df$Gene,collapse=", ")))

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

# 6. Save statistical results.
myfile <- file.path(rdatdir,"glm_stats.RData")
saveRDS(glm_results,myfile)

# Done!
message("\nDone!")
