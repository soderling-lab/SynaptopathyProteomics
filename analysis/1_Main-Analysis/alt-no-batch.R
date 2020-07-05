#!/usr/bin/env Rscript

#' ---
#' title:
#' description: Preprocessing of TMT data.
#' authors: Tyler W Bradshaw
#' ---

## Optional Parameters:
oldham_threshold = -2.5 # Threshold for detecting sample level outliers.

#---------------------------------------------------------------------
## Overview of Data Preprocessing:
#---------------------------------------------------------------------

# [INTRA-Experiment processing]
#     Peptide-level operations:
#     |* SL normalization
#     |* Remove QC outliers
#     |* Impute missing values 
#     Protein-level proessing:
#     |* Summarize proteins
#     |* SL normalization
#     |* Remove intra-batch batch-effect.
# [INTER-Experiment operations]
#     |* Remove QC sample outliers
#     |* IRS normalization
#     |* Filter proteins
#     |* Protein imputing

# Statistical analysis:
#    * Final TMM normalization
#    * GLM for differential abundance.

#---------------------------------------------------------------------
## Misc function - getrd().
#---------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
	in_root <- function(h=here, dir=dpat) { 
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE))) 
		return(check)
	}
	# Loop to find root.
	while (!in_root(here)) { 
		here <- dirname(here) 
	}
	root <- here
	return(root)
}

# Parse the command line arguments.
parse_args <- function(default="Cortex", args=commandArgs(trailingOnly=TRUE)){
	# Input must be Cortex or Striatum.
	msg <- c("Please specify a tissue type to be analyzed:\n",
	 "Choose either 'Cortex' or 'Striatum'.")
	# If interactive, return default tissue.
	if (interactive()) { 
		return("Cortex") 
	} else {
		# Check arguments.
		check <- !is.na(match(args[1], c("Cortex", "Striatum")))
		if (length(args == 1) & check) { 
			tissue  <- args[1]
			start <- Sys.time()
			message(paste("Starting analysis at:", start))
			message(paste0("Analyzing ", tissue,"..."))
		} else {
			stop(msg) 
		}
		return(tissue)
	}
}

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. Load custom functions and prepare
# the project directory for saving output files.

# Parse input.
tissue <- parse_args()

# Load the R env.
root <- getrd()
renv::load(root,quiet=TRUE)

# Load required packages.
suppressPackageStartupMessages({
  library(sva)
  library(grid)
  library(dplyr)
  library(WGCNA)
  library(edgeR)
  library(tibble)
  library(gtable)
  library(cowplot)
  library(ggplot2)
  library(gridExtra)
  library(flashClust)
  library(data.table)
})

# Load project specific data and functions.
suppressWarnings({ devtools::load_all() })

# Project directories.
datadir <- file.path(root,"data") # Input/Output data.
rdatdir <- file.path(root,"rdata") # Temporary data files.

#---------------------------------------------------------------------
## Load the raw data and sample info (traits).
#---------------------------------------------------------------------
# The raw peptide intensity data were exported from ProteomeDiscover (PD)
# version 2.2. Note that the default export from PD2.x is a unitless signal to
# noise ratio, and it is not recommended to use ths for quantification.

# Load the TMT data.
data_files <- c(
	      "Cortex" = "Cortex_Raw_Peptides.csv",
	      "Striatum" = "Striatum_Raw_Peptides.csv"
	      )

# Load sample information.
meta_files <- c(
		"Cortex" = "Cortex_Samples.csv",
		"Striatum" = "Striatum_Samples.csv"
		)

# Load the TMT data and sample meta data.
raw_peptide <- fread(file = file.path(datadir, data_files[tissue]))
sample_info <- fread(file = file.path(datadir, meta_files[tissue]))

#---------------------------------------------------------------------
## Create gene identifier map.
#---------------------------------------------------------------------

# Create gene map.
message("\nCreating gene identifier map.")
uniprot <- raw_peptide$Accession
entrez <- mgi_batch_query(ids=uniprot)
symbols <- getPPIs::getIDs(entrez,from="entrez",to="symbol",species="mouse")
gene_map <- as.data.table(keep.rownames="uniprot",entrez)
gene_map$symbol <- symbols[as.character(gene_map$entrez)]

# Save gene_map as rda.
myfile <- file.path(datadir,"gene_map.rda")
save(gene_map,file=myfile,version=2)

#---------------------------------------------------------------------
## Sample loading normalization within experiments.
#---------------------------------------------------------------------
# The function __normalize_SL__ performs sample loading (SL) normalization to
# equalize the run level intensity (column) sums. The data in each column are
# multiplied by a factor such that the mean of the column sums are are equal.
# Sample loading normalization is performed within an experiment under the
# assumption that equal amounts of protein were used for each of the 11 TMT
# channels.

# Define data columns for SL within experiments:
colID <- "Abundance"
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")

# Perform SL normalization.
SL_peptide <- normalize_SL(raw_peptide, colID, groups)

#---------------------------------------------------------------------
## Peptide level filtering based on QC samples.
#---------------------------------------------------------------------
# Peptides that were not quantified in all three qc replicates are removed.
# The data are binned by intensity, and measurments that are 4xSD from the mean
# ratio of the intensity bin are considered outliers and removed.

# Filter peptides based on QC precision.
message("\nRemoving outlier QC peptides...")
filter_peptide <- filter_QC(SL_peptide, groups, nbins = 5, threshold = 4)

# Generate table.
# Summarize number of Cortex and Striatum peptides removed.
out <- list("Cortex"=c(94, 78, 134, 59), 
	    "Striatum"=c(182, 67, 73, 75))[[tissue]] 
mytable <- data.frame(cbind(groups, out))
mytable <- tableGrob(mytable, rows = NULL, theme = ttheme_default())

# Save the table.
#myfile <- file.path(figsdir,
#		    paste0("N_Peptides_Removed_Table",image_format))
#ggsaveTable(mytable,prefix_file(myfile))

#---------------------------------------------------------------------
## Impute missing peptide values within an experiment.
#---------------------------------------------------------------------
# The function __impute_Peptides__ supports imputing missing values with the
# maximum likelyhood estimation (MLE) or KNN algorithms for missing not at
# random (MNAR) and missing at random data, respectively. Impution is performed
# with an experiment, and rows with more than 50% missing values are censored
# and will not be imputed. Peptides with more than 2 missing biological
# replicates or any missing quality control (QC) replicates will be
# censored and are not imputed.
message("\nImputing missing peptide values...")

# Impute missing values using KNN algorithm for MNAR data.
# Rows with missing QC replicates are ingored (qc_threshold=0).
# Rows with more than 2 (50%) missing biological replicates are
# ignored (bio_threshold=2).
data_impute <- impute_peptides(filter_peptide, groups, method = "knn")
impute_peptide <- data_impute$data_imputed

# Table of n imputed peptides.
n_out <- data_impute$n_out
mytable <- as.data.frame(do.call(rbind, n_out))
mytable <- tibble::add_column(mytable, rownames(mytable), .before = 1)
colnames(mytable) <- c("Experiment", "N Imputed")
mytable <- tableGrob(mytable, rows = NULL, theme = ttheme_default())

# Save table.
#myfile <- file.path(figsdir,paste0("N_Imputed_Peptides",image_format))
#ggsaveTable(mytable,prefix_file(myfile))

#---------------------------------------------------------------------
## Protein level summarization and normalization across all batches.
#---------------------------------------------------------------------
# Summarize to protein level by summing peptide intensities. Note that the
# peptide column in the returned data frame reflects the total number of
# peptides identified for a given protein across all experiments.

# Summarize to protein level:
raw_protein <- summarize_protein(impute_peptide)

# Normalize across all columns.
SL_protein <- normalize_SL(raw_protein, colID="Abundance", group="Abundance")

#---------------------------------------------------------------------
## IRS Normalization.
#---------------------------------------------------------------------
# Internal reference sclaing (IRS) normalization equalizes the protein-wise means
# of reference (QC) samples across all batches. Thus, IRS normalization accounts
# for the random sampling of peptides at the MS2 level which results in the
# identification/quantificaiton of proteins by different peptides in each
# experiment. IRS normalization was first described by __Plubell et al., 2017__.
message("\nPerforming IRS normalization...")

# Perform IRS normalization.
IRS_protein <- normalize_IRS(SL_protein, "QC", groups, robust = TRUE)

#---------------------------------------------------------------------
## Identify and remove QC sample outliers.
#---------------------------------------------------------------------
# IRS normalization utilizes QC samples as reference samples. Outlier QC
# measurements (caused by interference or other artifact) would influence the
# create unwanted variability. Thus, outlier QC samples are removed, if
# identified. The method used by __Oldham et al., 2016__ is used to identify
# QC sample outliers. A threshold of -2.5 is used.

# Data is...
data_in <- IRS_protein

# Illustrate Oldham's sample connectivity.
sample_connectivity <- ggplotSampleConnectivity(data_in,
  colID = "QC",
  threshold = oldham_threshold
)

# Check sample connectivity
tab <- sample_connectivity$table
df <- tibble::add_column(tab, "Name" = rownames(tab), .before = 1)
rownames(df) <- NULL
knitr::kable(df)

# Loop to identify Sample outliers using Oldham's connectivity method.
n_iter <- 5
data_in <- IRS_protein
out_samples <- list()
plots <- list()

# Loop:
for (i in 1:n_iter) {
  data_temp <- quiet(normalize_IRS(data_in, "QC", groups, robust = TRUE))
  oldham <- ggplotSampleConnectivity(data_temp, log = TRUE, colID = "QC")
  plots[[i]] <- oldham$connectivityplot +
    ggtitle(paste("Sample Connectivity (Iteration = ", i, ")", sep = ""))
  bad_samples <- rownames(oldham$table)[oldham$table$Z.Ki < oldham_threshold]
  message(paste(
    length(bad_samples), " outlier sample(s) identified in iteration ", i, ".",
    sep = ""
  ))
  if (length(bad_samples) == 0) bad_samples <- "none"
  out_samples[[i]] <- bad_samples
  out <- grepl(paste(unlist(out_samples), collapse = "|"), colnames(data_in))
  data_in <- quiet(normalize_IRS(data_in[, !out], "QC", groups, robust = TRUE))
}

# Outlier samples.
bad_samples <- unlist(out_samples)
message(paste("\nTotal number of outlier QC samples identified:", 
	      sum(bad_samples != "none")))

# Remove outliers from data.
samples_out <- paste(bad_samples, collapse = "|")
out <- grepl(samples_out, colnames(SL_protein))

# Redo IRS after outlier removal.
IRS_OutRemoved_protein <- normalize_IRS(SL_protein[, !out],
  "QC", groups,
  robust = TRUE
)

# Write over IRS_data
IRS_protein <- IRS_OutRemoved_protein

#---------------------------------------------------------------------
## Protein level filtering, and imputing.
#---------------------------------------------------------------------
# Proteins that are identified by only a single peptide are removed. Proteins
# that are identified in less than 50% of all samples are also removed. The
# nature of the remaining missng values are examined by density plot and
# imputed with the KNN algorithm for MNAR data.
message("\nRemoving irroproducibly quantified proteins...")

# Remove proteins that are identified by only 1 peptide as well as
# proteins identified in less than 50% of samples.
filter_protein <- filter_proteins(IRS_protein, "Abundance")

# Impute the remaining number of missing values with KNN.
message("\nImputing missing protein values...")
impute_protein <- impute_proteins(filter_protein, "Abundance", method = "knn")

#---------------------------------------------------------------------
## Save the proteomics data in tidy format.
#---------------------------------------------------------------------

# Tidy.
tidy_prot <- reshape2::melt(impute_protein,id.var=c("Accession","Peptides"),
			    value.var="Intensity", variable.name="Sample",
			    value.name = "Intensity") %>% as.data.table()

# Merge data and meta data by sample name.
tidy_prot <- left_join(tidy_prot,sample_info,by="Sample")

#---------------------------------------------------------------------
## Create protein networks.
#---------------------------------------------------------------------

message("\nBuilding protein networks.")

# Drop QC and do bicor.
dm <- tidy_prot %>% filter(!grepl("QC",Sample)) %>% 
	as.data.table() %>% 
	dcast(Sample ~ Accession, value.var="Intensity") %>% 
	as.matrix(rownames="Sample")
adjm <- WGCNA::bicor(log2(dm))

# Neten
ne_adjm <- neten::neten(adjm)

#---------------------------------------------------------------------
## Save output.
#---------------------------------------------------------------------
message("\nSaving data.")

## Output in root/rdata

# ne_adjm.csv - adjacency matrix saved as csv for Leidenalg.
myfile <- file.path(rdatdir,paste0(tolower(tissue),"_ne_adjm.csv"))
ne_adjm %>% as.data.table(keep.rownames="Accession") %>% fwrite(myfile)

## Output in root/data

# Save adjm as rda.
myfile <-file.path(datadir,paste0(tolower(tissue),"_adjm.rda"))
save(adjm,file=myfile,version=2)

# ne_adjm.rda
myfile <-file.path(datadir, paste0(tolower(tissue),"_ne_adjm.rda"))
save(ne_adjm,file=myfile,version=2)

# tidy_prot.rda
myfile <- file.path(datadir,paste0(tolower(tissue),".rda"))
tidy_prot <- tidy_prot %>% filter(!grepl("QC",Sample)) %>% as.data.table() # Drop QC!
save(tidy_prot,file=myfile,version=2)

#--------------------------------------------------------------------
## EdgeR glm
#--------------------------------------------------------------------

# Cast tp into data matrix for EdgeR. Don't log!
groups <- unique(tidy_prot$Genotype)

# Loop to perform intra-genotype WT v MUT comparisons:
results_list <- list()
for (geno in groups){
	dm <- tidy_prot %>% filter(Genotype==geno) %>%
		dcast(Accession ~ Sample, value.var="Intensity") %>% 
		as.matrix(rownames=TRUE)
	# Create dge object.
	dge <- DGEList(counts=dm)
	# Perform TMM normalization.
	dge <- calcNormFactors(dge)
	# Sample mapping.
	samples <- rownames(dge$samples)
	idx <- match(samples,tidy_prot$Sample)
	genotype <- tidy_prot$Genotype[idx]
	treatment <- tidy_prot$Treatment[idx]
	batch <- as.factor(tidy_prot$PrepDate[idx])
	dge$samples$group <- interaction(genotype,treatment)
	# Basic design matrix for GLM -- all groups treated seperately.
	design <- model.matrix(~ batch + group, data = dge$samples)
	# Estimate dispersion:
	dge <- estimateDisp(dge, design, robust = TRUE)
	# Fit a general linear model.
	fit <- glmQLFit(dge, design, robust = TRUE)
	# Create contrast.
	#wt <- colnames(design)[grepl("WT",colnames(design))]
	#mut <- colnames(design)[grepl("HET|KO",colnames(design))]
	#contr <- limma::makeContrasts(paste(mut,wt,sep="-"),levels=design)
	# Assess differences.
	#qlf <- glmQLFTest(fit,contrast=contr)
	qlf <- glmQLFTest(fit)
	# Call topTags to add FDR. Gather tabularized results.
	glm_results <- topTags(qlf, n = Inf, sort.by = "PValue")$table
	# Use gene map to annotate glm_results with entrez Ids and gene symbols.
	idx <- match(rownames(glm_results), gene_map$uniprot)
	glm_results <- add_column(glm_results, "Accession"=rownames(glm_results),
				  .before = 1)
	glm_results <- add_column(glm_results, "Entrez" = gene_map$entrez[idx], 
				  .after = 1)
	glm_results <- add_column(glm_results, "Symbol" = gene_map$symbol[idx], 
				  .after = 2)
	# Add expression data.
	wt_cols <- grep("WT",colnames(dm))
	mut_cols <- grep("HET|KO",colnames(dm))
	dt <- as.data.table(dm[,c(wt_cols,mut_cols)],keep.rownames="Accession")
	glm_results <- left_join(glm_results,dt,by="Accession")
	results_list[[geno]] <- glm_results
}

# Save as excel.
names(results_list) <- paste(names(results_list),tissue,"Results")
myfile <- file.path(root,"tables","TMT_Protein_GLM_Results.xlsx")
write_excel(results_list,myfile)
