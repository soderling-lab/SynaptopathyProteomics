#!/usr/bin/env Rscript

#' ---
#' title:
#' description: Preprocessing of TMT data.
#' authors: Tyler W Bradshaw
#' ---

# Parse command line input:
# Analysis (tissue) type: cortex (1) or striatum(2).
args <- commandArgs(trailingOnly = TRUE)
msg <- c("Please specify a tissue type to be analyzed:\n",
	 "Choose either 'Cortex' or 'Striatum'.")
check <- !is.na(match(args[1], c("Cortex", "Striatum")))
if (length(args == 1) & check) { 
	tissue  <- args[1]
	start <- Sys.time()
	message(paste("Starting analysis at:", start))
	message(paste0("Analyzing ", tissue,"..."))
} else {
	stop(msg) 
}

## Other Parameters:
save_plots = FALSE # Should plots be saved?
clean_figsdir = FALSE # Remove existing figures?
image_format = ".pdf" # Output figure format.
oldham_threshold = -2.5 # Threshold for detecting sample level outliers.
save_work = FALSE # Save workspace at end?
fig_width = 2.5 # Width in inches of figures.
fig_height = 2.5 # Height in inches of figures.
output_name = tissue # Prefix of output files.

#---------------------------------------------------------------------
## Overview of Data Preprocessing:
#---------------------------------------------------------------------

# All done within an experiment...
#    * SL normalization
#    * Peptide QC filtering (remove highly variable peptides)
#    * Peptide imputing (replace small number of missing values)

#    * Protein summarization
#    * SL normalization
#    * Intra-Batch ComBat

# All done across experiments...
#    * Outlier QC sample removal
#    * IRS normalization
#    * Protein filtering
#    * Protein imputing

# Combine cortex and striatum data:
#    * TAMPOR normalization
#    * Outlier sample removal

## Statistical analysis:
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

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. Load custom functions and prepare
# the project directory for saving output files.

# Load the R env.
rootdir <- getrd()
renv::load(rootdir,quiet=TRUE)

# Load required packages.
suppressPackageStartupMessages({
  library(sva)
  library(grid)
  library(dplyr)
  library(WGCNA)
  library(gtable)
  library(cowplot)
  library(ggplot2)
  library(gridExtra)
  library(flashClust)
  library(data.table)
})

# Load functions in root/R.
suppressWarnings({ devtools::load_all() })

# Directories:
fontdir <- file.path(rootdir, "fonts")
rdatdir <- file.path(rootdir, "rdata")
rawddir <- file.path(rootdir, "raw-data")
figsdir <- file.path(rootdir, "figs","Data-preprocessing",tissue)

# Create directory for figures if it doesn't exist.
if (!dir.exists(figsdir)) { dir.create(figsdir) }

# Remove any existing figures and tables.
if (clean_figsdir) {
	invisible(sapply(list.files(figsdir,full.names=TRUE),unlink))
}

# Globally set ggplots theme.
# Utilize arial font.
ggtheme(); set_font("Arial", font_path = fontdir)

#---------------------------------------------------------------------
## Load the raw data and sample info (traits).
#---------------------------------------------------------------------
# The raw peptide intensity data were exported from ProteomeDiscover (PD)
# version 2.2. Note that the default export from PD2.x is a unitless signal to
# noise ratio, and it is not recommended to use ths for quantification.

# Load the TMT data.
datafile <- c(
	      "Cortex" = "Cortex_Raw_Peptides.csv",
	      "Striatum" = "Striatum_Raw_Peptides.csv"
	      )

# Load sample information.
samplefile <- c(
		"Cortex" = "Cortex_Samples.csv",
		"Striatum" = "Striatum_Samples.csv"
		)

# Load the data from PD and sample info.
raw_peptide <- fread(file = file.path(rawddir, datafile[tissue]))
sample_info <- fread(file = file.path(rawddir, samplefile[tissue]))

# Insure traits are in matching order.
sample_info <- sample_info[order(sample_info$Order), ]

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

# Define experimental groups for checking QC variability:
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")

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

# Define experimental groups for checking QC variability:
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")

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
SL_protein <- normalize_SL(raw_protein, "Abundance", "Abundance")

#---------------------------------------------------------------------
## IntraBatch Protein-lavel ComBat.
#---------------------------------------------------------------------
# Each experimental cohort of 8 was prepared in two batches. This was necessary
# because the ultra-centrifuge rotor used to prepare purified synaptosomes
# holds a maximum of 6 samples. This intra-batch batch effect was recorded for
# 6/8 experiments. Here I will utilize the __ComBat()__ function from the `sva`
# package to remove this batch effect before correcting for inter-batch batch
# effects between batches with IRS normalization. Note that in
# the absence of evidence of a batch effect 
# (not annotated or cor(PCA,batch)<0.1),
# ComBat is not applied.

# Define experimental groups and column ID for expression data.
colID <- "Abundance"
data_in <- SL_protein

# Loop to perform ComBat on intraBatch batch effect (prep date).
# If there is no known batch effect (bicor(batch,PC1)<0.1) then
# the data is returned un-regressed.

# Note: The values of QC samples are not adjusted by ComBat.
# QC samples were prepared from a seperate batch of mice and
# represent a single batch.

data_out <- list() # ComBat data.
plot_list <- list() # MDS plots.
R <- list() # Bicor stats [bicor(batch,PC1)]

# Loop:
for (i in 1:length(groups)) {
  # Meta data.
  info_cols <- data_in[, !grepl(colID, colnames(data_in))]
  # Expression data.
  group <- groups[[i]]
  cols <- grepl(group, colnames(data_in))
  data_work <- as.matrix(data_in[, cols])
  rownames(data_work) <- paste(data_in$Accession,
    c(1:nrow(data_in)),
    sep = "_"
  )
  rows_out <- apply(data_work, 1, function(x) sum(is.na(x) > 0))
  data <- data_work[!rows_out, ]
  # Get Traits info.
  idx <- match(colnames(data), sample_info$ColumnName)
  traits_sub <- sample_info[idx, ]
  rownames(traits_sub) <- traits_sub$ColumnName
  # QC Samples will be ignored.
  ignore <- is.na(traits_sub$PrepDate)
  data_QC <- data[, ignore]
  CombatInfo <- traits_sub[!ignore, ]
  data <- data[, !ignore]
  # There should be no negative values.
  if (min(data, na.rm = TRUE) < 1) {
    data[data < 1] <- 1
    warning("Warning: Expression values less than 1 will be replaced with 1.")
  }
  # Check the correlation between batch and PC1.
  pc1 <- prcomp(t(log2(data)))$x[, 1]
  batch <- as.numeric(as.factor(CombatInfo$PrepDate)) - 1
  r1 <- suppressWarnings(WGCNA::bicor(batch, pc1))
  # Check that r2 is not NA.
  if (is.na(r1)) {
    r1 <- 0
  }
  # Check, in matching order?
  if (!all(colnames(data) == CombatInfo$ColumnName)) {
    warning("Warning: Names of traits and expression data do not match.")
  }
  # Check MDS plot prior to ComBat.
  traits_sub$Sample.Model <- paste("b",
    as.numeric(as.factor(traits_sub$PrepDate)),
    sep = "."
  )
  title <- paste(gsub(" ", "", unique(traits_sub$Model)), "pre-ComBat", sep = " ")
  plot1 <- ggplotMDS(log2(data),
    colID = "b",
    title = title, traits = traits_sub
  )$plot + theme(legend.position = "none")
  plot1 <- plot1 + scale_color_manual(values = unique(traits_sub$Color))
  # Apply ComBat.
  message(paste("\nPerforming", groups[i], "ComBat..."))
  if (length(unique(CombatInfo$PrepDate)) > 1 & abs(r1) > 0.1) {
    # Create ComBat model.
    model <- model.matrix(~ as.factor(CombatInfo$SampleType),
      data = as.data.frame(log2(data))
    )
    data_ComBat <- ComBat(
      dat = log2(data),
      batch = as.vector(CombatInfo$PrepDate), mod = model, mean.only = FALSE
    )
  } else {
    # No batch effect.
    if (abs(r1) < 0.1) {
      message(paste("Error: No quantifiable batch effect!",
        "\n", "The un-regressed data will be returned."))
      data_ComBat <- log2(data)
    } else {
      # No batch effect.
      message(paste("Error: ComBat can only be applied to factors with more than two levels!",
		    "\n", "The un-regressed data will be returned."))
      data_ComBat <- log2(data)
    }
  }
  # Correlation between batch and PC1 post-ComBat.
  pc1 <- prcomp(t(data_ComBat))$x[, 1]
  batch <- as.numeric(as.factor(CombatInfo$PrepDate)) - 1
  r2 <- suppressWarnings(WGCNA::bicor(batch, pc1))
  # Check that r2 is not NA.
  if (is.na(r2)) {
    r2 <- 0
  }
  R[[i]] <- cbind(r1, r2)
  # Check MDS plot after ComBat.
  title <- paste(gsub(" ", "", unique(traits_sub$Model)), "post-ComBat", sep = " ")
  plot2 <- ggplotMDS(
    data_in = data_ComBat, colID = "Abundance",
    title = title, traits = traits_sub
  )$plot + theme(legend.position = "none")
  plot2 <- plot2 + scale_color_manual(values = unique(traits_sub$Color))
  # Un-log.
  data_ComBat <- 2^data_ComBat
  # Recombine with QC data.
  data_out[[i]] <- cbind(info_cols[!rows_out, ], data_QC, data_ComBat)
  names(data_out)[[i]] <- group
  plots <- list(plot1, plot2)
  plot_list[[i]] <- plots
  names(plot_list[[i]]) <- paste(groups[i], c("preComBat", "postComBat"))
} # Ends ComBat loop.

# Merge the data frames with purrr::reduce()
combat_protein <- data_out %>% 
	purrr::reduce(left_join, by = c(colnames(data_in)[c(1, 2)]))

# Quantifying the batch effect.
# Check bicor correlation with batch before and after ComBat.
df <- do.call(rbind, lapply(R, function(x) x[1, ]))
df <- as.data.frame(t(apply(df, 1, function(x) round(abs(as.numeric(x)), 3))))
rownames(df) <- groups
df <- tibble::add_column(df, rownames(df), .before = 1)
colnames(df) <- c("Experiment", "preComBat", "postComBat")
mytable <- tableGrob(df, rows = NULL)

# Save table quantifying batch effects.
#myfile <- file.path(figsdir,
#		    paste0("Batch_Effect_Quant_Table",image_format))
#ggsaveTable(mytable,prefix_file(myfile))

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
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")
IRS_protein <- normalize_IRS(combat_protein, "QC", groups, robust = TRUE)

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

tab <- sample_connectivity$table
df <- tibble::add_column(tab, SampleName = rownames(tab), .before = 1)
rownames(df) <- NULL

# Loop to identify Sample outliers using Oldham's connectivity method.
n_iter <- 5
data_in <- combat_protein # normalized, regressed protein values.
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
message(paste("\nTotal number of outlier QC samples identified:", sum(bad_samples != "none")))

# Remove outliers from data.
samples_out <- paste(bad_samples, collapse = "|")
out <- grepl(samples_out, colnames(SL_protein))

# Redo IRS after outlier removal..
IRS_OutRemoved_protein <- normalize_IRS(combat_protein[, !out],
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
left_join(tidy_prot,samples,by="Sample")

#---------------------------------------------------------------------
## Create protein networks.
#---------------------------------------------------------------------

# Bicor
dm <- tidy_prot %>% filter(!grepl("QC",Sample)) %>% 
	as.data.table() %>% 
	dcast(Sample ~ Accession, value.var="Intensity") %>% 
	as.matrix(rownames="Sample") %>% log2()

adjm <- WGCNA::bicor(dm)

# Neten
ne_adjm <- neten::neten(adjm)

# Save as a simple matrix.
myfile <- file.path(rdatdir,"Cortex_NE_Adjm.csv")
ne_adjm %>% as.data.table(keep.rownames="Accession") %>% fwrite(myfile)

## DO LEIDEN

#---------------------------------------------------------------------
## Module-level statistical analysis
#---------------------------------------------------------------------

## 
# Load the partition.
data(cortex_partition)
partition = cortex_partition

# Percent clustered
modules = split(partition,partition)
pclustered <- sum(partition!=0)/length(partition)
pclustered

# module glm

