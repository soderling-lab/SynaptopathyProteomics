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
	message(paste0("\nAnalyzing ",tissue,"..."))
} else {
	stop(msg) 
}

## Other Parameters:
save_plots = TRUE # Should plots be saved?
clean_figsdir = TRUE # Remove existing figures?
image_format = ".pdf" # Output figure format.
oldham_threshold = -2.5 # Threshold for detecting sample level outliers.
save_work = TRUE # Save workspace at end?
output_name = tissue

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
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. Load custom functions and prepare
# the project directory for saving output files.

# Load the R env.
renv::load(getrd(),quiet=TRUE)

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
rootdir <- getrd()
datadir <- file.path(rootdir, "data")
fontdir <- file.path(rootdir, "fonts")
rdatdir <- file.path(rootdir, "rdata")
figsdir <- file.path(rootdir, "figs","Data-preprocessing",tissue)

# Remove any existing figures and tables.
if (clean_figsdir) {
	invisible(sapply(list.files(figsdir,full.names=TRUE),unlink))
}

# Globally set ggplots theme.
ggtheme()

# Utilize arial font.
set_font("Arial", font_path = fontdir)

#---------------------------------------------------------------------
## Load the raw data and sample info (traits).
#---------------------------------------------------------------------
# The raw peptide intensity data were exported from ProteomeDiscover (PD)
# version 2.2. Note that the default export from PD2.x is a unitless signal to
# noise ratio, and it is not recommended to use ths for quantification.


# Load the TMT data.
datafile <- c(
	      "Cortex" = "4227_TMT_Cortex_Combined_PD_Peptide_Intensity.csv",
	      "Striatum" = "4227_TMT_Striatum_Combined_PD_Peptide_Intensity.csv"
	      )

# Load sample information.
samplefile <- c(
		"Cortex" = "4227_TMT_Cortex_Combined_traits.csv",
		"Striatum" = "4227_TMT_Striatum_Combined_traits.csv"
		)

# Load the data from PD and sample info.
raw_peptide <- fread(file = file.path(datadir, datafile[tissue]))
sample_info <- fread(file = file.path(datadir, samplefile[tissue]))

# Insure traits are in matching order.
sample_info <- sample_info[order(sample_info$Order), ]

#---------------------------------------------------------------------
## Examine an example peptide.
#---------------------------------------------------------------------

# Generate a plot.
plot <- ggplotPeptideBarPlot(raw_peptide)
plot <- plot + theme(panel.background=element_blank())
plot <- plot + theme(panel.border = element_blank(), axis.line = element_line())
plot <- plot + scale_x_discrete(expand = c(0, 0))
plot <- plot + scale_y_continuous(expand = c(0, 0)) 

# Save.
myfile <- file.path(figsdir, paste0("Example_Peptide",image_format))
if (save_plots) { 
	ggsave(prefix_file(myfile,width=2), plot, height=5, width=5) 
}

#---------------------------------------------------------------------
## Examine peptide and protein level identification overalap.
#---------------------------------------------------------------------
# Approximately 40,000 unique peptides cooresponding to ~3,000 proteins are
# quantified across all four experiments. When comparing the peptides
# identified in each experiment, the overlap is only ~25%. ~20$ of peptides are
# identified in all four experiments. This means that in different experiments,
# the same protein will likely be quantified by different peptides.
# This is a major reason the internal reference scaling (IRS) normalizaiton
# approach is employed below.

# Determine the total number of unique peptides:
nPeptides <- format(length(unique(raw_peptide$Sequence)), big.mark = ",")
message(paste("\nTotal number of unique peptides identified:",nPeptides))

# Determine the total number of unique proteins:
nProteins <- format(length(unique(raw_peptide$Accession)), big.mark = ",")
message(paste("\nTotal number of unique proteins identified:",nProteins))

# Utilize gridExtra to create a table.
tab1 <- tableGrob(data.frame(nPeptides = nPeptides, nProteins = nProteins),
  rows = NULL,
  theme = ttheme_default()
)

# Examine the number of peptides per protein.
nPep <- raw_peptide %>%
  group_by(Accession) %>%
  summarize(nPeptides = length(Sequence))

# Remove one hit wonders!
# out <- nPep$nPeptides == 1
# message(paste(sum(out),"proteins identified by a single protein are removed."))
# nPep <- nPep[!out, ]

# Generate plot.
p1 <- ggplot(nPep, aes(nPeptides)) +
  geom_histogram(bins = 100, fill = "black") +
  ggtitle("Peptides per Protein") +
  xlab("Peptides") + ylab("Frequency")
p1 <- p1 + theme(panel.background=element_blank())
p1 <- p1 + theme(panel.border = element_blank(), axis.line = element_line())
p1 <- p1 + scale_x_continuous(expand = c(0, 0))
p1 <- p1 + scale_y_continuous(expand = c(0, 0)) 
p1 <- p1 + geom_vline(xintercept = median(p1$data$nPeptides),colour='red',linetype="dotted")
p1 <- p1 + annotate("text",x=median(p1$data$nPeptides)+20,y=600,label=paste("Median =",median(p1$data$nPeptides)))

# Peptide identification overlap per pairwise comparisons of experiments.
contrasts <- combn(c("Shank2", "Shank3", "Syngap1", "Ube3a"), 2)
info_cols <- c(1, 2, 3, 4, 5)
overlap <- peptide_overlap(raw_peptide, contrasts, info_cols)
all <- c(
  "NaN", "NaN", dim(na.omit(raw_peptide))[1], dim(raw_peptide)[1],
  round(100 * dim(na.omit(raw_peptide))[1] / dim(raw_peptide)[1], 4)
)
overlap <- cbind(overlap, all)
rownames(overlap) <- c("Exp1", "Exp2", "Intersection", "Union", "Percent")
colnames(overlap) <- c(paste(c(1, 2, 3, 4, 5, 6), sep = " "), "All")

# Create table showing pairwise comparisons.
mytable <- as.data.frame(t(overlap))
mytable <- tibble::add_column(mytable, rownames(mytable), .before = 1)
mytable$Intersection <- formatC(as.numeric(mytable$Intersection), 
				format = "d", big.mark = ",")
mytable$Union <- formatC(as.numeric(mytable$Union), format = "d", big.mark = ",")
mytable$Percent <- round(as.numeric(mytable$Percent), 2)
colnames(mytable)[1] <- "Comparison"
tab2 <- tableGrob(mytable, rows = NULL, theme = ttheme_default())

# Plot peptide identification overlap.
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")
p2 <- ggplotFreqOverlap(raw_peptide, "Abundance", groups) +
  labs(title = "Peptide Identification\nOverlap")
p2 <- p2 + theme(panel.background=element_blank())
p2 <- p2 + theme(panel.border = element_blank(), axis.line = element_line())
p2 <- p2 + scale_x_discrete(expand = c(0, 0))
p2 <- p2 + scale_y_discrete(expand = c(0, 0)) 

# Save tables and plots.
if (save_plots) {
	# Tab1.
	myfile <- prefix_file(file.path(figsdir,
					paste0("Raw_Peptide_Identification",
					       image_format)))
	ggsaveTable(tab1,myfile)
	# Tab2.
	myfile <- prefix_file(file.path(figsdir,
					paste0("Peptide_Identification_Table",
					       image_format)))
	ggsaveTable(tab2,myfile)
	# Plot1.
	myfile <- prefix_file(file.path(figsdir,
				paste0("Peptide_Histogram",
				       image_format)))
	ggsave(myfile,p1,height=5,width=5)
	# Plot2.
	myfile <- prefix_file(file.path(figsdir,
				paste0("Peptide_Overlap",
				       image_format)))
	ggsave(myfile,p2,height=5,width=5)
}

#---------------------------------------------------------------------
## Examine the raw data.
#---------------------------------------------------------------------
# The need for normalization is evident in the raw data. Note that in 
# the MDS plot, samples cluster by experiment--this is evidence of 
# a batch effect.

# Prepare the data.
data_in <- raw_peptide
title <- "Raw Peptide"
colors <- c(rep("green", 11), rep("purple", 11), 
	    rep("yellow", 11), rep("blue", 11))

# Generate boxplot
p1 <- ggplotBoxPlot(data_in, colID = "Abundance", colors, title)
l1 <- get_legend(p1) # extracts legend.
p1 <- p1 + theme(legend.position = "none")
p1 <- p1 + theme(axis.text.x = element_blank())
p1 <- p1 + theme(panel.background=element_blank())
p1 <- p1 + theme(panel.border = element_blank(), axis.line = element_line())
p1 <- p1 + scale_x_discrete(expand = c(0, 0))
p1 <- p1 + scale_y_continuous(expand = c(0, 0)) 

# Generate density plot.
p2 <- ggplotDensity(data_in, colID = "Abundance", title) +
  theme(legend.position = "none")
p2 <- p2 + theme(panel.background=element_blank())
p2 <- p2 + theme(panel.border = element_blank(), axis.line = element_line())
p2 <- p2 + scale_x_continuous(expand = c(0, 0))
p2 <- p2 + scale_y_continuous(expand = c(0, 0)) 

# Genotype specific colors must be specified in column order.
#colors <- rep(c("yellow", "blue", "green", "purple"), each = 11)
#p2 <- p2 + scale_color_manual(values = colors)

# Mean SD plot.
p3 <- ggplotMeanSdPlot(data_in, colID = "Abundance", title, log = TRUE)
p3 <- p3 + theme(panel.background=element_blank())
p3 <- p3 + theme(panel.border = element_blank(), axis.line = element_line())
p3 <- p3 + scale_x_continuous(expand = c(0, 0))
p3 <- p3 + scale_y_continuous(expand = c(0, 0)) 

# Generate PCA plot.
# FIXME: add percent to axes.
colors <- rep(c("yellow", "blue", "green", "purple"), each = 11)
p4 <- ggplotPCA(data_in, traits = sample_info, colors, title = "2D PCA Plot") +
  theme(legend.position = "none")
p4 <- p4 + theme(panel.background=element_blank())
p4 <- p4 + theme(panel.border = element_rect(fill=NA))
p4 <- p4 + scale_x_continuous(expand = c(0, 0))
p4 <- p4 + scale_y_continuous(expand = c(0, 0)) 

# Save plots.
if (save_plots) {
	file_names <- paste0("Raw_Peptide_",
			     c("Boxplot","Density_plot","MeanSD_plot","PCA_plot"))
	myfiles <- file.path(figsdir,paste0(file_names,image_format))
	plots <- list(p1,p2,p3,p4)
	ggsavePlots(plots,prefix_file(myfiles),height=5,width=5)
}

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
## Examine the SL Data.
#---------------------------------------------------------------------
# Sample loading normalization equalizes the run-level sums within
# an 11-plex TMT experiment.

# Generate boxplot.
title <- "SL Peptide"
data_in <- SL_peptide
colors <- c(rep("green", 11), rep("purple", 11), 
	    rep("yellow", 11), rep("blue", 11))
p1 <- ggplotBoxPlot(data_in, colID = "Abundance", colors, title) + 
	theme(legend.position = "none")
p1 <- p1 + theme(axis.text.x = element_blank())
p1 <- p1 + theme(legend.position = "none")
p1 <- p1 + theme(axis.text.x = element_blank())
p1 <- p1 + theme(panel.background=element_blank())
p1 <- p1 + theme(panel.border = element_blank(), axis.line = element_line())
p1 <- p1 + scale_x_discrete(expand = c(0, 0))
p1 <- p1 + scale_y_continuous(expand = c(0, 0)) 

# Generate density plot.
p2 <- ggplotDensity(data_in, colID = "Abundance", title) + 
	theme(legend.position = "none")
#colors <- c(rep("yellow", 11), rep("blue", 11), 
#	    rep("green", 11), rep("purple", 11))
#p2 <- p2 + scale_color_manual(values = colors)
p2 <- p2 + theme(panel.background=element_blank())
p2 <- p2 + theme(panel.border = element_blank(), axis.line = element_line())
p2 <- p2 + scale_x_continuous(expand = c(0, 0))
p2 <- p2 + scale_y_continuous(expand = c(0, 0)) 

# Generate meanSD plot.
p3 <- ggplotMeanSdPlot(data_in, colID = "Abundance", title, log = TRUE)
p3 <- p3 + theme(panel.background=element_blank())
p3 <- p3 + theme(panel.border = element_blank(), axis.line = element_line())
p3 <- p3 + scale_x_continuous(expand = c(0, 0))
p3 <- p3 + scale_y_continuous(expand = c(0, 0)) 

# Generate PCA plot.
colors <- rep(c("yellow", "blue", "green", "purple"), each = 11)
p4 <- ggplotPCA(data_in, traits = sample_info, colors, title = "2D PCA Plot") +
  theme(legend.position = "none")
p4 <- p4 + theme(panel.background=element_blank())
p4 <- p4 + theme(panel.border = element_rect(fill=NA))
p4 <- p4 + scale_x_continuous(expand = c(0, 0))
p4 <- p4 + scale_y_continuous(expand = c(0, 0)) 

# Save plots.
if (save_plots) {
	file_names <- paste0("SL_Peptide_",
			     c("Boxplot","Density_plot",
			       "MeanSD_plot","PCA_plot"))
	myfiles <- file.path(figsdir,paste0(file_names,image_format))
	ggsavePlots(list(p1,p2,p3,p4),prefix_file(myfiles))
}

#---------------------------------------------------------------------
## Illustrate the mean variance relationship of QC peptides.
#---------------------------------------------------------------------
# Quality control samples can be used to asses intra-experimental variability.
# Peptides that have highly variable QC measurements will increase protein
# level variability and should be removed. The peptide-level QC data are
# binned by intensity, and peptides whose mean ratio are more than four
# standard deviations away from the mean are considered outliers and
# removed.

# Generate QC correlation scatter plots for all experimental groups.
# FIXME: Error in min(xrange) : invalid 'type' (list) of argument
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")
plots <- ggplotCorQC(SL_peptide, groups, colID = "QC", nbins = 5)

# Generate intensity bin histograms.
# This will take a couple minutes.
hist_list <- lapply(as.list(groups), function(x) {
  ggplotQCHist(SL_peptide, x, nbins = 5, threshold = 4)
})
names(hist_list) <- groups

# Save plots.
if (save_plots) {
	file_names <- paste(rep(names(hist_list),each=5),
			    "QC_Hist",c(1,2,3,4,5),sep="_")
	myfiles <- file.path(figsdir,paste0(file_names,image_format))
	plots <- unlist(hist_list,recursive=FALSE)
	suppressWarnings({ggsavePlots(plots,prefix_file(myfiles),
		height=5,width=5)})
}

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
if (save_plots) {
	myfile <- file.path(figsdir,
			    paste0("N_Peptides_Removed_Table",image_format))
	ggsaveTable(mytable,prefix_file(myfile))
}

#---------------------------------------------------------------------
## Examine the nature of missing values.
#---------------------------------------------------------------------
# Missing values are inherent in high throughput experiments. There are two
# main classes of missing values, missing at random (MAR) and missing not at
# random (MNAR). The appropriate imputing algorithm should be chosen based on
# the nature of missing values. The distribution of peptides with missing
# values is examined by density plots. The left-shifted distribution of
# peptides with missing values indicates that peptides that have missing
# values are generally lower in abundance. Missing values are likely then to be
# missing not at random (MNAR), but missing because they are low-abundance
# and at or near the limit of detection. MNAR data can be imputed with the
# k-nearest neighbors (knn) algorithm in the next chunk. The __impute.knn__
# function from the package `impute` is used to impute MNAR data.

# Define groups for subseting the data.
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")

# Generate plots.
p1 <- ggplotDetect(filter_peptide, groups[1]) #+ ggtitle(NULL)
l1 <- cowplot::get_legend(p1)
p1 <- p1 + theme(legend.position = "none")
p1 <- p1 + theme(panel.background=element_blank())
p1 <- p1 + theme(panel.border = element_blank(), axis.line = element_line())
p1 <- p1 + scale_x_continuous(expand = c(0, 0))
p1 <- p1 + scale_y_continuous(expand = c(0, 0)) 

p2 <- ggplotDetect(filter_peptide, groups[2]) #+ ggtitle(NULL)
p2 <- p2 + theme(legend.position = "none")
p2 <- p2 + theme(panel.background=element_blank())
p2 <- p2 + theme(panel.border = element_blank(), axis.line = element_line())
p2 <- p2 + scale_x_continuous(expand = c(0, 0))
p2 <- p2 + scale_y_continuous(expand = c(0, 0)) 

p3 <- ggplotDetect(filter_peptide, groups[3]) #+ ggtitle(NULL)
p3 <- p3 + theme(legend.position = "none")
p3 <- p3 + theme(panel.background=element_blank())
p3 <- p3 + theme(panel.border = element_blank(), axis.line = element_line())
p3 <- p3 + scale_x_continuous(expand = c(0, 0))
p3 <- p3 + scale_y_continuous(expand = c(0, 0)) 

p4 <- ggplotDetect(filter_peptide, groups[4]) #+ ggtitle(NULL)
p4 <- p4 + theme(legend.position = "none")
p4 <- p4 + theme(panel.background=element_blank())
p4 <- p4 + theme(panel.border = element_blank(), axis.line = element_line())
p4 <- p4 + scale_x_continuous(expand = c(0, 0))
p4 <- p4 + scale_y_continuous(expand = c(0, 0)) 

# Save plots.
if (save_plots) {
	file_names <- paste0(groups,"_Protein_Missing_Value_Density")
	myfiles <- file.path(figsdir,paste0(file_names,image_format))
	ggsavePlots(list(p1,p2,p3,p4),prefix_file(myfiles),
		    height=5,width=5)
}

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
if (save_plots) {
	myfile <- file.path(figsdir,paste0("N_Imputed_Peptides",image_format))
	ggsaveTable(mytable,prefix_file(myfile))
}

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


#FIXME: plots do not look correct!

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
if (save_plots) {
	myfile <- file.path(figsdir,
			    paste0("Batch_Effect_Quant_Table",image_format))
	ggsaveTable(mytable,prefix_file(myfile))
}

#---------------------------------------------------------------------
##  Examine protein identification overlap.
#---------------------------------------------------------------------
# Approximately 80-90% of all proteins are identified in all experiments.

# Inspect the overlap in protein identifcation.
plot <- ggplotFreqOverlap(combat_protein, "Abundance", groups) +
  ggtitle("Protein Identification Overlap")
plot <- plot + theme(legend.position = "none")
plot <- plot + theme(panel.background=element_blank())
plot <- plot + theme(panel.border = element_blank(), axis.line = element_line())
plot <- plot + scale_x_discrete(expand = c(0, 0))
plot <- plot + scale_y_continuous(expand = c(0, 0)) 

# Save.
if (save_plots) {
	myfile <- file.path(figsdir,
			    paste0("Protein_Identification_Overlap",image_format))
	ggsave(prefix_file(myfile),plot,height=5,width=5)
}

#---------------------------------------------------------------------
## Examine the Normalized protein level data.
#---------------------------------------------------------------------

# Generate boxplot.
data_in <- combat_protein
title <- "Normalized protein"
colors <- c(rep("green", 11), rep("purple", 11), 
	    rep("yellow", 11), rep("blue", 11))
p1 <- ggplotBoxPlot(data_in, colID = "Abundance", colors, title)
p1 <- p1 + theme(legend.position = "none")
p1 <- p1 + theme(axis.text.x = element_blank())
p1 <- p1 + theme(panel.background=element_blank())
p1 <- p1 + theme(panel.border = element_blank(), axis.line = element_line())
p1 <- p1 + scale_x_discrete(expand = c(0, 0))
p1 <- p1 + scale_y_continuous(expand = c(0, 0)) 

# Generate density plot.
p2 <- ggplotDensity(data_in, colID = "Abundance", title) + 
	theme(legend.position = "none")
#colors <- c(rep("yellow", 11), rep("blue", 11), 
#	    rep("green", 11), rep("purple", 11))
#p2 <- p2 + scale_color_manual(values = colors)
p2 <- p2 + theme(panel.background=element_blank())
p2 <- p2 + theme(panel.border = element_blank(), axis.line = element_line())
p2 <- p2 + scale_x_continuous(expand = c(0, 0))
p2 <- p2 + scale_y_continuous(expand = c(0, 0)) 

# Generate meanSd plot.
p3 <- ggplotMeanSdPlot(data_in, colID = "Abundance", title, log = TRUE)
p3 <- p3 + theme(panel.background=element_blank())
p3 <- p3 + theme(panel.border = element_blank(), axis.line = element_line())
p3 <- p3 + scale_x_continuous(expand = c(0, 0))
p3 <- p3 + scale_y_continuous(expand = c(0, 0)) 

# Generate PCA plot.
colors <- rep(c("yellow", "blue", "green", "purple"), each = 11)
p4 <- ggplotPCA(data_in, traits = sample_info, colors, title = "2D PCA Plot") +
  theme(legend.position = "none")

# Save the plots.
if (save_plots) {
	file_names <- paste0("Norm_Protein_",
			     c("Boxplot","Density_plot","MeanSD_plot","PCA_plot"))
	myfiles <- file.path(figsdir,paste0(file_names,image_format))
	ggsavePlots(list(p1,p2,p3,p4),prefix_file(myfiles))
}

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

# Save plots for outliers.
if (save_plots) {
	idx <- which(unlist(out_samples) != "none")
	if (length(idx) > 0) {
		plots <- plots[which(unlist(out_samples) != "none")]
		file_names <- file.path(figsdir,
				paste0("Outlier_Sample_",c(1:length(plots))))
		myfiles <- paste0(file_names,image_format)
		ggsavePlots(plots,prefix_file(myfiles))
	}
}

#---------------------------------------------------------------------
## Examine the IRS Normalized protein level data.
#---------------------------------------------------------------------

# Colors for plots.
# Can handle chaning number of samples if any removed as outliers.
data_in <- IRS_protein
title <- "IRS Normalized protein"
n <- sapply(groups,function(x) length(grep(x,colnames(data_in))))
colors <- unlist(mapply(function(x,y) rep(x,each=y), 
			c("yellow", "blue", "green", "purple"), n,
			SIMPLIFY = FALSE),use.names=FALSE)

# Generate boxplot.
p1 <- ggplotBoxPlot(data_in, colID = "Abundance", colors, title)

# Generate density plot.
p2 <- ggplotDensity(data_in, colID = "Abundance", title) +
  theme(legend.position = "none")
p2 <- p2 + scale_color_manual(values = colors)

# Generate meanSd plot.
p3 <- ggplotMeanSdPlot(data_in, colID = "Abundance", title, log = TRUE)
p3 <- p3 + theme(panel.background=element_blank())
p3 <- p3 + theme(panel.border = element_blank(), axis.line = element_line())
p3 <- p3 + scale_x_continuous(expand = c(0, 0))
p3 <- p3 + scale_y_continuous(expand = c(0, 0)) 

# Generate PCA plot.
p4 <- ggplotPCA(data_in, traits = sample_info, colors, title = "2D PCA Plot") +
  theme(legend.position = "none")
p4 <- p4 + theme(panel.background=element_blank())
p4 <- p4 + theme(panel.border = element_rect(fill=NA))
p4 <- p4 + scale_x_continuous(expand = c(0, 0))
p4 <- p4 + scale_y_continuous(expand = c(0, 0)) 

# Save the plots.
if (save_plots) {
	file_names <- paste0("IRS_Protein_",
			     c("Boxplot","Density_plot","MeanSD_plot","PCA_plot"))
	myfiles <- file.path(figsdir,paste0(file_names,image_format))
	ggsavePlots(list(p1,p2,p3,p4),prefix_file(myfiles),height=5,width=5)
}

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

# Generate plot to examine distribution of remaining missing values.
plot <- ggplotDetect(filter_protein, "Abundance") +
  ggtitle("Protein missing value distribution")
plot <- plot + theme(panel.background=element_blank())
plot <- plot + theme(panel.border = element_blank(), axis.line = element_line())
plot <- plot + scale_x_continuous(expand = c(0, 0))
plot <- plot + scale_y_continuous(expand = c(0, 0)) 

# Impute the remaining number of missing values with KNN.
message("\nImputing missing protein values...")
impute_protein <- impute_proteins(filter_protein, "Abundance", method = "knn")

# Save plot.
if (save_plots) {
	myfile <- file.path(figsdir,paste0("Protein_Missing_Value_Density",
					   image_format))
	ggsave(prefix_file(myfile),plot,height=5,width=5)
}

#---------------------------------------------------------------------
# Reformat data for TAMPOR normalization script.
#---------------------------------------------------------------------
# Reformat final normalized, intra-batch regressed, IRS-normalized, filtered,
# and imputed data for TAMPOR Normalization.

# Data is...
data_in <- impute_protein
rownames(data_in) <- impute_protein$Accession
data_in$Accession <- NULL
data_in$Peptides <- NULL
data_in <- as.matrix(data_in)

# Extract Gene descriptions and uniprot accesssions for renaming rows.
idx <- match(rownames(data_in), SL_peptide$Accession)
descriptions <- SL_peptide$Description[idx]
uniprot <- SL_peptide$Accession[idx]

# Split at "GN=", extract the second element.
long_names <- sapply(strsplit(descriptions, "GN=", fixed = TRUE), "[", 2)

# Split at " ", extract the first element.
gene_names <- sapply(strsplit(long_names, "\\ "), "[", 1)

# Row names are Gene|Uniprot
row_names <- paste(gene_names, uniprot, sep = "|")

# Gather data, set row names
cols <- grep("Abundance", colnames(data_in))
data_out <- as.matrix(data_in[, cols])
rownames(data_out) <- row_names

# Column names are batch.channel
idx <- match(colnames(data_out), sample_info$ColumnName)
col_names <- sample_info$SampleID[idx]
colnames(data_out) <- col_names

# Reorder based on batch.
cleanDat <- data_out[, order(colnames(data_out))]

# Save raw data as Rdata.
myfile <- file.path(rdatdir, paste0(output_name, "_raw_peptide.RData"))
saveRDS(raw_peptide, myfile)

# Save cleanDat as RData.
myfile <- file.path(rdatdir, paste0(output_name, "_cleanDat.RData"))
saveRDS(cleanDat, myfile)

# Save the workspace?
if (save_work) {
	save(list = ls(all.names = TRUE), 
	     file=paste0(tissue,".RData"), envir=.GlobalEnv)
}

# Complete!
message("\nDone!")
