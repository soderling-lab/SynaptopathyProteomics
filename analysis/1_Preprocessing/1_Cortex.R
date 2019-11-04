#!/usr/bin/env Rscript

#' ---
#' title: 1_Cortex.R
#' description: Preprocessing of Cortex Data.
#' authors: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
## Prepare the workspace.
#-------------------------------------------------------------------------------
# Prepare the R workspace for the analysis. Load custom functions and prepare
# the project directory for saving output files.

# Load required packages.
suppressPackageStartupMessages({
	#library(readxl)
  library(data.table)
#  library(reshape2)
  library(WGCNA)
  library(dplyr)
  library(gridExtra)
  library(grid)
  library(gtable)
  library(ggplot2)
#  library(cowplot)
#  library(impute)
#  library(tibble)
  library(flashClust)
#  library(ggdendro)
  library(sva)
#  library(purrr)
#  library(ggrepel)
#  library(edgeR)
})

# Define tissue type for analysis: Cortex = 1; Striatum = 2.
type <- 1
tissue <- c("Cortex", "Striatum")[type]

# Directories:
here <- getwd()
rootdir <- dirname(dirname(here))
funcdir <- file.path(rootdir, "R")
datadir <- file.path(rootdir, "input")
Rdatadir <- file.path(rootdir, "data")
figsdir <- file.path(rootdir, "figures")
tabsdir <- file.path(rootdir, "tables")

# Load required custom functions.
source_myfun <- function() {
	myfun <- list.files(funcdir,pattern=".R",full.names=TRUE)
	invisible(sapply(myfun,source))
}
source_myfun()

# Define prefix for output figures and tables.
outputMatName <- paste0("1_", tissue)

# Globally set ggplots theme.
ggtheme()

# All plots will be stored in a list and saved to RData.
output_plots <- file.path(Rdatadir, paste0(outputMatName, "_plots.Rds"))
all_plots <- list()

#-------------------------------------------------------------------------------
## Load the raw data and sample info (traits).
#-------------------------------------------------------------------------------
# The raw peptide intensity data were exported from ProteomeDiscover (PD)
# version 2.2. Note that the default export from PD2.x is a unitless signal to
# noise ratio, and it is not recommended to use ths for quantification.

# Load the data.
datafile <- c(
  "4227_TMT_Cortex_Combined_PD_Peptide_Intensity.csv",
  "4227_TMT_Striatum_Combined_PD_Peptide_Intensity.csv"
)

# Load sample information.
samplefile <- c(
  "4227_TMT_Cortex_Combined_traits.csv",
  "4227_TMT_Striatum_Combined_traits.csv"
)

# Load the data from PD and sample info.
raw_peptide <- fread(file=file.path(datadir, datafile[type]))
sample_info <- fread(file=file.path(datadir, samplefile[type]))

# Insure traits are in matching order.
sample_info <- sample_info[order(sample_info$Order), ]

# Save to Rdata. To be saved to excel later.
myfile <- file.path(Rdatadir, paste0(outputMatName, "_raw_peptide.RData"))
save(raw_peptide,file=myfile)

myfile <- file.path(Rdatadir, paste0(outputMatName, "_sample_info.RData"))
save(sample_info,file=myfile)

#-------------------------------------------------------------------------------
## Examine an example peptide.
#-------------------------------------------------------------------------------

# Get a slice of the data.
set.seed(1)
dat <- as.data.frame(subset(raw_peptide, grepl("Dlg4", raw_peptide$Description)))
rownames(dat) <- paste(dat$Accession, dat$Sequence, c(1:nrow(dat)), sep = "_")
idy <- grepl("Shank2", colnames(dat))
dat <- dat[, idy]
dat <- na.omit(dat)

# Make bar plot for given peptide.
colIDs <- gsub(",", "", sapply(strsplit(colnames(dat), "\\ "), "[", 3))
geno <- gsub(",", "", sapply(strsplit(colnames(dat), "\\ "), "[", 5))
n <- sample(nrow(dat), 1)
df <- reshape2::melt(dat[n, ])
df$Channel <- colIDs
title <- paste0("Dlg4"," >",strsplit(rownames(dat)[n], "_")[[1]][2])

plot <- ggplot(df, aes(x = Channel, y = log2(value), fill = Channel)) +
  geom_bar(stat = "identity", width = 0.9, position = position_dodge(width = 1)) +
  xlab("TMT Channel") + ylab("Log2(Raw Intensity)") +
  ggtitle(title) + 
  theme(legend.position = "none")

# Store in list.
all_plots[[paste(tissue, "Example_TMT", sep = "_")]] <- plot

#-------------------------------------------------------------------------------
## Examine peptide and protein level identification overalap.
#-------------------------------------------------------------------------------
# Approximately 40,000 unique peptides cooresponding to ~3,000 proteins are
# quantified across all four experiments. When comparing the peptides
# identified in each experiment, the overlap is only ~25%. ~20$ of peptides are
# identified in all four experiments. This means that in different experiments,
# the same protein will likely be quantified by different peptides.
# This is a major reason the internal reference scaling (IRS) normalizaiton
# approach is employed below.

# Determine the total number of unique peptides:
nPeptides <- format(length(unique(raw_peptide$Sequence)), big.mark = ",")
print(paste(nPeptides, " unique peptides identified.", sep = ""))

# Determine the total number of unique proteins:
nProteins <- format(length(unique(raw_peptide$Accession)), big.mark = ",")
print(paste(nProteins, " unique proteins identified.", sep = ""))

# Utilize gridExtra to create a table.
tab1 <- tableGrob(data.frame(nPeptides = nPeptides, nProteins = nProteins),
  rows = NULL,
  theme = ttheme_default()
)

# Examine the number of peptides per protein.
nPep <- subset(raw_peptide) %>%
  dplyr::group_by(Accession) %>%
  dplyr::summarize(nPeptides = length(Sequence))

# Remove one hit wonders!
nPep <- nPep[!nPep$nPeptides == 1, ]

# Generate plot.
plot1 <- ggplot(nPep, aes(nPeptides)) + geom_histogram(bins = 100, fill = "black") +
  ggtitle("Peptides per Protein") + xlab("Peptides") + ylab("Frequency")

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
mytable <- add_column(mytable, rownames(mytable), .before = 1)
mytable$Intersection <- formatC(as.numeric(mytable$Intersection), format = "d", big.mark = ",")
mytable$Union <- formatC(as.numeric(mytable$Union), format = "d", big.mark = ",")
mytable$Percent <- round(as.numeric(mytable$Percent), 2)
colnames(mytable)[1] <- "Comparison"
tab2 <- tableGrob(mytable, rows = NULL, theme = ttheme_default())

# Plot peptide identification overlap.
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")
plot2 <- ggplotFreqOverlap(raw_peptide, "Abundance", groups) +
  labs(title = "Peptide Identification\nOverlap")

# Store plots in list.
all_plots[[paste(tissue, "n_pep_per_protein", sep = "_")]] <- plot1
all_plots[[paste(tissue, "pep_id_overlap", sep = "_")]] <- plot2
all_plots[[paste(tissue, "total_pep_and_prot_tab", sep = "_")]] <- tab1
all_plots[[paste(tissue, "pep_id_overlap_tab", sep = "_")]] <- tab2

#-------------------------------------------------------------------------------
## Examine the raw data.
#-------------------------------------------------------------------------------
# The need for normalization is evident in the raw data. Note that in the MDS
# plot, samples cluster by experiment--evidence of a batch effect.

# Prepare the data.
data_in <- raw_peptide
title <- "Raw Peptide"
colors <- c(rep("green", 11), rep("purple", 11), rep("yellow", 11), rep("blue", 11))

# Generate boxplot
p1 <- ggplotBoxPlot(data_in, colID = "Abundance", colors, title)
l1 <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")
p1 <- p1 + theme(axis.text.x = element_blank())

# Generate density plot.
p2 <- ggplotDensity(data_in, colID = "Abundance", title) + theme(legend.position = "none")

# Genotype specific colors must be specified in column order.
colors <- rep(c("yellow","blue","green","purple"),each=11)
p2 <- p2 + scale_color_manual(values = colors)
p3 <- ggplotMeanSdPlot(data_in, colID = "Abundance", title, log = TRUE)

# Generate PCA plot.
colors <- rep(c("yellow","blue","green","purple"),each=11)
p4 <- ggplotPCA(data_in, traits=sample_info, colors, title = "2D PCA Plot") +
	theme(legend.position = "none")

# Store plots in list.
all_plots[[paste(tissue, "raw_bp", sep = "_")]] <- p1
all_plots[[paste(tissue, "raw_dp", sep = "_")]] <- p2
all_plots[[paste(tissue, "raw_msd", sep = "_")]] <- p3
all_plots[[paste(tissue, "raw_mds", sep = "_")]] <- p4

#-------------------------------------------------------------------------------
## Sample loading normalization within experiments.
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
## Examine the SL Data.
#-------------------------------------------------------------------------------
# Sample loading normalization equalizes the run-level sums within
# an 11-plex TMT experiment.

# Generate boxplot.
title <- "SL Peptide"
data_in <- SL_peptide
colors <- c(rep("green", 11), rep("purple", 11), rep("yellow", 11), rep("blue", 11))
p1 <- ggplotBoxPlot(data_in, colID = "Abundance", colors, title) + theme(legend.position = "none")
p1 <- p1 + theme(axis.text.x = element_blank())

# Generate density plot.
p2 <- ggplotDensity(data_in, colID = "Abundance", title) + theme(legend.position = "none")
colors <- c(rep("yellow", 11), rep("blue", 11), rep("green", 11), rep("purple", 11))
p2 <- p2 + scale_color_manual(values = colors)

# Generate meanSd plot.
p3 <- ggplotMeanSdPlot(data_in, colID = "Abundance", title, log = TRUE)

# Generate PCA plot.
colors <- rep(c("yellow","blue","green","purple"),each=11)
p4 <- ggplotPCA(data_in, traits=sample_info, colors, title = "2D PCA Plot") +
	theme(legend.position = "none")

# Store plots.
all_plots[[paste(tissue, "sl_bp", sep = "_")]] <- p1
all_plots[[paste(tissue, "sl_dp", sep = "_")]] <- p2
all_plots[[paste(tissue, "sl_msd", sep = "_")]] <- p3
all_plots[[paste(tissue, "sl_mds", sep = "_")]] <- p4

#-------------------------------------------------------------------------------
## Examine the nature of missing values.
#-------------------------------------------------------------------------------
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
p1 <- ggplotDetect(SL_peptide, groups[1]) #+ ggtitle(NULL)
l1 <- cowplot::get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

p2 <- ggplotDetect(SL_peptide, groups[2]) #+ ggtitle(NULL)
p2 <- p2 + theme(legend.position = "none")

p3 <- ggplotDetect(SL_peptide, groups[3]) #+ ggtitle(NULL)
p3 <- p3 + theme(legend.position = "none")

p4 <- ggplotDetect(SL_peptide, groups[4]) #+ ggtitle(NULL)
p4 <- p4 + theme(legend.position = "none")

# Store plots.
all_plots[[paste(tissue, "missing_val_shank2", sep = "_")]] <- p1
all_plots[[paste(tissue, "missing_val_shank3", sep = "_")]] <- p2
all_plots[[paste(tissue, "missing_val_syngap1", sep = "_")]] <- p3
all_plots[[paste(tissue, "missing_val_ube3a", sep = "_")]] <- p4

#-------------------------------------------------------------------------------
## Impute missing peptide values within an experiment.
#-------------------------------------------------------------------------------
# The function __impute_Peptides__ supports imputing missing values with the
# maximum likelyhood estimation (MLE) or KNN algorithms for missing not at
# random (MNAR) and missing at random data, respectively. Impution is performed
# with an experiment, and rows with more than 50% missing values are censored
# and will not be imputed. Peptides with more than 2 missing biological
# replicates or any missing quality control (QC) replicates will be
# censored and are not imputed.

# Define experimental groups for checking QC variability:
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")

# Impute missing values using KNN algorithm for MNAR data.
# Rows with missing QC replicates are ingored (qc_threshold=0).
# Rows with more than 2 (50%) missing biological replicates are
# ignored (bio_threshold=2).
data_impute <- impute_peptide(SL_peptide, groups, method = "knn")
imp_peptide <- data_impute$data_imputed

# Table of n imputed peptides.
n_out <- data_impute$n_out
mytable <- as.data.frame(do.call(rbind, n_out))
mytable <- tibble::add_column(mytable, rownames(mytable), .before = 1)
colnames(mytable) <- c("Experiment", "N Imputed")
mytable <- tableGrob(mytable, rows = NULL, theme = ttheme_default())

# Store table.
all_plots[[paste(tissue, "n_imputed_pep_tab", sep = "_")]] <- mytable

#-------------------------------------------------------------------------------
## Illustrate the mean variance relationship of QC peptides.
#-------------------------------------------------------------------------------
# Quality control samples can be used to asses intra-experimental variability.
# Peptides that have highly variable QC measurements will increase protein
# level variability and should be removed. The peptide-level QC data are
# binned by intensity, and peptides whose mean ratio are more than four
# standard deviations away from the mean are considered outliers and
# removed.

# Generate QC correlation scatter plots for all experimental groups.
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")
plots <- ggplotcorQC(imp_peptide, groups, colID = "QC", nbins = 5)

# Generate intensity bin histograms. 
# This will take a couple minutes.
hist_list <- lapply(as.list(groups),function(x) 
		    ggplotQCHist(imp_peptide,x, nbins = 5, threshold = 4))
names(hist_list) <- groups

# Store plots.
all_plots[[paste(tissue, "corQC_list", sep = "_")]] <- plots
all_plots[[paste(tissue, "histQC_list", sep = "_")]] <- hist_list

#-------------------------------------------------------------------------------
## Peptide level filtering.
#-------------------------------------------------------------------------------
# Peptides that were not quantified in all three qc replicates are removed.
# The data are binned by intensity, and measruments that are 4xSD from the mean
# ratio of the intensity bin are considered outliers and removed.

# Define experimental groups for checking QC variability:
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")

# Filter peptides based on QC precision.
filt_peptide <- filter_QC(imp_peptide, groups, nbins = 5, threshold = 4)

# Generate table.
out <- list(c(94, 77, 133, 59), c(182, 67, 73, 75))[[type]] # Cox and Str peps removed.
mytable <- data.frame(cbind(groups, out))
mytable <- tableGrob(mytable, rows = NULL, theme = ttheme_default())

# Store table
all_plots[[paste(tissue, "n_pep_filtered_cortex", sep = "_")]] <- mytable

#-------------------------------------------------------------------------------
##  Protein level summarization and normalization across all batches.
#-------------------------------------------------------------------------------
# Summarize to protein level by summing peptide intensities. Note that the
# peptide column in the returned data frame reflects the total number of
# peptides identified for a given protein across all experiments.

# Summarize to protein level:
filt_protein <- summarize_protein(filt_peptide)

# Normalize across all columns (experiments).
SL_protein <- normalize_SL(filt_protein, "Abundance", "Abundance")

#-------------------------------------------------------------------------------
## IntraBatch Protein-lavel ComBat.
#-------------------------------------------------------------------------------
# Each experimental cohort of 8 was prepared in two batches. This was necessary
# because the ultra-centrifuge rotor used to prepare purified synaptosomes
# holds a maximum of 6 samples. This intra-batch batch effect was recorded for
# 6/8 experiments. Here I will utilize the __ComBat()__ function from the `sva`
# package to remove this batch effect before correcting for inter-batch batch
# effects between batches with IRS normalization and regression. Note that in
# the absence of evidence of a batch effect (not annotated or cor(PCA,batch)<0.1),
# ComBat is not applied.

# Define experimental groups and column ID for expression data.
colID <- "Abundance"
data_in <- SL_protein

# Loop to perform ComBat on intraBatch batch effect (prep date).
# If there is no known batch effect (bicor(batch,PC1)<0.1) then 
# the data is returned un-regressed.

# Note: The values of QC samples are not adjusted by ComBat. The QC 
# samples were prepared from a seperate batch of mice and represent a single batch.

data_out <- list() # ComBat data.
plot_list <- list() # MDS plots.
R <- list() # Bicor stats [bicor(batch,PC1)]

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
    print("Warning: Expression values less than 1 will be replaced with 1.")
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
    print("Warning: Names of traits and expression data do not match.")
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
  cat(paste("Performing", groups[i], "ComBat...", "\n"))
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
      cat(c(
        "Error: No quantifiable batch effect!",
        "\n", "The un-regressed data will be returned.", "\n"
      ))
      data_ComBat <- log2(data)
    } else {
      # No batch effect.
      cat(c(
        "Error: ComBat can only be applied to factors with more than two levels!",
        "\n", "The un-regressed data will be returned.", "\n"
      ))
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
data_return <- data_out %>% purrr::reduce(left_join, by = c(colnames(data_in)[c(1, 2)]))

# Quantifying the batch effect.
# Check bicor correlation with batch before and after ComBat.
df <- do.call(rbind, lapply(R, function(x) x[1, ]))
df <- as.data.frame(t(apply(df, 1, function(x) round(abs(as.numeric(x)), 3))))
rownames(df) <- groups
df <- tibble::add_column(df, rownames(df), .before = 1)
colnames(df) <- c("Experiment", "preComBat", "postComBat")
mytable <- tableGrob(df, rows = NULL)

# Store plots...
all_plots[[paste(tissue, "batch_effect_tab", sep = "_")]] <- mytable
all_plots[[paste(tissue, "shank2_combat_pca", sep = "_")]] <- plot_list[[1]]
all_plots[[paste(tissue, "shank3_combat_pca", sep = "_")]] <- plot_list[[2]]
all_plots[[paste(tissue, "syngap1_combat_pca", sep = "_")]] <- plot_list[[3]]
all_plots[[paste(tissue, "ube3a_combat_pca", sep = "_")]] <- plot_list[[4]]

#-------------------------------------------------------------------------------
##  Examine protein identification overlap.
#-------------------------------------------------------------------------------
# Approximately 80-90% of all proteins are identified in all experiments.

# Inspect the overlap in protein identifcation.
plot <- ggplotFreqOverlap(SL_protein, "Abundance", groups) +
  ggtitle("Protein Identification Overlap")

# Store plot.
all_plots[[paste(tissue, "prot_id_overlap", sep = "_")]] <- plot

#-------------------------------------------------------------------------------
## Examine the Normalized protein level data.
#-------------------------------------------------------------------------------

# Generate boxplot.
data_in <- SL_protein
title <- "Normalized protein"
colors <- c(rep("green", 11), rep("purple", 11), rep("yellow", 11), rep("blue", 11))
p1 <- ggplotBoxPlot(data_in, colID = "Abundance", colors, title)

# Generate density plot.
p2 <- ggplotDensity(data_in, colID = "Abundance", title) + theme(legend.position = "none")
colors <- c(rep("yellow", 11), rep("blue", 11), rep("green", 11), rep("purple", 11))
p2 <- p2 + scale_color_manual(values = colors)

# Generate meanSd plot.
p3 <- ggplotMeanSdPlot(data_in, colID = "Abundance", title, log = TRUE)

# Generate PCA plot.
colors <- rep(c("yellow","blue","green","purple"),each=11)
p4 <- ggplotPCA(data_in, traits=sample_info, colors, title = "2D PCA Plot") +
	theme(legend.position = "none")

# Store plots.
all_plots[[paste(tissue, "sl_prot_bp", sep = "_")]] <- p1
all_plots[[paste(tissue, "sl_prot_dp", sep = "_")]] <- p2
all_plots[[paste(tissue, "sl_prot_msd", sep = "_")]] <- p3
all_plots[[paste(tissue, "sl_prot_mds", sep = "_")]] <- p4

#-------------------------------------------------------------------------------
## IRS Normalization.
#-------------------------------------------------------------------------------
# Internal reference sclaing (IRS) normalization equalizes the protein-wise means
# of reference (QC) samples across all batches. Thus, IRS normalization accounts
# for the random sampling of peptides at the MS2 level which results in the
# identification/quantificaiton of proteins by different peptides in each
# experiment. IRS normalization was first described by __Plubell et al., 2017__.

# Perform IRS normaliztion.
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")
IRS_protein <- normalize_IRS(SL_protein, "QC", groups, robust = TRUE)

#-------------------------------------------------------------------------------
## Identify and remove QC outliers.
#-------------------------------------------------------------------------------
# IRS normalization utilizes QC samples as reference samples. Outlier QC
# measurements (caused by interference or other artifact) would influence the
# create unwanted variability. Thus, outlier QC samples are removed, if
# identified. The method used by __Oldham et al., 2016__ is used to identify
# QC sample outliers. A threshold of -2.5 is used.

# Data is...
data_in <- IRS_protein

# Illustrate Oldham's sample connectivity.
sample_connectivity <- ggplotSampleConnectivity(IRS_protein,
  colID = "QC",
  threshold = -2.5
)
tab <- sample_connectivity$table
df <- tibble::add_column(tab, SampleName = rownames(tab), .before = 1)
rownames(df) <- NULL

plot1 <- sample_connectivity$connectivityplot +
  ggtitle("QC Sample Connectivity")

# Loop to identify Sample outliers using Oldham's connectivity method.
n_iter <- 5
threshold <- -2.5
data_in <- SL_protein
out_samples <- list()
plots <- list()

# Loop:
for (i in 1:n_iter) {
  data_temp <- quiet(normalize_IRS(data_in, "QC", groups, robust = TRUE))
  oldham <- ggplotSampleConnectivity(data_temp, log = TRUE, colID = "QC")
  plots[[i]] <- oldham$connectivityplot +
    ggtitle(paste("Sample Connectivity (Iteration = ", i, ")", sep = ""))
  bad_samples <- rownames(oldham$table)[oldham$table$Z.Ki < threshold]
  print(paste(
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
bad_samples

# Remove outliers from data.
samples_out <- paste(bad_samples, collapse = "|")
out <- grepl(samples_out, colnames(SL_protein))

# Redo IRS after outlier removal..
IRS_OutRemoved_protein <- normalize_IRS(SL_protein[, !out], "QC", groups, robust = TRUE)

# Write over IRS_data
IRS_protein <- IRS_OutRemoved_protein

# Store plots.
all_plots[[paste(tissue, "sample_connectivity_list", sep = "_")]] <- plots

#-------------------------------------------------------------------------------
## Examine the IRS Normalized protein level data.
#-------------------------------------------------------------------------------

# Generate boxplot.
data_in <- IRS_protein
title <- "IRS Normalized protein"
colors <- c(
  rep("green", 11), rep("purple", 11),
  rep("yellow", 11), rep("blue", 11)
)
p1 <- ggplotBoxPlot(data_in, colID = "Abundance", colors, title)

# Generate density plot.
p2 <- ggplotDensity(data_in, colID = "Abundance", title) +
  theme(legend.position = "none")
colors <- c(
  rep("yellow", 11), rep("blue", 11),
  rep("green", 11), rep("purple", 11)
)
p2 <- p2 + scale_color_manual(values = colors)

# Generate meanSd plot.
p3 <- ggplotMeanSdPlot(data_in, colID = "Abundance", title, log = TRUE)

# Generate PCA plot.
colors <- rep(c("yellow","blue","green","purple"),each=11)
p4 <- ggplotPCA(data_in, traits=sample_info, colors, title = "2D PCA Plot") +
	theme(legend.position = "none")

# Store plots.
all_plots[[paste(tissue, "irs_bp", sep = "_")]] <- p1
all_plots[[paste(tissue, "irs_dp", sep = "_")]] <- p2
all_plots[[paste(tissue, "irs_msd", sep = "_")]] <- p3
all_plots[[paste(tissue, "irs_mds", sep = "_")]] <- p4

#-------------------------------------------------------------------------------
## Protein level filtering, imputing, and final TMM normalization.
#-------------------------------------------------------------------------------
# Proteins that are identified by only a single peptide are removed. Proteins
# that are identified in less than 50% of all samples are also removed. The
# nature of the remaining missng values are examined by density plot and
# imputed with the KNN algorithm for MNAR data. Finally, TMM normalization is
# applied to correct for any biases introduced by these previous steps.

# Remove proteins that are identified by only 1 peptide as well as
# proteins identified in less than 50% of samples.
filt_protein <- filter_proteins(IRS_protein, "Abundance")

# Generate plot to examine distribution of remaining missing values.
plot <- ggplotDetect(filt_protein, "Abundance") +
  ggtitle("Protein missing value distribution")

# Impute the remaining number of missing values with KNN.
imp_protein <- impute_protein(filt_protein, "Abundance", method="knn")

# Final normalization with TMM.
TMM_protein <- normalize_TMM(imp_protein, "Abundance")

# Store plots.
all_plots[[paste(tissue, "prot_missing_val_dp", sep = "_")]] <- plot

#-------------------------------------------------------------------------------
## Examine the TMM Normalized protein level data.
#-------------------------------------------------------------------------------

# Generate boxplot.
# Adjust color vector if samples were removed.
# Cortex outliers = 1x Ube3a, and 1x Syngap1
data_in <- TMM_protein
title <- "TMM Normalized protein"
colors <- c(rep("green", 11), rep("purple", 11), rep("yellow", 11), rep("blue", 11))
p1 <- ggplotBoxPlot(data_in, colID = "Abundance", colors, title)

# Generate density plot.
p2 <- ggplotDensity(data_in, colID = "Abundance", title) + theme(legend.position = "none")
colors <- c(rep("yellow", 11), rep("blue", 11), rep("green", 11), rep("purple", 11))
p2 <- p2 + scale_color_manual(values = colors)

# Generate meanSd plot.
p3 <- ggplotMeanSdPlot(data_in, colID = "Abundance", title, log = TRUE)

# Generate MDS plot.
colors <- c(rep("yellow", 3), rep("blue", 3), rep("green", 3), rep("purple", 3))
p4 <- ggplotMDS(data_in, colID = "Abundance", colors, title, sample_info, labels = TRUE) +
  theme(legend.position = "none")

# Store plots.
all_plots[[paste(tissue, "tmm_bp", sep = "_")]] <- p1
all_plots[[paste(tissue, "tmm_dp", sep = "_")]] <- p2
all_plots[[paste(tissue, "tmm_msd", sep = "_")]] <- p3
all_plots[[paste(tissue, "tmm_mds", sep = "_")]] <- p4

#-------------------------------------------------------------------------------
## Reformat final normalized, regressed data for TAMPOR Normalization.
#-------------------------------------------------------------------------------
#' Data are reformatted for TAMPOR normalization in the `2_TMT_Analysis.R` script.

# Data is...
data_in <- TMM_protein
rownames(data_in) <- TMM_protein$Accession
data_in$Accession <- NULL
data_in$Peptides <- NULL
data_in <- as.matrix(data_in)

# data_in <- 2^t(data.fit)
# data_in[1:5, 1:5] # un-log
# dim(data_in)

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

# Save as cleanDat
file <- paste0(Rdatadir, "/", outputMatName, "_CleanDat_TAMPOR_Format.Rds")
saveRDS(cleanDat, file)

# Save plot list.
myfile <- file.path(Rdatadir, paste0(outputMatName, "_plots.Rds"))
saveRDS(all_plots, myfile)

## ENDOFILE
#------------------------------------------------------------------------------
