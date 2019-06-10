#' ---
#' title: TMT Analysis part 1.
#' author: Tyler W Bradshaw
#' urlcolor: blue
#' header-includes:
#' - \usepackage{float}
#' - \floatplacement{figure}{H}
#' output:
#'    pdf_document:
#'      fig_caption: true
#'      toc: true
#'      number_sections: false
#'      highlight: tango
#' ---
#' 
#-------------------------------------------------------------------------------
#' ## Prepare the workspace.
#-------------------------------------------------------------------------------
#' Prepare the R workspace for the analysis. Load all required packages and 
#' custom functions. Set-up directories for saving output files. 
#'
#+ eval = TRUE, echo = FALSE, error = FALSE, results = 'hide'
 

# Use ctl+alt+T to execute an entire code chunk.
   
# Run this chunk before doing anything!
rm(list = ls())
if(!is.null(dev.list())) {
  dev.off()
  }
cat("\014") # alternative is cat("\f")
options(stringsAsFactors = FALSE)

# Sometimes, if you have not cleared the workspace of all loaded packages,
# you man incounter problems.
# To remove all packages, you can call the following:
library(magrittr)
library(JGmisc)
detachAllPackages(keep = NULL)

#  Load required packages.
suppressWarnings({
  suppressPackageStartupMessages({
    library(JGmisc)
    library(readxl)
    library(knitr)
    library(readr)
    library(dplyr)
    library(reshape2)
    library(DEP)
    library(tibble)
    library(SummarizedExperiment)
    library(ggplot2)
    library(hexbin)
    library(vsn)
    library(BurStMisc)
    library(dplyr)
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    library(edgeR)
    library(openxlsx)
    library(stringr)
    library(imp4p)
    library(Cairo)
    library(pryr)
    library(qvalue)
    library(gridExtra)
    library(cowplot)
    library(WGCNA)
    library(impute)
    library(ggrepel)
    library(sva)
    library(anRichment)
    library(ggdendro)
    library(flashClust)
    library(purrr)
    library(ggpubr)
    library(doParallel)
    library(NMF)
    library(FSA)
    library(plyr)
    library(RColorBrewer)
    library(gtable)
    library(grid)
    library(ggplotify)
    library(TBmiscr)
  })
})

# To install TBmiscr:
#library(devtools)
#devtools::install_github("twesleyb/TBmiscr")

# Define version of the code.
CodeVersion <- "Final_TMT_Analysis_part1"

# Define tisue type: cortex = 1; striatum = 2.
type <- 2
tissue <- c("Cortex", "Striatum")[type]

# Set the working directory.
rootdir <- "D:/Documents/R/Synaptopathy-Proteomics"
setwd(rootdir)

# Set any other directories.
functiondir <- paste(rootdir, "Code", sep = "/")
datadir <- paste(rootdir, "Input", sep = "/")
Rdatadir <- paste(rootdir, "RData", sep = "/")

# Create code-version specific figure and table folders if they do not already 
# exist.
# Creat otuput directory for figures.
outputfigs <- paste(rootdir, "Figures", tissue, sep = "/")
outputfigsdir <- paste(outputfigs, CodeVersion, sep = "/")
if (!file.exists(outputfigsdir)) {
  dir.create(file.path(outputfigsdir))
} else {
  msg <- c("This directory already exists.",
           "Warning: Some files may be overwritten when running this script.")
  print(msg)
}

# Create output directory for tables.
outputtabs <- paste(rootdir, "Tables", tissue, sep = "/")
outputtabsdir <- paste(outputtabs, CodeVersion, sep = "/")
if (!file.exists(outputtabsdir)) {
  dir.create(file.path(outputtabsdir))
} else {
  print(msg)
}
# Create output directory for reports.
outputreports <- paste(rootdir, "Reports", tissue, sep = "/")
outputrepsdir <- paste(outputreports, CodeVersion, sep = "/")
if (!file.exists(outputrepsdir)) {
  dir.create(file.path(outputrepsdir))
} else {
  print(msg)
}

# Load required custom functions.
my_functions <- paste(functiondir, "0_TMT_Preprocess_Functions.R", sep = "/")
source(my_functions)

# Define prefix for output figures and tables.
outputMatName <- paste(tissue, "_TMT_Analysis", sep = "")

# Save script.
rstudioapi::documentSave()

# Globally set ggplots theme.
ggplot2::theme_set(theme_gray())

#-------------------------------------------------------------------------------
#' ## Load the raw data and sample info (traits) from excel.
#-------------------------------------------------------------------------------
#' The raw peptide intensity data were exported from ProteomeDiscover (PD) 
#' version 2.2. Note that the default export from PD2.x is a unitless signal to
#' noise ratio, and it is not recommended to use ths for quantification. 
#' 
#+ eval = TRUE, message = FALSE

# Load the data from excel using readxl::read_excel
datafile <- c(
  "4227_TMT_Cortex_Combined_PD_Peptide_Intensity.xlsx",
  "4227_TMT_Striatum_Combined_PD_Peptide_Intensity.xlsx"
)

# Load sample information.
samplefile <- c(
  "4227_TMT_Cortex_Combined_PD_Protein_Intensity_EBD_traits.csv",
  "4227_TMT_Striatum_Combined_PD_Protein_Intensity_EBD_traits.csv"
)

# Load the data from PD and sample info.
data_PD <- read_excel(paste(datadir, "/", datafile[type], sep = ""), 1)
sample_info <- read.csv(paste(datadir, "/", samplefile[type], sep = ""))

# Insure traits are in matching order.
sample_info <- sample_info[order(sample_info$Order), ]

#-------------------------------------------------------------------------------
#' ## Cleanup and reorganize the data from PD.
#-------------------------------------------------------------------------------
#' The raw data are cleaned up with the custom function __cleanPD__.
#' 
#+ eval = TRUE

# Load the data from PD.
raw_peptide <- cleanPD(data_PD, sample_info)

#-------------------------------------------------------------------------------
#' ## Examine an example peptide.
#-------------------------------------------------------------------------------

dat <- subset(raw_peptide,grepl("Dlg4",raw_peptide$Description))
rownames(dat) <- paste(dat$Accession,dat$Sequence,c(1:nrow(dat)),sep="_")
idy <- grepl("Shank2",colnames(dat))
dat <- dat[,idy]
dat <- na.omit(dat)

# Make bar plot for given peptide.
colIDs <- gsub(",","",sapply(strsplit(colnames(dat),"\\ "),"[",3))

geno <- gsub(",","",sapply(strsplit(colnames(dat),"\\ "),"[",5))


n <- sample(nrow(dat),1)
df <- melt(dat[n,])
df$Channel <- colIDs
title <- strsplit(rownames(dat)[n],"_")[[1]][2]

plot <- ggplot(df, aes(x = Channel, y = value, fill = Channel)) + 
  geom_bar(stat="identity", width = 0.9, position = position_dodge(width = 1)) + 
  xlab("TMT Channel") + ylab("Intensity") + 
  ggtitle(title) + 
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(color = "black", size = 11, face = "bold"))
plot

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_Example_TMT.tiff")
ggsave(file,plot, width = 3, height = 2.5, units = "in")

#-------------------------------------------------------------------------------
#' ## Examine peptide and protein level identification overalap.
#-------------------------------------------------------------------------------
#' Approximately 40,000 unique peptides cooresponding to ~3,000 proteins are 
#' quantified across all four experiments. When comparing the peptides 
#' identified in each experiment, the overlap is only ~25%. ~20$ of peptides are
#' identified in all four experiments. This means that in different experiments,
#' the same protein will likely be quantified by different peptides. 
#' This is a major reason the internal reference scaling (IRS) normalizaiton 
#' approach is employed below.
#' 
#+ eval = TRUE

# Total number of unique peptides:
nPeptides <- format(length(unique(raw_peptide$Sequence)), big.mark = ",")
print(paste(nPeptides, " unique peptides identified.", sep = ""))

# Total number of unique proteins
nProteins <- format(length(unique(raw_peptide$Accession)), big.mark = ",")
print(paste(nProteins, " unique proteins identified.",sep = ""))

# Table
mytable <- data.frame(nPeptides = nPeptides,
                         nProteins = nProteins)
table <- tableGrob(mytable, rows = NULL, theme = ttheme_default())
grid.arrange(table)

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_Raw_nPeptides_nProteins.tiff")
ggsave(file,table)

# Examine the number of peptides per protein.
nPep <- subset(raw_peptide) %>%
  group_by(Accession) %>%
  dplyr::summarize(nPeptides = length(Sequence))

# Creat table with stats.
stats <- as.matrix(summary(nPep$nPeptides))
df <- add_column(as.data.frame(stats), rownames(stats), .before = 1)
tt <- ttheme_default(base_size = 11, core = list(bg_params = list(fill = "white")))
tab <- tableGrob(df, cols = NULL, rows = NULL, theme = tt)
g <- gtable_add_grob(tab,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
  t = 1, b = nrow(tab), l = 1, r = ncol(tab)
)

# Generate plot.
plot <- ggplot(nPep, aes(nPeptides)) + geom_histogram(bins = 100, fill = "black") +
  ggtitle("Number of peptides per protein") +
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold")
  )

# Add annotation.
p <- ggranges(plot)$TopRight # ggranges calculates position of annotation.
plot <- plot + annotation_custom(g, p$xmin, p$xmax, p$ymin, p$ymax)

#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = "Number of peptides identified per protein."
plot

# Save figures.
#plots <- list(plot)
#file <- paste0(outputfigsdir, "/", outputMatName, "_nPeptide_per_Protein.pdf")
#ggsavePDF(plots, file)

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_nPeptide_per_Protein.tiff")
ggsave(file,plot)

# Peptide identification overlap per pairwise comparisons of experiments.
contrasts <- combn(c("Shank2", "Shank3", "Syngap1", "Ube3a"), 2)
info_cols <- c(1, 2, 3, 4, 5)
overlap <- peptide_overlap_TMT(raw_peptide, contrasts, info_cols)
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
table <- tableGrob(mytable, rows = NULL, theme = ttheme_default())

#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = "Peptide identification overlap. All experimental pairwise comparisons."
grid.arrange(table)

# Plot peptide identification overlap.
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")
plot <- ggplotFreqOverlap(raw_peptide, "Abundance", groups) + ggtitle("Peptide Identification Overlap")

#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = "Peptide identification overlap."
plot

# Save figures.
#plots <- list(plot, table)
#file <- paste0(outputfigsdir, "/", outputMatName, "_Peptide_identification_overlap.pdf")
#ggsavePDF(plots, file)

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_Peptide_identification_overlap_Barplot.tiff")
ggsave(file,plot)

file <- paste0(outputfigsdir, "/", outputMatName, "_Peptide_identification_overlap_Table.tiff")
ggsave(file,table)

#-------------------------------------------------------------------------------
#' ## Examine the raw data.
#-------------------------------------------------------------------------------
#' The need for normalization is evident in the raw data. Note that in the MDS
#' plot, samples cluster by experiment--evidence of a batch effect.
#'   
#+ eval = TRUE

data_in <- raw_peptide
title <- NULL
# Colors for boxplot must be specified in ggplot order for boxplot.
colors <- c(rep("green", 11), rep("purple", 11), rep("yellow", 11), rep("blue", 11))

# Generate boxplot.
p1 <- ggplotBoxPlot(data_in, colID = "Abundance", colors, title)

# Generate density plot.
p2 <- ggplotDensity(data_in, colID = "Abundance", title) + theme(legend.position = "none")
# Genotype specific colors must be specified in column order.
colors <- c(rep("yellow", 11), rep("blue", 11), rep("green", 11), rep("purple", 11))
p2 <- p2 + scale_color_manual(values = colors)

# Generate meanSd plot.
p3 <- ggplotMeanSdPlot(data_in, colID = "Abundance", title, log = TRUE)

# Generate MDS plot.
colors <- c(rep("yellow", 3), rep("blue", 3), rep("green", 3), rep("purple", 3))
p4 <- ggplotMDS(data_in, colID = "Abundance", colors, title, sample_info, labels = TRUE) +
  theme(legend.position = "none")

# Figure.
caption <- strwrap("Raw peptide data. A. Boxplot B. Density plot. 
                   C. Mean SD plot. D. MDS plot.", width = Inf, simplify = TRUE)
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = caption
fig <- plot_grid(p1, p2, p3, p4, labels = "auto")
fig 

# Save plots.
#plot_list <- list(p1, p2, p3, p4)
#file <- paste0(outputfigsdir, "/", outputMatName, "_Raw_peptide_plots.pdf")
#ggsavePDF(plot_list, file)

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_Raw_Peptide.tiff")
ggsave(file,fig)

#-------------------------------------------------------------------------------
#' ## Sample loading normalization within experiments.
#-------------------------------------------------------------------------------
#' The function __normalize_SL__ performs sample loading (SL) normalization to
#' equalize the run level intensity (column) sums. The data in each column are 
#' multiplied by a factor such that the mean of the column sums are are equal.
#' Sample loading normalization is performed within an experiment under the 
#' assumption that equal amounts of protein were used for each of the 11 TMT
#' channels.
#' 
#+ eval = TRUE

# Define data columns for SL within experiments:
colID <- "Abundance"
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")

# Perform SL normalization.
SL_peptide <- normalize_SL(raw_peptide, colID, groups)

#-------------------------------------------------------------------------------
#' ## Examine the SL Data.
#-------------------------------------------------------------------------------
#' Sample loading normalization equalizes the run-level sums within 
#' an 11-plex TMT experiment. 
#' 
#+ eval = TRUE

data_in <- SL_peptide
title <- NULL

# Generate boxplot.
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

# Figure.
caption <- strwrap("Sample loading normalized peptide data. A. Boxplot B. Density plot. 
                   C. Mean SD plot. D. MDS plot.", width = Inf, simplify = TRUE)
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = caption
fig <- plot_grid(p1, p2, p3, p4, labels = "auto")
fig 

# Save plots.
#plot_list <- list(p1, p2, p3, p4)
#file <- paste0(outputfigsdir, "/", outputMatName, "_SL_Normalized_Peptide.pdf")
#ggsavePDF(plot_list, file)

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_SL_Peptide.tiff")
ggsave(file,fig)

#-------------------------------------------------------------------------------
#' ## Examine the nature of missing values.
#-------------------------------------------------------------------------------
#' Missing values are inherent in high throughput experiments. There are two 
#' main classes of missing values, missing at random (MAR) and missing not at
#' random (MNAR). The appropriate imputing algorithm should be chosen based on 
#' the nature of missing values. The distribution of peptides with missing 
#' values is examined by density plots. The left-shifted distribution of 
#' peptides with missing values indicates that peptides that have missing 
#' values are generally lower in abundance. Missing values are likely then to be
#' missing not at random (MNAR), but missing because they are low-abundance 
#' and at or near the limit of detection. MNAR data can be imputed with the 
#' k-nearest neighbors (knn) algorithm in the next chunk. The __impute.knn__
#' from the package `impute` is used to impute MNAR data.
#' 
#+ eval = TRUE

# Define groups for subseting the data.
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")

# Generate plots.
p1 <- ggplotDetect(SL_peptide, groups[1]) + ggtitle(NULL)
p2 <- ggplotDetect(SL_peptide, groups[2]) + ggtitle(NULL)
p3 <- ggplotDetect(SL_peptide, groups[3]) + ggtitle(NULL)
p4 <- ggplotDetect(SL_peptide, groups[4]) + ggtitle(NULL)

# Figure.
caption <- strwrap("Peptide-level missing value distributions. A. Shank2 B. Shank3. 
                   C. Syngap1. D. Ube3a.", width = Inf, simplify = TRUE)
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = caption
fig <- plot_grid(p1, p2, p3, p4, labels = "auto")
fig

# Save plots.
#plot_list <- list(p1, p2, p3, p4)
#file <- paste0(outputfigsdir, "/", outputMatName, "_Peptide_Missing_Value_Density_plots.pdf")
#ggsavePDF(plot_list, file)

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_Peptide_Missing_Values.tiff")
ggsave(file,fig)

#-------------------------------------------------------------------------------
#' ## Impute missing peptide values within an experiment.
#-------------------------------------------------------------------------------
#' The function __impute_Peptides__ supports imputing missing values with the 
#' maximum likelyhood estimation (MLE) or KNN algorithms for missing not at 
#' random (MNAR) and missing at random data, respectively. Impution is performed 
#' with an experiment, and rows with more than 50% missing values are censored 
#' and will not be imputed. Peptides with more than 2 missing biological 
#' replicates or any missing quality control (QC) replicates will be 
#' censored and are not imputed. 
#' 
#+ eval = TRUE

# Define experimental groups for checking QC variability:
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")

# Impute missing values using KNN algorithm for MNAR data.
# Rows with missing QC replicates are ingored (qc_threshold=0).
# Rows with more than 2 (50%) missing biological replicates are 
# ignored (bio_threshold=2).
data_impute <- impute_Peptides(SL_peptide, groups, method = "knn")
impute_peptide <- data_impute$data_imputed

# Table of n imputed peptides.
n_out <- data_impute$n_out
table <- as.data.frame(do.call(rbind, n_out))
table <- add_column(table, rownames(table), .before = 1)
colnames(table) <- c("Experiment", "N Imputed")
table <- tableGrob(table, rows = NULL, theme = ttheme_default())

# Table.
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = "Summary of imputed values." 
grid.arrange(table)

# Save figures.
#file <- paste0(outputfigsdir, "/", outputMatName, "_N_Imputed_Peptides.pdf")
#ggsavePDF(table, file)

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_N_Imputed_Peptides.tiff")
ggsave(file,table)

#-------------------------------------------------------------------------------
#' ## Illustrate the mean variance relationship of QC peptides.
#-------------------------------------------------------------------------------
#' Quality control samples can be used to asses intra-experimental variability.
#' Peptides that have highly variable QC measurements will increase protein
#' level variability and should be removed. The peptide-level QC data are 
#' binned by intensity, and peptides whose mean ratio are more than four
#' standard deviations away from the mean are considered outliers and
#' removed. 
#' 
#+ eval = TRUE, warning = FALSE

# Generate QC correlation scatter plots for all experimental groups.
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")
plots <- ggplotcorQC(impute_peptide, groups, colID = "QC", nbins = 5)

# Generate intensity bin histograms. Example, Shank2.
hist_list <- list()
hist_list[["Shank2"]] <- ggplotQCHist(impute_peptide, "Shank2", nbins = 5, threshold = 4)
hist_list[["Shank3"]] <- ggplotQCHist(impute_peptide, "Shank3", nbins = 5, threshold = 4)
hist_list[["Syngap1"]] <- ggplotQCHist(impute_peptide, "Syngap1", nbins = 5, threshold = 4)
hist_list[["Ube3a"]] <- ggplotQCHist(impute_peptide, "Ube3a", nbins = 5, threshold = 4)

# Figure.
#fixme: should optimize size of figure for faster render on pdf. also, annotation layers are not scaled correctly. 
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = "Shank2."
p1 <- plots$Shank2 + theme(legend.position = "none")
p1
p2 <- hist_list$shank2[[1]]
p2

# Save all scatter plots.
#file <- paste0(outputfigsdir, "/", outputMatName, "_QC_ScatterPlots.pdf")
#ggsavePDF(plots, file)

# Save histograms.
#file <- paste0(outputfigsdir, "/", outputMatName, "_QC_Ratio_Histograms.pdf")
#ggsavePDF(hist_list, file)

# Save tiffs. 
genos <- c("Shank2","Shank3","Syngap1","Ube3a")
file <- as.list(
  paste0(outputfigsdir, "/", outputMatName, 
               "_",genos,"_QC_ScatterPlot.tiff"))
# Use mapply to save plots.
quiet(mapply(ggsave, file, plots))

# Histograms.
x <- c(1:5)
file <- paste0(outputfigsdir, "/", outputMatName, "_Shank2","_QC_hist",x,".tiff")
quiet(mapply(ggsave, file, hist_list$Shank2))

#-------------------------------------------------------------------------------
#' ## Peptide level filtering.
#-------------------------------------------------------------------------------
#' Peptides that were not quantified in all three qc replicates are removed.
#' The data are binned by intensity, and measruments that are 4xSD from the mean
#' ratio of the intensity bin are considered outliers and removed.
#'
#+ eval = TRUE

# Define experimental groups for checking QC variability:
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")

# Filter peptides based on QC precision.
filter_peptide <- filterQCv2(impute_peptide, groups, nbins = 5, threshold = 4)

# Generate table.
out <- list(c(94,77,133,59),c(182,67,73,75))[[type]] # Cox and Str peps removed.
mytable <- data.frame(cbind(groups,out))
table <- tableGrob(mytable, rows = NULL, theme = ttheme_default())
grid.arrange(table)

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_nProteins_Filtered.tiff")
ggsave(file,table)

#-------------------------------------------------------------------------------
#' ##  Protein level summarization and normalization across all batches.
#-------------------------------------------------------------------------------
#' Summarize to protein level by summing peptide intensities. Note that the
#' peptide column in the returned data frame reflects the total number of
#' peptides identified for a given protein across all experiments. 
#'
#+ eval = TRUE

# Summarize to protein level:
SL_protein <- summarize_Protein(filter_peptide)

# Normalize across all columns (experiments).
SL_protein <- normalize_SL(SL_protein, "Abundance", "Abundance")

#-------------------------------------------------------------------------------
#' ## IntraBatch Protein-lavel ComBat.
#-------------------------------------------------------------------------------
#' Each experimental cohort of 8 was prepared in two batches. This was necessary
#' because the ultra-centrifuge rotor used to prepare purified synaptosomes
#' holds a maximum of 6 samples. This intra-batch batch effect was recorded for 
#' 6/8 experiments. Here I will utilize the __ComBat()__ function from the `sva` 
#' package to remove this batch effect before correcting for inter-batch batch 
#' effects between batches with IRS normalization and regression. Note that in 
#' the absence of evidence of a batch effect (not annotated or cor(PCA,batch)<0.1), 
#' ComBat is not applied.
#'
#+ eval = TRUE, message = FALSE

# Define experimental groups and column ID for expression data.
groups
colID <- "Abundance"
data_in <- SL_protein

# Loop to perform ComBat on intraBatch batch effect (prep date).
# If there is no known batch effect, the data is returned un-regressed.
# Note: QC samples are not adjusted by ComBat.

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
  r1 <- suppressWarnings(bicor(batch, pc1))

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
                                   as.numeric(as.factor(traits_sub$PrepDate)), sep = ".")
  title <- paste(gsub(" ", "", unique(traits_sub$Model)), "pre-ComBat", sep = " ")
  plot1 <- ggplotMDSv2(log2(data),
    colID = "b",
    title = title, traits = traits_sub
  )$plot + theme(legend.position = "none")
  plot1 <- plot1 + scale_color_manual(values = unique(traits_sub$Color))

  # Apply ComBat.
  cat(paste("Performing", groups[i], "ComBat...", "\n"))
  if (length(unique(CombatInfo$PrepDate)) > 1 & abs(r1) > 0.1) {
    # Create ComBat model.
    model <- model.matrix(~ as.factor(CombatInfo$SampleType), 
                          data = as.data.frame(log2(data)))
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
  r2 <- suppressWarnings(bicor(batch, pc1))

  # Check that r2 is not NA.
  if (is.na(r2)) {
    r2 <- 0
  }

  R[[i]] <- cbind(r1, r2)

  # Check MDS plot after ComBat.
  title <- paste(gsub(" ", "", unique(traits_sub$Model)), "post-ComBat", sep = " ")
  plot2 <- ggplotMDSv2(
    data_in = data_ComBat, colID = "Abundance",
    title = title, traits = traits_sub
  )$plot + theme(legend.position = "none")
  plot2 <- plot2 + scale_color_manual(values = unique(traits_sub$Color))
  # Add annotation layer.
  # xpos <- sum(unlist(ggplot_build(plot2)$layout$panel_params[[1]][1]))/2
  # ypos <- sum(unlist(ggplot_build(plot2)$layout$panel_params[[1]][8]))/2
  # lab <- paste("bicor(batch,PC1) = ",round(r2$bicor,3))
  # plot2 <- plot2 + annotate("text",  x = xpos, y = ypos, label = lab)

  # Un-log.
  data_ComBat <- 2^data_ComBat
  # Recombine with QC data.
  data_out[[i]] <- cbind(info_cols[!rows_out, ], data_QC, data_ComBat)
  names(data_out)[[i]] <- group

  plots <- list(plot1, plot2)
  plot_list[[i]] <- plots
  names(plot_list[[i]]) <- paste(groups[i], c("preComBat", "postComBat"))
}

# Merge the data frames with reduce()
data_return <- data_out %>% reduce(left_join, by = c(colnames(data_in)[c(1, 2)]))

# Quantifying the batch effect.
# Check bicor correlation with batch before and after ComBat.
df <- do.call(rbind, lapply(R, function(x) x[1, ]))
df <- as.data.frame(t(apply(df, 1, function(x) round(abs(as.numeric(x)), 3))))
rownames(df) <- groups
df <- add_column(df, rownames(df), .before = 1)
colnames(df) <- c("Experiment", "preComBat", "postComBat")
table <- tableGrob(df, rows = NULL)

# Table and figures.
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = "Quantification of the intra-batch batch effect."
grid.newpage()
grid.arrange(table)

# Save fig.
file <- paste0(outputfigsdir, "/", outputMatName, "_IntraBatchEffect_Quantificaiton.tiff")
ggsave(file,table)

# Plots to check MDS:
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = "Shank2 ComBat."
ggpubr::ggarrange(plotlist = plot_list[[1]])
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = "Shank3 ComBat."
ggpubr::ggarrange(plotlist = plot_list[[2]])
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = "Syngap1 ComBat."
ggpubr::ggarrange(plotlist = plot_list[[3]])
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = "Ube3a ComBat."
ggpubr::ggarrange(plotlist = plot_list[[4]])

# Save Plots.
#plot_list <- c(list(table), unlist(plot_list, recursive = FALSE))
#file <- paste0(outputfigsdir, "/", outputMatName, "_Protein_level_ComBat_PCA.pdf")
#ggsavePDF(plot_list, file)

# Write over SL_protein with regressed data.
#SL_protein <- data_return

#-------------------------------------------------------------------------------
#' ##  Examine protein identification overlap.
#-------------------------------------------------------------------------------
#' Approximately 80-90% of all proteins are identified in all experiments. 
#' 
#+ eval = TRUE

# Inspect the overlap in protein identifcation.
plot <- ggplotFreqOverlap(SL_protein, "Abundance", groups) +
  ggtitle("Protein Identification Overlap")

#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = "Protein identification overlap."
plot

# Save plot.
#file <- paste0(outputfigsdir, "/", outputMatName, "_Protein_identification_overlap.pdf")
#ggsavePDF(plot, file)

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_Protein_identification_overlap.tiff")
ggsave(file,plot)

#-------------------------------------------------------------------------------
#' ## Examine the Normalized protein level data.
#-------------------------------------------------------------------------------
#+ eval = TRUE

data_in <- SL_protein
title <- "Normalized protein"

# Generate boxplot.
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

# Figure.
caption <- strwrap("Normalized protein level data. A. Boxplot. B. Density plot. 
                   C. Mean SD plot. D. MDS plot.", width = Inf, simplify = TRUE)
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = caption
fig <- plot_grid(p1, p2, p3, p4, labels = "auto")
fig 

# Save plots.
#file <- paste0(outputfigsdir, "/", outputMatName, "_SL_Normalized_Protein.pdf")
#ggsavePDF(plots = list(p1, p2, p3, p4), file)

# Save tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_Normalized_Protein.tiff")
ggsave(file,fig)

#-------------------------------------------------------------------------------
#' ## IRS Normalization.
#-------------------------------------------------------------------------------
#' Internal reference sclaing (IRS) normalization equalizes the protein-wise means 
#' of reference (QC) samples across all batches. Thus, IRS normalization accounts 
#' for the random sampling of peptides at the MS2 level which results in the 
#' identification/quantificaiton of proteins by different peptides in each 
#' experiment. IRS normalization was first described by __Plubell et al., 2017__. 
#'
#+ eval = TRUE

# Perform IRS normaliztion.
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")
IRS_protein <- normalize_IRS(SL_protein, "QC", groups, robust = TRUE)

#-------------------------------------------------------------------------------
#' ## Identify and remove QC outliers.
#-------------------------------------------------------------------------------
#' IRS normalization utilizes QC samples as reference samples. Outlier QC 
#' measurements (caused by interference or other artifact) would influence the 
#' create unwanted variability. Thus, outlier QC samples are removed, if 
#' identified. The method used by __Oldham et al., 2016__ is used to identify 
#' QC sample outliers. A threshold of -2.5 is used. 
#' 
#+ eval = TRUE

# Data is...
data_in <- IRS_protein

# Illustrate Oldham's sample connectivity.
sample_connectivity <- ggplotSampleConnectivityv2(IRS_protein, colID = "QC", 
                                                  threshold = -2.5)
tab <- sample_connectivity$table
df <- add_column(tab,SampleName = rownames(tab),.before = 1)
rownames(df) <- NULL
knitr::kable(df)
plot <- sample_connectivity$connectivityplot + 
  ggtitle("QC Sample Connectivity")

# Figure.
caption <- strwrap("QC Sample Connectivity. Examination of QC samples Z-Score 
                  normalized connectivity as a means of identifying outlier 
                  QC samples.", width = Inf, simplify = TRUE)
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = caption
plot

# Save as tiff. 
file <- paste0(outputfigsdir, "/", outputMatName, "_QC_Sample_Outlier.pdf")
ggsave(file,plot)

# Loop to identify Sample outliers using Oldham's connectivity method.
n_iter <- 5
threshold <- -2.5
data_in <- SL_protein
out_samples <- list()
plots <- list()

# Loop:
for (i in 1:n_iter) {
  data_temp <- quiet(normalize_IRS(data_in, "QC", groups, robust = TRUE))
  oldham <- ggplotSampleConnectivityv2(data_temp, log = TRUE, colID = "QC")
  plots[[i]] <- oldham$connectivityplot +
    ggtitle(paste("Sample Connectivity (Iteration = ", i, ")", sep = ""))
  bad_samples <- rownames(oldham$table)[oldham$table$Z.Ki < threshold]
  print(paste(
    length(bad_samples), " outlier sample(s) identified in iteration ", i, ".", sep = "")
    )
  if (length(bad_samples) == 0) bad_samples <- "none"
  out_samples[[i]] <- bad_samples
  out <- grepl(paste(unlist(out_samples), collapse = "|"), colnames(data_in))
  data_in <- quiet(normalize_IRS(data_in[, !out], "QC", groups, robust = TRUE))
}

# Outlier samples.
bad_samples <- unlist(out_samples)
bad_samples

# Remove from data, re-do IRS, and check sample connectivity.
samples_out <- paste(bad_samples, collapse = "|")
out <- grepl(samples_out, colnames(SL_protein))

# Redo IRS after removal of outlier.
IRS_OutRemoved_protein <- normalize_IRS(SL_protein[, !out], "QC", groups, robust = TRUE)

# Illustrate Oldham's sample connectivity after outlier removal.
sample_connectivity <- ggplotSampleConnectivityv2(IRS_OutRemoved_protein, colID = "QC")
plot <- sample_connectivity$connectivityplot

# Figure.
caption <- c("QC Sample Connectivity after outlier removal.")
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = caption
plot

# Save plots
#file <- paste0(outputfigsdir, "/", outputMatName, "_QC_Sample_Connectivity.pdf")
#ggsavePDF(plots, file)

# Write over IRS_data
IRS_protein <- IRS_OutRemoved_protein

#-------------------------------------------------------------------------------
#' ## Examine the IRS Normalized protein level data.
#-------------------------------------------------------------------------------
#+ eval = TRUE

data_in <- IRS_protein
title <- "IRS Normalized protein"

# Generate boxplot.
colors <- c(rep("green", 11), rep("purple", 11), 
            rep("yellow", 11), rep("blue", 11))
p1 <- ggplotBoxPlot(data_in, colID = "Abundance", colors, title)

# Generate density plot.
p2 <- ggplotDensity(data_in, colID = "Abundance", title) + 
  theme(legend.position = "none")
colors <- c(rep("yellow", 11), rep("blue", 11), 
            rep("green", 11), rep("purple", 11))
p2 <- p2 + scale_color_manual(values = colors)

# Generate meanSd plot.
p3 <- ggplotMeanSdPlot(data_in, colID = "Abundance", title, log = TRUE)

# Generate MDS plot.
colors <- c(rep("yellow", 3), rep("blue", 3), rep("green", 3), rep("purple", 3))
p4 <- ggplotMDS(data_in, colID = "Abundance", colors, title, sample_info, labels = TRUE) +
  theme(legend.position = "none")

# Figure.
caption <- strwrap("IRS normalized protein. A. Boxplot. B. Density plot.
                   C. Mean SD plot. D. MDS plot", width = Inf, simplify = TRUE)
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = caption
fig <- plot_grid(p1, p2, p3, p4, labels = "auto")
fig 

# Save plots.
#file <- paste0(outputfigsdir, "/", outputMatName, "_IRS_Normalized_Protein.pdf")
#ggsavePDF(plots = list(p1, p2, p3, p4), file)

# Save tiff. 
file <- paste0(outputfigsdir, "/", outputMatName, "_IRS_Normalized_Protein.tiff")
ggsave(file,fig)

#-------------------------------------------------------------------------------
#' ## Protein level filtering, imputing, and final TMM normalization.
#-------------------------------------------------------------------------------
#' Proteins that are identified by only a single peptide are removed. Proteins 
#' that are identified in less than 50% of all samples are also removed. The 
#' nature of the remaining missng values are examined by density plot and 
#' imputed with the KNN algorithm for MNAR data. Finally, TMM normalization is
#' applied to correct for any biases introduced by these previous steps. 
#' 
#+ eval = TRUE

# Remove proteins that are identified by only 1 peptide.
# Remove proteins identified in less than 50% of samples.
filt_protein <- filter_proteins(IRS_protein, "Abundance")

# Generate plot to examine distribution of remaining missing values.
plot <- ggplotDetect(filt_protein, "Abundance") +
  ggtitle("Protein missing value distribution")

# Figure
caption <- c("Protein level missing value distribution.")
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = caption

# Impute the remaining number of missing values with KNN.
imp_protein <- impute_KNN(filt_protein, "Abundance")

# Final normalization with TMM.
TMM_protein <- normalize_TMM(imp_protein, "Abundance")

dim(TMM_protein)

# Save plot.
#file <- paste0(outputfigsdir, "/", outputMatName, 
#               "_Protein_Missing_Value_Density_plots.pdf")
#ggsavePDF(plot, file)

# Save as tiff. 
file <- paste0(outputfigsdir, "/", outputMatName, "_Missing_Protein_Values.pdf")
ggsave(file,plot)

#-------------------------------------------------------------------------------
#' ## Examine the TMM Normalized protein level data.
#-------------------------------------------------------------------------------
#+ eval = TRUE

data_in <- TMM_protein
title <- "TMM Normalized protein"

# Generate boxplot.
# Adjust color vector if samples were removed.
# Cortex outliers = 1x Ube3a, and 1x Syngap1
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

# Figure.
caption <- strwrap("TMM normalized protein. A. Boxplot. B. Density plot.
                   C. Mean SD plot. D. MDS plot", width = Inf, simplify = TRUE)
#+ warning = FALSE, echo = FALSE, fig.align = 'center', fig.cap = caption
fig <- plot_grid(p1, p2, p3, p4, labels = "auto")
fig 

# Save plots.
#file <- paste0(outputfigsdir, "/", outputMatName, "_TMM_Normalized_Protein.pdf")
#ggsavePDF(plots = list(p1, p2, p3, p4), file)

# Save tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_TMM_Normalized_Protein.pdf")
ggsave(file,fig)

#-------------------------------------------------------------------------------
#' ## IntraBatch statistical testing: EdgeR GLM.
#-------------------------------------------------------------------------------
#' Differential protein expression is evaluated using a generalized linear model
#' (GLM) implemented by `EdgeR's` __glmQLFit()__ and __glmQLFTest()__ functions.
#' Comparisons are made within a batch (WT vs KO or HET).
#' 
#+ eval = TRUE

# data is...
data_in <- TMM_protein # normalized, imputed, protein level data.
data_in[1:5, 1:5]
dim(data_in)

# Format data for EdgeR. Remove QC Samples!
cols <- grepl("Abundance", colnames(TMM_protein))
dm <- as.matrix(data_in[, cols])
rownames(dm) <- data_in$Accession
out <- grepl("QC", colnames(dm))
dm <- dm[, !out]
dim(dm)

# Check, there should be no missing values.
sum(is.na(dm)) == 0

# Check, traits and data in matching order?
traits <- sample_info
traits_sub <- subset(traits, !traits$SampleType == "QC")
traits_sub <- subset(traits_sub, traits_sub$ColumnName %in% colnames(dm))
all(colnames(dm) == traits_sub$ColumnName)

# Create DGEList object with mapping to genotype (group).
group <- as.factor(traits_sub$Sample.Model)
y_DGE <- DGEList(counts = dm, group = group)

# Basic design matrix for GLM.
design <- model.matrix(~ 0 + group, data = y_DGE$samples)
colnames(design)[1:length(unique(traits_sub$Sample.Model))] <- levels(y_DGE$samples$group)
design

# Estimate dispersion:
y_DGE <- estimateDisp(y_DGE, design, robust = TRUE)

# PlotBCV
plot <- ggplotBCV(y_DGE)
plot

# Fit a general linear model.
fit <- glmQLFit(y_DGE, design, robust = TRUE)

# Examine the QL fitted dispersion.
plotQLDisp(fit)

# Create a list of contrasts for pairwise comparisons.
contrasts <- list(
  WTvKO.Shank2 <- makeContrasts(KO.Shank2 - WT.Shank2, levels = design),
  WTvKO.Shank3 <- makeContrasts(KO.Shank3 - WT.Shank3, levels = design),
  WTvHET.Syngap1 <- makeContrasts(HET.Syngap1 - WT.Syngap1, levels = design),
  WTvKO.Ube3a <- makeContrasts(KO.Ube3a - WT.Ube3a, levels = design)
)

# Call glmQLFTest() to evaluate differences in contrasts.
qlf <- lapply(contrasts, function(x) glmQLFTest(fit, contrast = x))

## Determine number of significant results with decideTests().
#  Default FDR is 0.05.
summary_table <- lapply(qlf, function(x) summary(decideTests(x)))
overall <- t(matrix(unlist(summary_table), nrow = 3, ncol = 4))
rownames(overall) <- unlist(lapply(contrasts, function(x) colnames(x)))
colnames(overall) <- c("Down", "NotSig", "Up")
overall <- as.data.frame(overall)
overall <- add_column(overall, Contrast = rownames(overall), .before = 1)
overall <- overall[, c(1, 3, 2, 4)]
overall$TotalSig <- rowSums(overall[, c(3, 4)])
# Table of DE candidates.
table <- tableGrob(overall, rows = NULL)
grid.arrange(table)

# Save table.
file <- paste0(outputfigsdir, "/", outputMatName, "_InraBatch_GLM_Table.tiff")
ggsave(file,table)

# Call topTags to add FDR. Gather tablurized results.
results <- lapply(qlf, function(x) topTags(x, n = Inf, sort.by = "none")$table)

# Convert logCPM column to percent WT.
# Annotate with candidate column.
results <- lapply(results, function(x) annotateTopTags(x))

# Annotate with Gene names and Entrez IDS.
results <- lapply(results, function(x) annotate_Entrez(x))

# Add fitted data to results. Sort by p-value.
names(results) <- groups
for (i in 1:length(groups)) {
  cols <- grepl(groups[i], colnames(fit$fitted.values))
  results[[i]] <- merge(results[[i]],
    log2(fit$fitted.values[, cols]),
    by = "row.names"
  )
  results[[i]]$Row.names <- NULL
  results[[i]] <- results[[i]][order(results[[i]]$FDR), ]
}

# Pvalue histograms.
p1 <- ggplotPvalHist(results[[1]], "gold1", "Shank2")
p2 <- ggplotPvalHist(results[[2]], "blue", "Shank3")
p3 <- ggplotPvalHist(results[[3]], "green", "Syngap1")
p4 <- ggplotPvalHist(results[[4]], "purple", "Ube3a")
fig <- plot_grid(p1, p2, p3, p4, labels = "auto")
fig

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_IntraBatch_GLM_PvalueHist.tiff")
ggsave(file,fig)

# Add summary table to results.
overall <- list(overall)
names(overall) <- "Summary"
results <- c(overall, results)
results_intraBatch <- results

# Save workbook.
file <- paste0(outputtabsdir, "/", outputMatName, "_IntraBatch_GLM_results.xlsx")
write.excel(results, file)

# Save figures.
#file <- paste0(outputfigsdir, "/", outputMatName, "_IntraBatch_GLM_DE_Summary.pdf")
#ggsavePDF(plots = list(table, p1, p2, p3, p4), file)

# Save to RDS.
file <- paste0(Rdatadir, "/", outputMatName, "_IntraBatch_Results.RDS")
saveRDS(results_intraBatch, file)

#-------------------------------------------------------------------------------
#' ## IntraBatch GO and KEGG enrichment testing with EdgeR.
#-------------------------------------------------------------------------------
#' GO and KEGG enrichment analyis are performed using `EdgeR's`` __goana()__ 
#' and __kegga()__ functions.
#' 
#+ eval = FALSE

## Perform GO and KEGG testing.
# edgeR_GSE() is a wrapper around the goana() and kegga() functions from EdgeR.
# This function operates on the QLF or ET objects.Rownames should be UniprotIDs.
# Proteins with an FDR <0.1 are be considered significant.

# Use lapply to generate GSE results.
# This will take a few minutes.
GSE_results <- lapply(qlf, function(x) edgeR_GSE(x, FDR = 0.1, filter = TRUE))

# Name GSE_results.
GSE_results <- unlist(GSE_results, recursive = FALSE)
names(GSE_results) <- paste(rep(groups, each = 2), names(GSE_results))

# Save workbook.
file <- paste0(outputtabsdir, "/", outputMatName, "_IntraBatch_GLM_GSE_Results.xlsx")
write.excel(GSE_results, file)

#-------------------------------------------------------------------------------
#' ## Perform moderated EB regression of genetic strain as a covariate.
#-------------------------------------------------------------------------------
#' Moderated Empirical Bayes (EB) regression as implemented by the `WGCNA` 
#' function is __empiricalBayesLM()__ is performed to remove the affect of
#' genetic background.
#' 
#+ eval = FALSE

# Data is...
# Check, there should be no missing values.
sum(is.na(TMM_protein)) == 0
data_in <- na.omit(TMM_protein)
dim(data_in)
colID <- "Abundance"
traits <- sample_info

# Remove the QC Data
out <- grepl("QC", colnames(data_in))
data <- data_in[, !out]
traits <- subset(traits, !traits$SampleType == "QC")
traits <- subset(traits, traits$ColumnName %in% colnames(data_in))

# Format data for EBLM
cols <- grep("Abundance", colnames(data))
data <- log2(as.matrix(data[, cols]))
rownames(data) <- data_in$Accession
data <- t(data)
data[1:5, 1:5]
dim(data)

# Check, data and traits are in matching order?
all(rownames(data) == traits$ColumnName)

# Define covariates.
status <- traits$SampleType
sex <- as.factor(traits$Sex)
age <- as.numeric(traits$Age)
batch <- as.factor(traits$PrepDate)
strain <- as.factor(traits$Model)

# Design, we will perform regression on strain (mouse genetic background).
design <- as.data.frame(cbind(status, strain, batch, sex, age))
covariates <- cbind(design$strain)

# Eblm regression.
fit.eblm <- empiricalBayesLM(data,
  removedCovariates = covariates,
  fitToSamples = design$status == "WT"
)

# Get fitted data.
data.fit <- fit.eblm$adjustedData

# Examine overall PCA before and after EB regression.
colors <- traits$Color
plot1 <- ggplotPCA(t(data), traits, colors, title = "2D PCA Plot (Pre-Regression)")
plot2 <- ggplotPCA(t(data.fit), traits, colors, title = "2D PCA Plot (Post-Regression)")
fig <- plot_grid(plot1, plot2)
fig 

# Just post-regression:
plot2

# Save plots.
#file <- paste0(outputfigsdir, "/", outputMatName, "_InterBatch_eBLM_Regression_PCA.pdf")
#ggsavePDF(plots = list(plot1, plot2), file)

# Save tiffs. 
files <- paste0(outputfigsdir, "/", outputMatName, c("pre","post"),"_InterBatch_eBLM_Regression_PCA.tiff")
quiet(mapply(ggsave,files,list(plot1,plot2)))

#-------------------------------------------------------------------------------
#' ## Reformat final normalized, regressed data for TAMPOR processing.
#-------------------------------------------------------------------------------
#' Data are reformatted for TAMPOR normalization in the `2_TMT_Analysis.R` script.
#' 
#+ eval = FALSE 

# Data is...
data_in <- 2^t(data.fit)
data_in[1:5, 1:5]  # un-log
dim(data_in)

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
head(row_names)

# Gather data, set row names
cols <- grep("Abundance", colnames(data_in))
data_out <- as.matrix(data_in[, cols])
rownames(data_out) <- row_names

# Column names are batch.channel
idx <- match(colnames(data_out), sample_info$ColumnName)
idx
col_names <- sample_info$SampleID[idx]
colnames(data_out) <- col_names
head(col_names)

# Reorder based on batch.
data_out <- data_out[, order(colnames(data_out))]

# Save as cleanDat
cleanDat <- data_out
cleanDat[1:5, 1:5]
dim(cleanDat)

# Save to Rdata
file <- paste(Rdatadir, tissue, "CleanDat_IRS_eBLM_TAMPOR_format.Rds", sep = "/")
saveRDS(cleanDat, file)
# cleanDat <- readRDS(file)

# The cortex and striatum data from this chunk will be combined for TAMPOR
# normalization and statistical testing in the
# TMT_Analysis_EBD_v14_Combined_Cortex_Striatum.R script.

#-------------------------------------------------------------------------------
#' ## InterBatch Statistical comparisons with EdgeR GLM.
#-------------------------------------------------------------------------------
#' After regression of genetic strain, differential protein expression is 
#' evaluated across batches, using all WT samples. 
#' 
#+ eval = FALSE

# Create DGEList object with EB adjusted data.
data <- 2^t(data.fit)
data[1:5, 1:5]
y_DGE <- DGEList(counts = data)

# Example, checking the the normalization with plotMD.
plotMD(cpm(y_DGE, log = TRUE), column = 2)
abline(h = 0, col = "red", lty = 2, lwd = 2)

# Create sample mapping.
traits <- sample_info
rownames(traits) <- traits$ColumnName
traits <- subset(traits, rownames(traits) %in% colnames(data))
all(traits$ColumnName == colnames(data))
group <- traits$Sample.Model
strain <- traits$Model
sex <- traits$Sex
unique(group)
group[grepl("WT", group)] <- "WT"
unique(group)
y_DGE$samples$group <- as.factor(group)

# Basic design matrix for GLM.
design <- model.matrix(~ 0 + group, data = y_DGE$samples)
colnames(design) <- levels(y_DGE$samples$group)
design

# Estimate dispersion:
y_DGE <- estimateDisp(y_DGE, design, robust = TRUE)

# PlotBCV
plot <- ggplotBCV(y_DGE)
plot

# Fit a general linear model.
fit <- glmQLFit(y_DGE, design, robust = TRUE)

# Examine the QL fitted dispersion.
plotQLDisp(fit)

# Which genes are DE among all contrasts?
# Create contrast matrix for ANOVA-like comparison.
aov_contrasts <- makeContrasts(
  WTvKO.Shank2 = KO.Shank2 - WT,
  WTvKO.Shank3 = KO.Shank3 - WT,
  WTvHET.Syngap1 = HET.Syngap1 - WT,
  WTvKO.Ube3a = KO.Ube3a - WT, levels = design
)
aov_contrasts

# Calculate ANOVA-like results.
# The QL F-test is applied to identify genes that are DE among all four groups. This
# combines the four pairwise comparisons into a single F-statistic and p-value. The top set of
# significant genes can be displayed with topTags.
aov <- glmQLFTest(fit, contrast = aov_contrasts)
aov_tt <- topTags(aov, n = Inf, sort.by = "none")

# Extract the results, and annotated with Entrez IDS and gene symbols.
res <- annotate_Entrez(aov_tt$table)
# Discard the CPM column.
res$logCPM <- NULL
results_GLMoverall <- res

# Create a list of contrasts for pairwise comparisons.
contrasts <- list(
  WTvKO.Shank2 <- makeContrasts(KO.Shank2 - WT, levels = design),
  WTvKO.Shank3 <- makeContrasts(KO.Shank3 - WT, levels = design),
  WTvHET.Syngap1 <- makeContrasts(HET.Syngap1 - WT, levels = design),
  WTvKO.Ube3a <- makeContrasts(KO.Ube3a - WT, levels = design)
)

# Call glmQLFTest() to evaluate differences in contrasts.
qlf <- lapply(contrasts, function(x) glmQLFTest(fit, contrast = x))

## Determine number of significant results with decideTests().
summary_table <- lapply(qlf, function(x) summary(decideTests(x)))
overall <- t(matrix(unlist(summary_table), nrow = 3, ncol = 4))
rownames(overall) <- unlist(lapply(contrasts, function(x) colnames(x)))
colnames(overall) <- c("Down", "NotSig", "Up")
overall <- as.data.frame(overall)
overall <- add_column(overall, Contrast = rownames(overall), .before = 1)
overall <- overall[, c(1, 3, 2, 4)]
overall$TotalSig <- rowSums(overall[, c(3, 4)])
# Table of DE candidates.
table <- tableGrob(overall, rows = NULL)
grid.arrange(table)

# Save table.
file <- paste0(outputfigsdir, "/", outputMatName, "_InterBatch_eBLM_Table.tiff")
ggsave(file,table)

# Call topTags to add FDR. Gather tablurized results.
results <- lapply(qlf, function(x) topTags(x, n = Inf, sort.by = "none")$table)

# Function to annotate DE candidates:
annotateTopTags <- function(y_TT) {
  y_TT$logCPM <- 100 * (2^y_TT$logFC)
  colnames(y_TT)[2] <- "%WT"
  colnames(y_TT)[3] <- "F Value"
  y_TT$candidate <- "no"
  y_TT[which(y_TT$FDR <= 0.10 & y_TT$FDR > 0.05), dim(y_TT)[2]] <- "low"
  y_TT[which(y_TT$FDR <= 0.05 & y_TT$FDR > 0.01), dim(y_TT)[2]] <- "med"
  y_TT[which(y_TT$FDR <= 0.01), dim(y_TT)[2]] <- "high"
  y_TT$candidate <- factor(y_TT$candidate, levels = c("high", "med", "low", "no"))
  return(y_TT)
}

# Convert logCPM column to percent WT.
# Annotate with candidate column.
results <- lapply(results, function(x) annotateTopTags(x))

# Annotate with Gene names and Entrez IDS.
results <- lapply(results, function(x) annotate_Entrez(x))

## Add GLM overall stats and fitted data to results.
GLM_stats <- results_GLMoverall[, c(8, 9, 10)]
colnames(GLM_stats)[1] <- "F Value"
colnames(GLM_stats) <- paste("Overall", colnames(GLM_stats))
names(results) <- groups
data_glm <- log2(qlf[[1]]$fitted.values)
data_glm <- merge(GLM_stats, data_glm, by = "row.names", all = TRUE)
rownames(data_glm) <- data_glm$Row.names
data_glm$Row.names <- NULL
data_glm[1:5, 1:5]

# Loop to add fitted data + GLM overall stats to pairwise comparisons in results.
for (i in 1:length(groups)) {
  cols <- c(1, 2, 3, grep(paste(groups[i], "WT", sep = "|"), colnames(data_glm)))
  results[[i]] <- merge(results[[i]], data_glm[, cols], by = "row.names")
  results[[i]]$Row.names <- NULL
  results[[i]] <- results[[i]][order(results[[i]]$PValue), ]
}

# Loop to Re-organize columns.
for (i in 1:length(groups)) {
  col_names <- colnames(results[[i]])
  colsWT <- grep("WT", col_names)[-1]
  colsKO <- grep("KO|HET", col_names)
  colsOverall <- grep("Overall", col_names)
  colsElse <- (1:length(col_names))[-c(1, 2, 3, colsWT, colsKO, colsOverall)]
  idx <- c(1, 2, 3, colsOverall, colsElse, colsKO, colsWT)
  col_names[idx]
  results[[i]] <- results[[i]][col_names[idx]]
}

# Results
results_interBatch <- results

# Pvalue Histograms (need to ingore "P overall" column):
p1 <- ggplotPvalHist(results[[1]][-c(1:6)], "gold1", "Shank2")
p2 <- ggplotPvalHist(results[[2]][-c(1:6)], "blue", "Shank3")
p3 <- ggplotPvalHist(results[[3]][-c(1:6)], "green", "Syngap1")
p4 <- ggplotPvalHist(results[[4]][-c(1:6)], "purple", "Ube3a")
fig <- plot_grid(p1, p2, p3, p4, labels = "auto")
fig

# Save plots.
#file <- paste0(outputfigsdir, "/", outputMatName, "_InterBatch_GLM_PvalHist.pdf")
#ggsavePDF(plots = list(table, p1, p2, p3, p4), file)

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_InterBatch_GLM_PvalHist.tiff")
ggsave(file,fig)

## Write results to file.
# Add summary table to results.
overall <- list(overall)
names(overall) <- "Summary"
results <- c(overall, results)

# Save workbook.
file <- paste0(outputtabsdir, "/", outputMatName, "_InterBatch_eBLM_GLM_Results.xlsx")
write.excel(results, file)

# Save to RDS.
file <- paste0(Rdatadir, "/", outputMatName, "_InterBatch_Results.RDS")
saveRDS(results_interBatch, file)

#-------------------------------------------------------------------------------
#' ## InterBatch GO and KEGG enrichment testing.
#-------------------------------------------------------------------------------
#+ eval = FALSE

## Perform GO and KEGG testing.
# qlf_GSE() is a wrapper around the goana() and kegga() functions from EdgeR.
# This function operates on the QLF object.
# Proteins with FDR <0.05 will be considered differentially expressed.

# Use lapply to generate GSE results with custom function edgeR_GSE().
# This will take a few minutes.
GSE_results <- lapply(qlf, function(x) edgeR_GSE(x, FDR = 0.05, filter = TRUE))

# Name GSE_results.
GSE_results <- unlist(GSE_results, recursive = FALSE)
names(GSE_results) <- paste(rep(groups, each = 2), names(GSE_results))

# Initiate an excel workbook.
wb <- createWorkbook()

# Loop to add a worksheets:
for (i in 1:length(GSE_results)) {
  df <- GSE_results[[i]]
  addWorksheet(wb, sheetName = names(GSE_results[i]))
  writeData(wb, sheet = i, df)
}

# Save workbook.
file <- paste0(outputtabsdir, "/", outputMatName, "_InterBatch_eBLM_GLM_GSE_Results.xlsx")
saveWorkbook(wb, file, overwrite = TRUE)

#-------------------------------------------------------------------------------
#' ## Generate protein boxplots for significantly DE proteins (Interbatch).
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Subset the data. Keep proteins with p.adj overall <0.05.
data_sub <- results$Shank2[results_interBatch$Shank2$`Overall PValue` < 0.05, ]
keep <- data_sub$Uniprot
exprDat <- log2(y_DGE$counts)
idx <- rownames(exprDat) %in% keep
exprDat <- exprDat[idx, ]

# Annotate exprDat rows as gene|uniprot
exprDat <- log2(y_DGE$counts)
Uniprot <- rownames(exprDat)
Gene <- mapIds(
  org.Mm.eg.db,
  keys = Uniprot,
  column = "SYMBOL",
  keytype = "UNIPROT",
  multiVals = "first"
)
rownames(exprDat) <- paste(Gene, rownames(exprDat), sep = "|")

# Insure all WT samples are annotated as WT in traits.
traits_temp <- traits
traits_temp$Sample.Model[grepl("WT", traits_temp$Sample.Model)] <- "WT"
traits_temp <- subset(traits_temp, rownames(traits_temp) %in% colnames(exprDat))

# Generate plots.
plot_list <- ggplotProteinBoxes(
  data_in = exprDat,
  interesting.proteins = rownames(exprDat),
  dataType = "Relative Abundance",
  traits = traits_temp,
  order = c(1, 2, 3, 4, 5),
  scatter = TRUE
)
# Add custom colors.
colors <- c("gray", "yellow", "blue", "green", "purple")
plot_list <- lapply(plot_list, function(x) x + scale_fill_manual(values = colors))

# Example plot.
plot_list[[1]]

## Add significance stars.
# Build a df with statistical results.
stats <- lapply(results_interBatch, function(x)
  as.data.frame(cbind(Uniprot = x$Uniprot, FDR = x$FDR)))
names(stats) <- c("KO.Shank2", "KO.Shank3", "HET.Syngap1", "KO.Ube3a")
df <- stats %>% reduce(left_join, by = "Uniprot")
colnames(df)[c(2:ncol(df))] <- names(stats)

# Annotate rows as gene|uniprot
Uniprot <- df$Uniprot
Gene <- mapIds(
  org.Mm.eg.db,
  keys = Uniprot,
  column = "SYMBOL",
  keytype = "UNIPROT",
  multiVals = "first"
)
rownames(df) <- paste(Gene, Uniprot, sep = "|")
df$Uniprot <- NULL
stats <- df

# Example plot.
plot <- plot_list[[1]]
annotate_stars(plot, stats)

# Loop to add stars.
plot_list <- lapply(plot_list, function(x) annotate_stars(x, stats))

# Top proteins.
p1 <- plot_list$`Shank2|Q80Z38`
p2 <- plot_list$`Shank3|Q4ACU6`
p3 <- plot_list$`Syngap1|F6SEU4`
p4 <- plot_list$`Ube3a|O08759`
fig <- plot_grid(p1,p2,p3,p4)
fig

# Save to pdf.
#file <- paste0(outputfigsdir, "/", tissue, "_WGCNA_Analysis_InterBatch_ProteinBoxPlots.pdf")
#ggsavePDF(plots = plot_list, file)

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "TopProteins_BoxPlots.tiff")
ggsave(file,fig)

#-------------------------------------------------------------------------------
#' ## Render RMarkdown report.
#-------------------------------------------------------------------------------
#' This script is formatted for automated rendering of an RMarkdown report.
#'
#+ eval = FALSE

# Code directory. 
dir <- paste(rootdir,"Code",sep="/")
file <- paste(dir,"1_TMT_Analysis.R",sep="/")

# Save and render.
rstudioapi::documentSave()
rmarkdown::render(file)
