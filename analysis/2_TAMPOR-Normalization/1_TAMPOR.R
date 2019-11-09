#!/usr/bin/env Rscript

#' ---
#' title: 1_TAMPOR.R
#' description: TAMPOR Normalization of preprocessed TMT data.
#' authors: Tyler W Bradshaw, Eric B Dammer.
#' ---

#-------------------------------------------------------------------------------
## Prepare the workspace.
#-------------------------------------------------------------------------------
# Prepare the R workspace for the analysis. Load custom functions and prepare
# the project directory for saving output files.

rm(list = ls())
if (!is.null(dev.list())) {
  dev.off()
}
cat("\f")
options(stringsAsFactors = FALSE)

# Load required packages.
suppressPackageStartupMessages({
  library(dplyr)
  library(reshape2)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(edgeR)
  library(gridExtra)
  library(gtable)
  library(grid)
  library(anRichment)
  library(openxlsx)
  library(org.Mm.eg.db)
})

# Define tisue type:
tissue <- "Combined"

# Set the working directory.
here <- getwd()
rootdir <- dirname(dirname(here))

# Set any other directories.
functiondir <- paste(rootdir, "R", sep = "/")
datadir <- paste(rootdir, "data", sep = "/")
Rdatadir <- paste(rootdir, "rdata", sep = "/")
outputfigs <- paste(rootdir, "figs", tissue, sep = "/")
outputtabs <- paste(rootdir, "tables", sep = "/")

# Load required custom functions.
source_myfun <- function() {
  myfun <- list.files(functiondir, pattern = ".R", full.names = TRUE)
  invisible(sapply(myfun, source))
}
source_myfun()

# Define prefix for output figures and tables.
outputMatName <- paste0("2_", tissue)

# Globally set ggplots theme.
ggtheme()

# Store any plots in list.
all_plots <- list()

#-------------------------------------------------------------------------------
## Merge cortex and striatum data.
#-------------------------------------------------------------------------------
# We will utilize TAMPOR to combine the Cortex and Striatum datasets. Merge the
# preprocessed data and traits files.

# Merge traits data.
# Load the cortex and striatum traits files.
inputTraitsCSV <- c(
  "4227_TMT_Cortex_Combined_traits.csv",
  "4227_TMT_Striatum_Combined_traits.csv"
)

# Load the sample info into a list, traits.
files <- paste(datadir, inputTraitsCSV, sep = "/")
traits <- lapply(files, function(x) read.csv(file = x, header = TRUE))

# Bind the elements of the list.
traits <- do.call(rbind, traits)

# Rename SampleIDs
idx <- c(45:nrow(traits))
vars <- strsplit(traits$SampleID, "\\.")[idx]
new_batch <- c(rep("b5", 11), rep("b6", 11), rep("b7", 11), rep("b8", 11))
channel <- sapply(vars, "[", 2)
ids <- paste(new_batch, channel, sep = ".")
traits$SampleID[idx] <- ids

# Insure that rownames are new sampleIDS
rownames(traits) <- traits$SampleID

# Add column for tissue type.
traits$Tissue <- c(rep("Cortex", 44), rep("Striatum", 44))

# Add column for batch.
traits$Batch <- sapply(strsplit(rownames(traits), "\\."), "[", 1)
alltraits <- traits

## Merge expression data.
# Load the Cortex and Striatum cleanDat.
myfiles <- file.path(Rdatadir, c(
  "1_Cortex_cleanDat.RData",
  "1_Striatum_cleanDat.RData"
))
data <- list(
  "Cortex"  = readRDS(myfiles[1]),
  "Striatum" = readRDS(myfiles[2])
)

# Fortify and add accession column
data_fort <- lapply(
  data,
  function(x) add_column(as.data.frame(x), ID = rownames(x), .before = 1)
)

# Bind by Accession
data_merge <- data_fort %>% purrr::reduce(left_join, by = "ID")
rownames(data_merge) <- data_merge$ID
data_merge$ID <- NULL

# Remove rows with missing values.
allDat <- na.omit(data_merge)

## Clean-up formatting for TAMPOR.
# Batch b1-b4 are cortex. Batches b5-b8 are straitum.
col_names <- colnames(allDat)
# Cortex = group1...
group1 <- col_names[grepl(".x", col_names)]
group1 <- gsub(".x", "", group1)
# Striatum = group2...
group2 <- col_names[!grepl(".x", col_names)]
group2 <- gsub(".y", "", group2)
group2 <- gsub("b1", "b5", group2)
group2 <- gsub("b2", "b6", group2)
group2 <- gsub("b3", "b7", group2)
group2 <- gsub("b4", "b8", group2)

# Change column names to batch.channel.
colnames(allDat) <- c(group1, group2)

# GIS index is all WT samples.
# WT Cortex and WT Striatum scaled by TAMPOR to be the same.
controls <- alltraits$SampleID[grepl("WT", alltraits$SampleType)]

# Save merged traits file.
myfile <- paste0(Rdatadir, "/", outputMatName, "_traits.RData")
saveRDS(alltraits, file = myfile)

#-------------------------------------------------------------------------------
## Perform TAMPOR normalization.
#-------------------------------------------------------------------------------

# Insure than any samples that were removed from cleanDat are removed from
# traits (any outliers identified in previous scripts).
traits <- alltraits[rownames(alltraits) %in% colnames(allDat), ]

# Perform TAMPOR.
results <- TAMPOR(
  dat = allDat,
  traits = traits,
  batchPrefixInSampleNames = TRUE,
  samplesToIgnore = "None",
  GISchannels = controls,
  parallelThreads = 8
)

# Collect normalize relative abundance data.
cleanDat <- results$cleanRelAbun

#-------------------------------------------------------------------------------
## Remove any sample outliers.
#-------------------------------------------------------------------------------

# Remove QC samples.
out <- colnames(cleanDat) %in% rownames(traits)[traits$SampleType == "QC"]
data_in <- log2(cleanDat[, !out])

# Illustrate Oldham's sample connectivity.
sample_connectivity <- ggplotSampleConnectivity(data_in, log = TRUE, colID = "b.")
plot <- sample_connectivity$connectivityplot +
  ggtitle("Sample Connectivity post-TAMPOR")

# Store plot.
all_plots[["TAMPOR_Oldham_Outliers"]] <- plot

# Loop to identify Sample outliers using Oldham's connectivity method.
n_iter <- 5
threshold <- -3.0
out_samples <- list()

# Loop:
for (i in 1:n_iter) {
  data_temp <- data_in
  oldham <- ggplotSampleConnectivity(data_temp, log = TRUE, colID = "b", threshold = -3)
  bad_samples <- rownames(oldham$table)[oldham$table$Z.Ki < threshold]
  print(paste(length(bad_samples), " outlier sample(s) identified in iteration ", i, ".", sep = ""))
  if (length(bad_samples) == 0) bad_samples <- "none"
  out_samples[[i]] <- bad_samples
  out <- grepl(paste(unlist(out_samples), collapse = "|"), colnames(data_in))
  data_in <- data_in[, !out]
}

# Outlier samples.
bad_samples <- unlist(out_samples)
message(paste("Outlier samples:", traits$Sample.Model[rownames(traits) %in% bad_samples]))

# Save data with QC samples, but outliers removed to file.
cleanDat <- cleanDat[, !colnames(cleanDat) %in% bad_samples]
myfile <- file.path(Rdatadir, "2_Combined_cleanDat.RData")
saveRDS(cleanDat, myfile)

#-------------------------------------------------------------------------------
## Examine sample clustering with MDS and PCA post-TAMPOR Normalization.
#-------------------------------------------------------------------------------

# Insure that any outlier samples have been removed.
traits <- alltraits[rownames(alltraits) %in% colnames(cleanDat), ]

# Remove QC data.
out <- colnames(cleanDat) %in% rownames(traits)[traits$SampleType == "QC"]
data_in <- log2(cleanDat[, !out])

# Check, traits and cleanDat should match data.
traits <- traits[match(colnames(data_in), rownames(traits)), ]
if (!all(rownames(traits) == colnames(data_in))) {
  stop("data do not match traits.")
}

## PCA Plots.
colors <- traits$Color
traits$ColumnName <- rownames(traits)

# Cortex and striatum.
idx <- traits$Tissue == "Cortex"
idy <- traits$Tissue == "Striatum"
plot1 <- ggplotPCA(log2(data_in[, idx]), traits,
  colID = "b.",
  colors[idx], title = "Cortex"
)
plot2 <- ggplotPCA(log2(data_in[, idy]), traits,
  colID = "b.",
  colors[idy], title = "Striatum"
)
plot3 <- ggplotPCA(log2(data_in), traits, colors,
  colID = "b.",
  title = "Combined"
)

# Store in list.
all_plots[["Cortex_postTAMPOR_PCA"]] <- plot1
all_plots[["Striatum_postTAMPOR_PCA"]] <- plot2
all_plots[["Combined_postTAMPOR_PCA"]] <- plot3

#-------------------------------------------------------------------------------
## Create protein identifier map.
#-------------------------------------------------------------------------------

ids <- rownames(cleanDat)
symbol <- sapply(strsplit(ids, "\\|"), "[", 1)
uniprot <- sapply(strsplit(ids, "\\|"), "[", 2)

# Map uniprot to entrez.
suppressMessages({
  entrez <- mapIds(org.Mm.eg.db,
    keys = uniprot,
    column = "ENTREZID", keytype = "UNIPROT", multiVals = "first"
  )
})
# If missing map gene to entrez.
is_missing <- is.na(entrez)
suppressMessages({
  entrez[is_missing] <- mapIds(org.Mm.eg.db,
    keys = symbol[is_missing],
    column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"
  )
})
is_missing <- is.na(entrez)

# Map remainder by hand.
# ids[is_missing]
entrez[is_missing] <- c(102631912, 14070, 18517)

# Check
if (sum(is.na(entrez)) > 0) {
  stop("Not all genes mapped to entrez!")
}

# Use entrez IDs to get consistent gene names.
suppressMessages({
  gene <- mapIds(org.Mm.eg.db,
    keys = entrez,
    column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"
  )
})

# Check
if (sum(is.na(gene)) > 0) {
  stop("Not all genes mapped to symbols!")
}

# Protein identifier map.
prot_map <- data.frame(ids, uniprot, entrez, gene)

# Save to Rdata.
myfile <- file.path(Rdatadir, "2_Prot_Map.RData")
saveRDS(prot_map, myfile)

#-------------------------------------------------------------------------------
## EdgeR statistical comparisons post-TAMPOR.
#-------------------------------------------------------------------------------

# Remove QC samples prior to passing data to EdgeR.
out <- alltraits$SampleType[match(colnames(cleanDat), rownames(alltraits))] == "QC"
data_in <- cleanDat[, !out]

# Summarize proteins...
nprots <- dim(data_in)[1]
nsamples <- dim(data_in)[2]
message(paste(nprots, "proteins identified in", nsamples, "samples."))

# Create DGEList object...
y_DGE <- DGEList(counts = data_in)

# TMM Normalization.
y_DGE <- calcNormFactors(y_DGE)

# Create sample mapping to Tissue.Genotype.
# Group WT Cortex samples and WT Striatum samples together.
traits <- subset(alltraits, rownames(traits) %in% colnames(data_in))
traits <- traits[match(colnames(data_in), rownames(traits)), ]
all(traits$SampleID == colnames(data_in))
group <- paste(traits$Tissue, traits$Sample.Model, sep = ".")
group[grepl("Cortex.WT", group)] <- "Cortex.WT"
group[grepl("Striatum.WT", group)] <- "Striatum.WT"
traits$group <- group
y_DGE$samples$group <- as.factor(group)

# Basic design matrix for GLM.
design <- model.matrix(~ 0 + group, data = y_DGE$samples)
colnames(design) <- levels(y_DGE$samples$group)

# Estimate dispersion:
y_DGE <- estimateDisp(y_DGE, design, robust = TRUE)

# PlotBCV
# plot <- ggplotBCV(y_DGE)
# all_plots[["TAMPOR_BCV"]] <- plot

# Fit a general linear model.
fit <- glmQLFit(y_DGE, design, robust = TRUE)

# Generate contrasts.
g1 <- colnames(design)[grepl("Cortex", colnames(design))][-5]
g2 <- colnames(design)[grepl("Striatum", colnames(design))][-5]
cont1 <- makePairwiseContrasts(list(g1), list("Cortex.WT"))
cont2 <- makePairwiseContrasts(list(g2), list("Striatum.WT"))

# Make contrasts for EdgeR.
# For some reason loops or lapply dont work with the makeContrasts function.
contrasts <- list(
  makeContrasts(cont1[1], levels = design),
  makeContrasts(cont1[2], levels = design),
  makeContrasts(cont1[3], levels = design),
  makeContrasts(cont1[4], levels = design),
  makeContrasts(cont2[1], levels = design),
  makeContrasts(cont2[2], levels = design),
  makeContrasts(cont2[3], levels = design),
  makeContrasts(cont2[4], levels = design)
)
names(contrasts) <- unlist({
  lapply(contrasts, function(x) sapply(strsplit(colnames(x), " "), "[", 1))
})

# Call glmQLFTest() to evaluate differences in contrasts.
qlf <- lapply(contrasts, function(x) glmQLFTest(fit, contrast = x))
names(qlf) <- names(contrasts)

## Determine number of significant results with decideTests().
summary_table <- lapply(qlf, function(x) summary(decideTests(x)))
overall <- t(matrix(unlist(summary_table), nrow = 3, ncol = 8))
rownames(overall) <- unlist(lapply(contrasts, function(x) colnames(x)))
colnames(overall) <- c("Down", "NS", "Up")
overall <- as.data.frame(overall)
row_names <- sapply(strsplit(rownames(overall), " - "), "[", 1)
row_names <- gsub(".KO.|.HET.", " ", row_names)
overall <- add_column(overall, Experiment = row_names, .before = 1)
overall <- overall[, c(1, 3, 2, 4)]
overall$"Total Sig" <- rowSums(overall[, c(3, 4)])
overall <- overall[c(2, 6, 3, 7, 1, 5, 4, 8), ] # Reorder.

# Table of DA candidates.
# Modify tables theme to change font size.
# Cex is a scaling factor relative to the defaults.
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params = list(cex = 0.75)),
  colhead = list(fg_params = list(cex = 0.75)),
  rowhead = list(fg_params = list(cex = 0.75))
)

# Create table and add borders.
mytable <- tableGrob(overall, rows = NULL, theme = mytheme)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 2, b = nrow(mytable), l = 1, r = ncol(mytable)
)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 1, l = 1, r = ncol(mytable)
)

# Store in list.
all_plots[["TAMPOR_DE_Table"]] <- mytable

# Call topTags to add FDR. Gather tabularized results.
glm_results <- lapply(qlf, function(x) topTags(x, n = Inf, sort.by = "none")$table)

# Convert logCPM column to percent WT and annotate with candidate column.
glm_results <- lapply(glm_results, function(x) annotateTopTags(x))

# Use prot_map to annotate glm_results with entrez Ids and gene symbols.
for (i in 1:length(glm_results)) {
  x <- glm_results[[i]]
  idx <- match(rownames(x), prot_map$ids)
  x <- add_column(x, "Gene|Uniprot" = prot_map$ids[idx], .before = 1)
  x <- add_column(x, "Uniprot" = prot_map$uniprot[idx], .after = 1)
  x <- add_column(x, "Entrez" = prot_map$entrez[idx], .after = 2)
  x <- add_column(x, "Symbol" = prot_map$gene[idx], .after = 3)
  glm_results[[i]] <- x
}

# Add expression data.
for (i in 1:length(glm_results)) {
  namen <- names(glm_results)[i]
  df <- glm_results[[i]]
  comparison <- contrasts[[namen]]
  groups <- rownames(comparison)[!comparison == 0]
  samples <- traits$SampleID[traits$group %in% groups]
  dat <- cleanDat[, samples]
  colnames(dat) <- traits$ColumnName[match(colnames(dat), traits$SampleID)]
  dat <- dat[, c(grep("HET|KO", colnames(dat)), grep("WT", colnames(dat)))]
  out <- merge(df, log2(dat), by = "row.names")
  glm_results[[i]] <- out
}

# Sort by pvalue.
glm_results <- lapply(glm_results, function(x) x[order(x$PValue), ])

# Reorder by genotype.
idx <- c(
  grep("Shank2", names(glm_results)),
  grep("Shank3", names(glm_results)),
  grep("Syngap1", names(glm_results)),
  grep("Ube3a", names(glm_results))
)
glm_results <- glm_results[idx]

# Final renaming.
namen <- unlist({
  lapply(
    lapply(strsplit(gsub("HET.|KO.", "", names(glm_results)), "\\."), rev),
    function(x) paste(x, collapse = " ")
  )
})
names(glm_results) <- namen

# Remove Row.names column.
f <- function(x) {
  x$Row.names <- NULL
  return(x)
}
glm_results <- lapply(glm_results, f)

# Save results to file.
myfile <- file.path(outputtabs, paste0(outputMatName, "_GLM_Results.xlsx"))
write_excel(glm_results, myfile)

#-------------------------------------------------------------------------------
## Volcano plots for each genotype.
#-------------------------------------------------------------------------------

# Add column for genotype and unique ID to results in list.
output <- list()
for (i in 1:length(glm_results)) {
  Experiment <- names(glm_results)[i]
  df <- glm_results[[i]]
  df <- add_column(df, Experiment, .before = 1)
  output[[i]] <- df[, c(1:11)]
}

# Merge the results.
df <- do.call(rbind, output)

# Add column for genotype specific colors.
colors <- as.list(c("#FFF200", "#00A2E8", "#22B14C", "#A349A4"))
names(colors) <- c("Shank2", "Shank3", "Syngap1", "Ube3a")
df$Color <- unlist(colors[sapply(strsplit(df$Experiment, "\\ "), "[", 1)])

# Split by genotype (color).
results <- split(df, df$Color)
names(results) <- names(colors)[match(names(results), colors)]

# Function for producing volcano plots.
ggplotVolcanoPlot <- function(df) {
  df$x <- df[, grep("FC", colnames(df))]
  df$y <- -log10(df[, grep("PValue", colnames(df))])
  logic <- df$FDR < 0.05
  df$Color[!logic] <- "gray"
  df$Color <- as.factor(df$Color)
  y_int <- -1 * log10(max(df$PValue[df$FDR < 0.05]))
  plot <- ggplot(data = df, aes(x = x, y = y, color = Color)) +
    geom_point(size = 3, alpha = 0.5) + scale_color_manual(values = levels(df$Color)) +
    geom_hline(yintercept = y_int, linetype = "dashed", color = "black", size = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.6) +
    # geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 0.6) +
    xlab(expression(bold(Log[2](Fold ~ Change)))) +
    ylab(expression(bold(-Log[10](Pvalue)))) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 12, face = "bold"),
      axis.title.y = element_text(color = "black", face = "bold", size = 11, angle = 90, vjust = 0.5),
      axis.title.x = element_text(color = "black", face = "bold", size = 11, angle = 0, hjust = 0.5, vjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = "none"
    )
  # Add annotation.
  ypos <- unlist(ggplot_build(plot)$layout$panel_params[[1]][8])
  xpos <- unlist(ggplot_build(plot)$layout$panel_params[[1]][1])
  plot <- plot + annotate("text",
    x = xpos[1] + 0.3 * (xpos[2] - xpos[1]),
    y = y_int + 0.04 * (ypos[2] - ypos[1]),
    label = "FDR < 0.05", size = 4
  )
  return(plot)
}

# Generate plots.
plots <- lapply(as.list(names(colors)), function(x) ggplotVolcanoPlot(results[[x]]))
names(plots) <- names(colors)

# Add titles.
for (i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + ggtitle(names(plots)[i])
}

# Store plots in list...
vp1 <- plots$Shank2
vp2 <- plots$Shank3
vp3 <- plots$Syngap1
vp4 <- plots$Ube3a

all_plots[["Shank2_VP"]] <- vp1
all_plots[["Shank3_VP"]] <- vp2
all_plots[["Syngap1_VP"]] <- vp3
all_plots[["Ube3a_VP"]] <- vp4

#-------------------------------------------------------------------------------
## GO enrichment analaysis for DEPs from each genotype.
#-------------------------------------------------------------------------------

# Build a df with the combined statistical results.
stats <- lapply(glm_results, function(x) data.frame(Uniprot = x$Uniprot, FDR = x$FDR))
names(stats) <- names(glm_results)
df <- stats %>% purrr::reduce(left_join, by = "Uniprot")
colnames(df)[c(2:ncol(df))] <- names(stats)

## Prepare a matrix of class labels (colors) to pass to enrichmentAnalysis().
labels <- data.frame(
  Shank2 = df$"Shank2 Cortex" < 0.05 | df$"Shank2 Striatum" < 0.05,
  Shank3 = df$"Shank3 Cortex" < 0.05 | df$"Shank3 Striatum" < 0.05,
  Syngap1 = df$"Syngap1 Cortex" < 0.05 | df$"Syngap1 Striatum" < 0.05,
  Ube3a = df$"Ube3a Cortex" < 0.05 | df$"Ube3a Striatum" < 0.05
)
rownames(labels) <- df$Uniprot

# Convert TRUE to column names.
logic <- labels == TRUE # 1 will become TRUE, and 0 will become FALSE.
# Loop through each column to replace 1 with column header (color).
for (i in 1:ncol(logic)) {
  col_header <- colnames(labels)[i]
  labels[logic[, i], i] <- col_header
}

# Map Uniprot IDs to Entrez.
entrez <- prot_map$entrez[match(rownames(labels), prot_map$uniprot)]

# Insure that labels is a matrix.
labels <- as.matrix(labels)

# The labels matrix and vector of cooresponding entrez IDs
# will be passed to enrichmentAnalysis().

# Build a GO annotation collection:
myfile <- file.path(Rdatadir, "musGOcollection.RData")
if (file.exists(myfile)) {
  musGOcollection <- readRDS(myfile)
} else {
  musGOcollection <- buildGOcollection(organism = "mouse")
  saveRDS(musGOcollection, file.path(Rdatadir, "musGOcollection.RData"))
}

# Perform GO analysis for each module using hypergeometric (Fisher.test) test.
# As implmented by the WGCNA function enrichmentAnalysis().
# FDR is the BH adjusted p-value.
# Insure that the correct background (used as reference for enrichment)
# has been selected!
# useBackgroud = "given" will use all given genes as reference background.

GOenrichment <- enrichmentAnalysis(
  classLabels = labels,
  identifiers = entrez,
  refCollection = musGOcollection,
  useBackground = "given",
  threshold = 0.05,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "FALSE"
)

# Collect the results.
results_GOenrichment <- list()
for (i in 1:length(GOenrichment$setResults)) {
  results_GOenrichment[[i]] <- GOenrichment$setResults[[i]]$enrichmentTable
}
names(results_GOenrichment) <- colnames(labels)

# Save as excel workbook.
myfile <- file.path(outputtabs, paste0(outputMatName, "_GO_Analysis.xlsx"))
write_excel(results_GOenrichment, myfile)

#-------------------------------------------------------------------------------
## Condition overlap plot.
#-------------------------------------------------------------------------------

# Load statistical results..
results <- glm_results

# Combine by FDR.
stats <- lapply(results, function(x) {
  data.frame(Uniprot = x$Uniprot, Gene = x$Gene, FDR = x$FDR)
})
names(stats) <- names(results)
df <- stats %>% purrr::reduce(left_join, by = c("Uniprot", "Gene"))
colnames(df)[c(3:ncol(df))] <- paste("FDR", names(stats), sep = ".")
rownames(df) <- paste(df$Gene, df$Uniprot, sep = "|")
df$Uniprot <- NULL
df$Gene <- NULL

# Gather sigProts.
sigProts <- list()
for (i in 1:ncol(df)) {
  idx <- df[, i] < 0.05
  sigProts[[i]] <- rownames(df)[idx]
}
names(sigProts) <- names(stats)

# Build a matrix showing overlap.
col_names <- names(stats)
row_names <- names(stats)

# All possible combinations.
contrasts <- expand.grid(col_names, row_names)
contrasts$GenoA <- as.vector(contrasts$GenoA)
contrasts$GenoB <- as.vector(contrasts$GenoB)
colnames(contrasts) <- c("GenoA", "GenoB")

# Loop
int <- list()
for (i in 1:dim(contrasts)[1]) {
  a <- unlist(as.vector(sigProts[contrasts$GenoA[i]]))
  b <- unlist(as.vector(sigProts[contrasts$GenoB[i]]))
  int[[i]] <- intersect(a, b)
}

# Add to contrasts.
contrasts$Intersection <- unlist(lapply(int, function(x) length(x)))

# Make overlap matrix.
dm <- matrix(contrasts$Intersection, nrow = 8, ncol = 8)
rownames(dm) <- colnames(dm) <- row_names

# Remove upper tri and melt.
dm[lower.tri(dm)] <- NA

# Calculate percent overlap.
dm2 <- sweep(dm, 1, apply(dm, 1, function(x) max(x, na.rm = TRUE)), FUN = "/")

# Melt
df <- melt(dm, na.rm = TRUE)
df$percent <- round(melt(dm2, na.rm = TRUE)$value, 2)

# Generate plot
plot <- ggplot(df, aes(Var2, Var1, fill = percent)) +
  geom_tile(color = "black", size = 0.5) +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  scale_fill_gradient2(name = "Percent Overlap") +
  theme(
    # axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    # axis.title.y = element_text(color = "black", size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.7),
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    barwidth = 7, barheight = 1,
    title.position = "top", title.hjust = 0.5
  )) +
  coord_fixed()

# Store plot
all_plots[["TAMPOR_Condition_Overlap"]] <- plot

#--------------------------------------------------------------------------------
## Generate protein boxplots.
#-------------------------------------------------------------------------------

# Remove QC from traits. Group WT cortex and WT striatum.
out <- alltraits$SampleType == "QC"
traits <- alltraits[!out, ]
traits$Tissue.Sample.Model <- paste(traits$Tissue, traits$Sample.Model, sep = ".")
traits$Condition <- traits$Tissue.Sample.Model
traits$Condition[grepl("Cortex.WT", traits$Condition)] <- "Cortex.WT"
traits$Condition[grepl("Striatum.WT", traits$Condition)] <- "Striatum.WT"

# Levels for boxplots (order of the boxes):
lvls <- c(
  "Cortex.WT", "Striatum.WT",
  "Cortex.KO.Shank2", "Striatum.KO.Shank2",
  "Cortex.KO.Shank3", "Striatum.KO.Shank3",
  "Cortex.HET.Syngap1", "Striatum.HET.Syngap1",
  "Cortex.KO.Ube3a", "Striatum.KO.Ube3a"
)

# Generate plots.
plot_list <- ggplotProteinBoxPlot(
  data_in = log2(cleanDat),
  interesting.proteins = rownames(cleanDat),
  traits = traits,
  order = lvls,
  scatter = TRUE
)

# Add custom colors.
colors <- c(
  "gray", "gray",
  "#FFF200", "#FFF200",
  "#00A2E8", "#00A2E8",
  "#22B14C", "#22B14C",
  "#A349A4", "#A349A4"
)
plot_list <- lapply(plot_list, function(x) x + scale_fill_manual(values = colors))

## Add significance stars.
# Build a df with statistical results.
stats <- lapply(glm_results, function(x) data.frame(Uniprot = x$Uniprot, FDR = x$FDR))
stats <- stats %>% purrr::reduce(left_join, by = "Uniprot")
colnames(stats)[c(2:ncol(stats))] <- names(glm_results)
rownames(stats) <- stats$Uniprot
stats$Uniprot <- NULL

# Loop to add stars.
plot_list <- lapply(plot_list, function(x) annotate_stars(x, stats))

# Store boxplots.
all_plots[["all_box_plots"]] <- plot_list

# Top proteins.
p1 <- plot_list$`Shank2|Q80Z38`
p2 <- plot_list$`Shank3|Q4ACU6`
p3 <- plot_list$`Syngap1|F6SEU4`
p4 <- plot_list$`Ube3a|O08759`

# Modify x-axis labels-- significant bar axis labels are in red.
a <- c("black", "black", "red", "red", "black", "black", "black", "black", "black", "black")
b <- c("black", "black", "red", "black", "red", "red", "black", "black", "black", "red")
c <- c("black", "black", "black", "black", "black", "red", "red", "red", "black", "black")
d <- c("black", "black", "black", "black", "black", "black", "black", "red", "red", "red")
p1 <- p1 + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
p2 <- p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = b))
p3 <- p3 + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = c))
p4 <- p4 + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = d))

# Store plots in list.
all_plots[[paste(tissue, "Shank2_BP", sep = "_")]] <- p1
all_plots[[paste(tissue, "Shank3_BP", sep = "_")]] <- p2
all_plots[[paste(tissue, "Syngap1_BP", sep = "_")]] <- p3
all_plots[[paste(tissue, "Ube3a_BP", sep = "_")]] <- p4

#-------------------------------------------------------------------------------
## Write data to excel spreadsheet.
#-------------------------------------------------------------------------------

# Load data.
files <- list(
  traits = paste(Rdatadir, "2_Combined_traits.Rdata", sep = "/"),
  raw_cortex = paste(Rdatadir, "1_Cortex_raw_peptide.RData", sep = "/"),
  raw_striatum = paste(Rdatadir, "1_Striatum_raw_peptide.RData", sep = "/"),
  cleanDat = paste(Rdatadir, "2_Combined_cleanDat.RData", sep = "/")
)
data <- lapply(files, function(x) readRDS(x))

# Clean up traits.
traits <- data$traits
rownames(traits) <- NULL
colnames(traits)[1] <- "Batch.Channel"
traits$Color <- NULL
traits$Order <- NULL
colnames(traits)[2] <- "LongName"

# Gather raw data.
raw_cortex <- data$raw_cortex
raw_striatum <- data$raw_striatum
idx <- grepl("Abundance", colnames(raw_cortex))
colnames(raw_cortex)[idx] <- paste(colnames(raw_cortex)[idx], "Cortex", sep = ", ")
idx <- grepl("Abundance", colnames(raw_striatum))
colnames(raw_striatum)[idx] <- paste(colnames(raw_striatum)[idx], "Striatum", sep = ", ")

# Gather normalized data.
norm_data <- as.data.frame(log2(data$cleanDat))
idx <- match(colnames(norm_data), traits$Batch.Channel)
colnames(norm_data) <- paste(traits$LongName[idx], traits$Tissue[idx], sep = ", ")
norm_data <- add_column(norm_data, "Gene|Uniprot" = rownames(norm_data), .before = 1)
rownames(norm_data) <- NULL

# Write to excel workbook.
wb <- createWorkbook()
addWorksheet(wb, sheetName = "sample_info")
addWorksheet(wb, sheetName = "raw_cortex")
addWorksheet(wb, sheetName = "raw_striatum")
addWorksheet(wb, sheetName = "combined_normalized_data")
writeData(wb, sheet = 1, keepNA = TRUE, traits)
writeData(wb, sheet = 2, keepNA = TRUE, raw_cortex)
writeData(wb, sheet = 3, keepNA = TRUE, raw_striatum)
writeData(wb, sheet = 4, keepNA = TRUE, norm_data)
file <- paste(rootdir, "tables", "2_Combined_TMT_Data.xlsx", sep = "/")
saveWorkbook(wb, file, overwrite = TRUE)

# Save all plots.
myfile <- file.path(Rdatadir, paste0(outputMatName, "_plots.RData"))
saveRDS(all_plots, myfile)

message("Done!")
