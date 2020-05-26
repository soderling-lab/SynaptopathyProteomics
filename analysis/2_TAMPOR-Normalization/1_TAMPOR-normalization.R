#!/usr/bin/env Rscript

#' ---
#' title: 1_TAMPOR.R
#' description: TAMPOR Normalization of preprocessed TMT data.
#' authors: Tyler W Bradshaw, Eric B Dammer.
#' ---

## Parameters
save_plots = TRUE
clear_plots = TRUE
save_work = TRUE
fig_width = 2.5
fig_height = 2.5

## Prefix for output files.
output_name = "Combined"

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. 
# Load custom functions and prepare the project directory for saving 
# output files.

# Load renv.
renv::load(getrd(),quiet=TRUE)

# Load required packages.
suppressPackageStartupMessages({
  library(grid)
  library(dplyr)
  library(edgeR)
  library(purrr)
  library(tibble)
  library(gtable)
  library(ggplot2)
  library(openxlsx)
  library(reshape2)
  library(gridExtra)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

# Load additional functions in root/R.
devtools::load_all()

# Set any other directories.
root <- getrd()
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
fontdir <- file.path(root, "fonts")
tabsdir <- file.path(root, "tables")
figsdir <- file.path(root, "figs","TAMPOR")

# If necessary, create output directories for figs and tables.
if (!dir.exists(figsdir)) { dir.create(figsdir) }
if (!dir.exists(tabsdir)) { dir.create(tabsdir) }

# Remove any existing figures.
if (clear_plots) {
	invisible(sapply(list.files(figsdir,full.names=TRUE),unlink))
}

# Globally set ggplots theme.
ggtheme()
set_font("Arial",font_path=fontdir)

#---------------------------------------------------------------------
## Merge cortex and striatum data.
#---------------------------------------------------------------------
# We will utilize TAMPOR to combine the Cortex and Striatum datasets.
# Merge the preprocessed data and traits files.

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
myfiles <- file.path(rdatdir, c(
  "Cortex_cleanDat.RData",
  "Striatum_cleanDat.RData"
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
myfile <- file.path(rdatdir, "Combined_traits.RData")
saveRDS(alltraits, file = myfile)

#---------------------------------------------------------------------
## Perform TAMPOR normalization.
#---------------------------------------------------------------------

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

#---------------------------------------------------------------------
## Identify and remove any sample outliers.
#---------------------------------------------------------------------

message("\nExaming data for sample level outliers...")

# Remove QC samples from the data.
out <- colnames(cleanDat) %in% rownames(traits)[traits$SampleType == "QC"]
data_in <- log2(cleanDat[, !out])

# Illustrate Oldham's sample connectivity.
sample_connectivity <- ggplotSampleConnectivity(data_in, 
						log = TRUE, colID = "b.")
plot <- sample_connectivity$connectivityplot +
  ggtitle("Sample Connectivity post-TAMPOR")
plot <- plot + theme(panel.background=element_blank())
plot <- plot + theme(panel.border =  element_blank(), axis.line= element_line())
plot <- plot + scale_x_continuous(expand=c(0,0))
plot <- plot + scale_y_continuous(expand=c(0,0))

# Save.
if (save_plots) {
	myfile <- file.path(figsdir,"Sample_Outliers.pdf")
	ggsave(myfile,plot,height=fig_height,width=fig_width)
}

# Loop to identify Sample outliers using Oldham's connectivity method.
n_iter <- 5
threshold <- -3.0
out_samples <- list()

# Loop:
for (i in 1:n_iter) {
  data_temp <- data_in
  oldham <- ggplotSampleConnectivity(data_temp, log = TRUE, colID = "b", threshold = -3)
  bad_samples <- rownames(oldham$table)[oldham$table$Z.Ki < threshold]
  message(paste(length(bad_samples), " outlier sample(s) identified in iteration ", i, ".", sep = ""))
  if (length(bad_samples) == 0) bad_samples <- "none"
  out_samples[[i]] <- bad_samples
  out <- grepl(paste(unlist(out_samples), collapse = "|"), colnames(data_in))
  data_in <- data_in[, !out]
}

# Outlier samples.
bad_samples <- unlist(out_samples)
message(paste("\nOutlier sample(s) removed:", 
	      traits$Sample.Model[rownames(traits) %in% bad_samples]))

# Save data with QC samples, but outliers removed to file.
cleanDat <- cleanDat[, !colnames(cleanDat) %in% bad_samples]
myfile <- file.path(rdatdir, "Combined_cleanDat.RData")
saveRDS(cleanDat, myfile)

#---------------------------------------------------------------------
## Examine sample clustering post-TAMPOR Normalization.
#---------------------------------------------------------------------

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

plot1 <- ggplotPCA(log2(data_in[, idx]), traits, colID = "b.",
  colors[idx], title = "Cortex")
plot1 <- plot1 + theme(panel.background=element_blank())
plot1 <- plot1 + theme(panel.border =  element_rect(fill="NA"))
#plot1 <- plot1 + scale_x_continuous(expand=c(0,0))
#plot1 <- plot1 + scale_y_continuous(expand=c(0,0))

plot2 <- ggplotPCA(log2(data_in[, idy]), traits, colID = "b.",
  colors[idy], title = "Striatum")
plot2 <- plot2 + theme(panel.background=element_blank())
plot2 <- plot2 + theme(panel.border =  element_rect(fill="NA"))

plot3 <- ggplotPCA(log2(data_in), traits, colors, colID = "b.",
  title = "Combined")
plot3 <- plot3 + theme(panel.background=element_blank())
plot3 <- plot3 + theme(panel.border =  element_rect(fill="NA"))

# Save.
if (save_plots) {
	myfile <- file.path(figsdir,"Cortex_PCA.pdf")
	ggsave(myfile,plot1,height=fig_height,width=fig_width)
	myfile <- file.path(figsdir,"Striatum_PCA.pdf")
	ggsave(myfile,plot2,height=fig_height,width=fig_width)
	myfile <- file.path(figsdir,"Combined_PCA.pdf")
	ggsave(myfile,plot3,height=fig_height,width=fig_width)
}

#---------------------------------------------------------------------
## Create protein identifier map.
#---------------------------------------------------------------------

# Get uniprot ids from rownames.
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

# Map remainder of missing ids by hand.
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

# Check:
if (sum(is.na(gene)) > 0) {
  stop("Not all genes mapped to symbols!")
}

# Protein identifier map.
protmap <- data.frame(ids, uniprot, entrez, gene)

# Save to Rdata.
myfile <- file.path(rdatdir, "Protein_ID_Map.RData")
saveRDS(protmap, myfile)

#---------------------------------------------------------------------
## EdgeR statistical comparisons post-TAMPOR.
#---------------------------------------------------------------------

message("\nPerforming statistical testing with EdgeR glm...")

# Statistical comparisons are KO/HET versus all WT of a tissue type.

# Prepare data for EdgeR.
# Data should NOT be log2 transformed.
# Remove QC samples prior to passing data to EdgeR.
out <- alltraits$SampleType[match(colnames(cleanDat), 
				  rownames(alltraits))] == "QC"
data_in <- cleanDat[, !out]

# Number of proteins...
nprots <- formatC(dim(data_in)[1],big.mark=",")
nsamples <- dim(data_in)[2]
message(paste("\nQuantified", nprots, "proteins from", 
	      nsamples, "samples."))

# Create DGEList object...
y_DGE <- DGEList(counts = data_in)

# TMM Normalization.
y_DGE <- calcNormFactors(y_DGE)

# Create sample mapping to Tissue.Genotype.
# Group WT Cortex samples and WT Striatum samples together.
traits <- subset(alltraits, rownames(traits) %in% colnames(data_in))
traits <- traits[match(colnames(data_in), rownames(traits)), ]
if (!all(traits$SampleID == colnames(data_in))) { stop() }
group <- paste(traits$Tissue, traits$Sample.Model, sep = ".")
group[grepl("Cortex.WT", group)] <- "Cortex.WT"
group[grepl("Striatum.WT", group)] <- "Striatum.WT"
traits$group <- group
y_DGE$samples$group <- as.factor(group)

# Basic design matrix for GLM -- all groups treated seperately.
design <- model.matrix(~ 0 + group, data = y_DGE$samples)
colnames(design) <- levels(y_DGE$samples$group)

# Estimate dispersion:
y_DGE <- estimateDisp(y_DGE, design, robust = TRUE)

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

# Pretty summary:
message("\nSummary of differentially abundant proteins (FDR<0.1):")
knitr::kable(overall,row.names=FALSE)

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

# Check the table.
plot <- cowplot::plot_grid(mytable)

# Save.
if (save_plots) {
	myfile <- file.path(figsdir,"Sig_Prots_Summary.pdf")
	ggsaveTable(mytable,myfile)
}

# Call topTags to add FDR. Gather tabularized results.
glm_results <- lapply(qlf, function(x){ 
			      topTags(x, n = Inf, sort.by = "none")$table })

# Convert logCPM column to percent WT and annotate with candidate column.
glm_results <- lapply(glm_results, function(x) annotateTopTags(x))

# Use protmap to annotate glm_results with entrez Ids and gene symbols.
for (i in 1:length(glm_results)) {
  x <- glm_results[[i]]
  idx <- match(rownames(x), protmap$ids)
  x <- add_column(x, "Gene|Uniprot" = protmap$ids[idx], .before = 1)
  x <- add_column(x, "Uniprot" = protmap$uniprot[idx], .after = 1)
  x <- add_column(x, "Entrez" = protmap$entrez[idx], .after = 2)
  x <- add_column(x, "Symbol" = protmap$gene[idx], .after = 3)
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
myfile <- file.path(rdatdir,paste0(output_name,"_GLM_Results.RData"))
saveRDS(glm_results, myfile)

# Save results to file as spreadsheet.
myfile <- file.path(tabsdir,paste0(output_name,"_GLM_Results.xlsx"))
write_excel(glm_results, myfile)

#---------------------------------------------------------------------
## Collect GLM statistics in a list.
#---------------------------------------------------------------------

# Names of relevant columns.
colNames <- colnames(glm_results[[1]])[c(2, 5:9)]
stats <- colNames[c(2:length(colNames))]

# Collect data.
subGLM <- lapply(glm_results, function(x) x[, colNames])

# Combine into a single df.
df <- subGLM %>% purrr::reduce(left_join, by = "Uniprot")

# Rename columns.
newNames <- paste(
  rep(names(glm_results), each = length(colNames) - 1),
  sapply(strsplit(colnames(df)[c(2:ncol(df))], "\\."), "[", 1)
)
colnames(df)[c(2:ncol(df))] <- newNames

# Collect each statistic into a single df in a list.
glm_stats <- sapply(stats, function(x) df[, c(1, grep(x, colnames(df)))])

# Clean up data a little...
glm_stats <- lapply(glm_stats, function(x) {
  x <- x[order(x$Uniprot), ]
  idx <- match(x$Uniprot, protmap$uniprot)
  rownames(x) <- protmap$ids[idx]
  x$Uniprot <- NULL
  return(x)
})

# Save GLM stats.
myfile <- file.path(rdatdir, "GLM_Stats.RData")
saveRDS(glm_stats, myfile)

#---------------------------------------------------------------------
## Volcano plots for each genotype.
#---------------------------------------------------------------------

message("\nGenerating volcano plots...")

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

# Generate plots.
# FIXME: not working!
plots <- lapply(as.list(names(colors)), function(x) { 
			ggplotVolcanoPlot(results[[x]])
})
names(plots) <- names(colors)

# Add titles.
for (i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + ggtitle(names(plots)[i])
}

# Save.
if (save_plots) {
	myfile <- file.path(figsdir,"Shank2_Volcano.pdf")
	ggsave(myfile,plots$Shank2,height=fig_height,width=fig_width)
	myfile <- file.path(figsdir,"Shank3_Volcano.pdf")
	ggsave(myfile,plots$Shank3,height=fig_height,width=fig_width)
	myfile <- file.path(figsdir,"Syngap1_Volcano.pdf")
	ggsave(myfile,plots$Syngap1,height=fig_height,width=fig_width)
	myfile <- file.path(figsdir,"Ube3a_Volcano.pdf")
	ggsave(myfile,plots$Ube3a,height=fig_height,width=fig_width)
}

#---------------------------------------------------------------------
## Condition overlap plot.
#---------------------------------------------------------------------

message("\nGenerating condition overlap plot...")

# Load statistical results..
results <- glm_results

# Combine by FDR.
stats <- lapply(results, function(x) {
  data.frame(Protein = x$Gene, FDR = x$FDR)
})
names(stats) <- names(results)
df <- stats %>% purrr::reduce(left_join, by = c("Protein"))
colnames(df)[c(2:ncol(df))] <- paste("FDR", names(stats), sep = ".")

# Data frame of stats to be used to annotate plots.
stats_df <- df
rownames(stats_df) <- stats_df$Protein
stats_df$Protein <- NULL

# Gather sigProts.
sigProts <- list()
for (i in c(2:ncol(df))) {
  idx <- df[, i] < 0.05
  sigProts[[i-1]] <- df$Protein[idx]
}
names(sigProts) <- names(stats)
all_sigProts <- unique(unlist(sigProts))

# sigProts by tissue type:
idx <- grep("Cortex",names(sigProts))
idy <- grep("Striatum",names(sigProts))
tissue_sigProts <- list("Cortex" = unique(unlist(sigProts[idx])),
			"Striatum" = unique(unlist(sigProts[idy])))

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
df <- melt(dm, na.rm = FALSE)
df$percent <- round(melt(dm2)$value, 2)
df <- df[!is.na(df$value) | !is.na(df$percent),]
df$percent[is.na(df$percent)] <- 0

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

# Save.
if (save_plots) {
	myfile <- file.path(figsdir,"Condition_Overlap_Plot.pdf")
	ggsave(myfile,plot,height=fig_height,width=fig_width)
}

#---------------------------------------------------------------------
## Generate faceted protein boxplots.
#---------------------------------------------------------------------

message("\nGenerating faceted protein plots...")

# Remove QC from traits. Group WT cortex and WT striatum.
out <- alltraits$SampleType == "QC"
traits <- alltraits[!out, ]
traits$Tissue.Sample.Model <- paste(traits$Tissue, 
				    traits$Sample.Model, sep = ".")
traits$Condition <- traits$Tissue.Sample.Model
traits$Condition[grepl("Cortex.WT", traits$Condition)] <- "Cortex.WT"
traits$Condition[grepl("Striatum.WT", traits$Condition)] <- "Striatum.WT"

# Levels for boxplots (order of the boxes):
box_order <- c(
  "Cortex.WT",
  "Cortex.KO.Shank2",
  "Cortex.KO.Shank3",
  "Cortex.HET.Syngap1",
  "Cortex.KO.Ube3a",
  "Striatum.WT",
  "Striatum.KO.Shank2",
  "Striatum.KO.Shank3",
  "Striatum.HET.Syngap1",
  "Striatum.KO.Ube3a"
)

# Generate all faceted plots.
# This takes a couple minutes.
plot_list <- ggplotProteinBoxPlot(
  data_in = log2(cleanDat),
  interesting.proteins = rownames(cleanDat),
  traits = traits,
  box_order,
  scatter = TRUE
)

# Add custom colors.
colors <- rep(c("gray", "#FFF200", "#00A2E8", "#22B14C", "#A349A4"), 2)
plot_list <- lapply(plot_list, function(x) x + 
		    scale_fill_manual(values = colors))

# Facet plots, add significance stars, and reformat x.axis labels.
plot_list <- lapply(plot_list, function(x) annotate_plot(x, stats_df))

# Collect significant plots.
# any significance in both tissues!
prots <- Reduce(intersect,tissue_sigProts)
plot_list <- plot_list[prots]

# Customization.
plot_list <- lapply(plot_list, function(plot) {
	       plot <- plot + theme(panel.background=element_blank())
	       plot <- plot + theme(panel.border =  element_blank(), 
				    axis.line= element_line())
	       return(plot)
		    })

# Save sig plots.
plotdir <- file.path(figsdir,"Faceted-Protein-Boxplots")

# Create directory if necessary.
if (!dir.exists(plotdir)) { dir.create(plotdir) }

# Remove existing plots.
if (clear_plots) {
	invisible(sapply(list.files(plotdir,full.names=TRUE),unlink))
}

# Loop to save plots.
if (save_plots) {
	message("\nSaving plots, this will take several minutes...")
	for (i in 1:length(plot_list)) {
		namen <- gsub("\\|","_",names(plot_list)[i])
		myfile <- file.path(plotdir,paste0(namen,".pdf"))
		ggsave(myfile,plot_list[[i]],
		       height=fig_height,width=fig_width)
	}
}

#---------------------------------------------------------------------
## Save Cortex plots seperately.
#---------------------------------------------------------------------

message("\nGenerating Cortex protein plots...")

# Remove QC from traits. Group WT cortex and WT striatum.
out <- alltraits$SampleType == "QC"
traits <- alltraits[!out, ]
traits$Tissue.Sample.Model <- paste(traits$Tissue, 
				    traits$Sample.Model, sep = ".")
traits$Condition <- traits$Tissue.Sample.Model
traits$Condition[grepl("Cortex.WT", traits$Condition)] <- "Cortex.WT"
traits$Condition[grepl("Striatum.WT", traits$Condition)] <- "Striatum.WT"

# Levels for boxplots (order of the boxes):
box_order <- c(
  "Cortex.WT",
  "Cortex.KO.Shank2",
  "Cortex.KO.Shank3",
  "Cortex.HET.Syngap1",
  "Cortex.KO.Ube3a"
)

# Generate plots.
plot_list <- ggplotProteinBoxPlot(
  data_in = log2(cleanDat),
  interesting.proteins = rownames(cleanDat),
  traits = traits,
  box_order,
  scatter = TRUE
)

# Add custom colors.
colors <- c("gray", "#FFF200", "#00A2E8", "#22B14C", "#A349A4")
plot_list <- lapply(plot_list, function(x) x + 
		    scale_fill_manual(values = colors))

# Add significance stars, and reformat x.axis labels.
plot_list <- lapply(plot_list, function(x) annotate_plot(x, stats_df))

# Collect significant plots.
sigCortex <- unique(unlist(sigProts[grep("Cortex",names(sigProts))]))
plot_list <- plot_list[sigCortex]

# Additional customization.
plot_list <- lapply(plot_list, function(plot) {
		df <- plot$data
		df <- plot$data %>% 
			group_by(Group) %>% 
			summarize(Median = median(Intensity))
		wt_median <- df$Median[grepl("WT",df$Group)]
		plot <- plot +
			geom_hline(yintercept=wt_median,
				   linetype="dotted",colour="black")
	       plot <- plot + theme(panel.background=element_blank())
	       plot <- plot + theme(panel.border =  element_blank(), 
				    axis.line= element_line())
	       return(plot)
		    })

# Save sig plots.
plotdir <- file.path(figsdir,"Cortex-Protein-Boxplots")

# Create directory if necessary.
if (!dir.exists(plotdir)) { dir.create(plotdir) }

# Remove existing plots.
if (clear_plots) {
	invisible(sapply(list.files(plotdir,full.names=TRUE),unlink))
}

# Save sig cortex plots.
if (save_plots) {
	message("\nSaving plots, this will take several minutes...")
	for (i in 1:length(plot_list)) {
		namen <- gsub("\\|","_",names(plot_list)[i])
		myfile <- file.path(plotdir,paste0(namen,".pdf"))
		ggsave(myfile,plot_list[[i]],
		       height=fig_height,width=fig_width)
	}
}

#---------------------------------------------------------------------
## Save Striatum plots seperately.
#---------------------------------------------------------------------

message("\nGenerating Striatum protein plots...")

# Remove QC from traits. Group WT cortex and WT striatum.
out <- alltraits$SampleType == "QC"
traits <- alltraits[!out, ]
traits$Tissue.Sample.Model <- paste(traits$Tissue, 
				    traits$Sample.Model, sep = ".")
traits$Condition <- traits$Tissue.Sample.Model
traits$Condition[grepl("Cortex.WT", traits$Condition)] <- "Cortex.WT"
traits$Condition[grepl("Striatum.WT", traits$Condition)] <- "Striatum.WT"

# Levels for boxplots (order of the boxes):
box_order <- c(
  "Striatum.WT",
  "Striatum.KO.Shank2",
  "Striatum.KO.Shank3",
  "Striatum.HET.Syngap1",
  "Striatum.KO.Ube3a"
)

# Generate plots.
plot_list <- ggplotProteinBoxPlot(
  data_in = log2(cleanDat),
  interesting.proteins = rownames(cleanDat),
  traits = traits,
  box_order,
  scatter = TRUE
)

# Add custom colors.
colors <- c("gray", "#FFF200", "#00A2E8", "#22B14C", "#A349A4")
plot_list <- lapply(plot_list, function(x) x + 
		    scale_fill_manual(values = colors))

# Add significance stars, and reformat x.axis labels.
plot_list <- lapply(plot_list, function(x) annotate_plot(x, stats_df))

# Collect significant plots.
sigStriatum <- unique(unlist(sigProts[grep("Striatum",names(sigProts))]))
plot_list <- plot_list[sigStriatum]

# Custumization.
plot_list <- lapply(plot_list, function(plot) {
	       plot <- plot + theme(panel.background=element_blank())
	       plot <- plot + theme(panel.border =  element_blank(), 
				    axis.line= element_line())
	       return(plot)
		    })

# Save sig plots.
plotdir <- file.path(figsdir,"Striatum-Protein-Boxplots")

# Create directory if necessary.
if (!dir.exists(plotdir)) { dir.create(plotdir) }

# Remove existing plots.
if (clear_plots) {
	invisible(sapply(list.files(plotdir,full.names=TRUE),unlink))
}

# Save sig striatum plots.
if (save_plots) {
	message("\nSaving plots, this will take several minutes...")
	for (i in 1:length(plot_list)) {
		namen <- gsub("\\|","_",names(plot_list)[i])
		myfile <- file.path(plotdir,paste0(namen,".pdf"))
		ggsave(myfile,plot_list[[i]],
		       height=fig_height,width=fig_width)
	}
}

#---------------------------------------------------------------------
## Write data to excel spreadsheet.
#---------------------------------------------------------------------

message("\nSaving data to file!")

# Load data.
files <- list(
  traits = paste(rdatdir, "Combined_traits.RData", sep = "/"),
  raw_cortex = paste(rdatdir, "Cortex_raw_peptide.RData", sep = "/"),
  raw_striatum = paste(rdatdir, "Striatum_raw_peptide.RData", sep = "/"),
  cleanDat = paste(rdatdir, "Combined_cleanDat.RData", sep = "/")
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
colnames(raw_cortex)[idx] <- paste(colnames(raw_cortex)[idx], 
				   "Cortex", sep = ", ")
idx <- grepl("Abundance", colnames(raw_striatum))
colnames(raw_striatum)[idx] <- paste(colnames(raw_striatum)[idx], 
				     "Striatum", sep = ", ")

# Gather normalized data.
norm_data <- as.data.frame(log2(data$cleanDat))
idx <- match(colnames(norm_data), traits$Batch.Channel)
colnames(norm_data) <- paste(traits$LongName[idx], 
			     traits$Tissue[idx], sep = ", ")
norm_data <- add_column(norm_data, 
			"Gene|Uniprot" = rownames(norm_data), .before = 1)
rownames(norm_data) <- NULL

# Write to excel workbook.
wb <- createWorkbook()
addWorksheet(wb, sheetName = "sample-info")
addWorksheet(wb, sheetName = "raw-cortex-peptide")
addWorksheet(wb, sheetName = "raw-striatum-peptide")
addWorksheet(wb, sheetName = "combined-normalized-protein")
writeData(wb, sheet = 1, keepNA = TRUE, traits)
writeData(wb, sheet = 2, keepNA = TRUE, raw_cortex)
writeData(wb, sheet = 3, keepNA = TRUE, raw_striatum)
writeData(wb, sheet = 4, keepNA = TRUE, norm_data)
myfile <- file.path(tabsdir, paste0(output_name,"_TMT_Data.xlsx"))
saveWorkbook(wb, myfile, overwrite = TRUE)

# Save the workspace?
if (save_work) {
	save(list = ls(all.names = TRUE), 
	     file=paste0(output_name,".RData"), envir=.GlobalEnv)
}

# Complete!
message("\nDone!")
