GO enrichment analaysis for DEPs from each genotype.
#-------------------------------------------------------------------------------

# Build a df with the combined statistical results.
stats <- lapply(glm_results, function(x) data.frame(Uniprot = x$Uniprot, FDR = x$FDR))
names(stats) <- names(glm_results)
df <- stats %>% purrr::reduce(left_join, by = "Uniprot")
colnames(df)[c(2:ncol(df))] <- names(stats)

# Save to Rdata.
myfile <- file.path(outputtabs,"2_GLM_Results.RData")
saveRDS(df, myfile)

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
