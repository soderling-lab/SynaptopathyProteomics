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
myfile <- file.path(datadir, "combined_protein.rda")
combined_protein <- cleanDat
save(combined_protein, file = myfile, version = 2)

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

# Check the table.
plot <- cowplot::plot_grid(mytable)

# Save.
if (save_plots) {
	myfile <- file.path(figsdir,"Sig_Prots_Summary.pdf")
	ggsaveTable(mytable,myfile)
}


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
