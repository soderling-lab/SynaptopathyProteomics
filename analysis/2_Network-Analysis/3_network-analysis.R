#!/usr/bin/env Rscript

#' ---
#' title: Network Analysis
#' description: Protein co-expression network analysis
#' authors: Tyler W Bradshaw
#' ---

## User parameters to change:
analysis_type = "Cortex"
root = "/mnt/d/projects/SynaptopathyProteomics" 

## Input data in root/rdata:
input_data <- list("Cortex" = list(
				   adjm = "Cortex_Adjm.csv",
				   netw = "Cortex_NE_Adjm.csv",
				   gmap = "Cortex_gene_map.RData",
			   	   data = "Cortex_norm_protein.csv",
				   part = "Cortex_NE_SurpriseVertexPartition.csv",
				   part = "Cortex_partition_self_preservation_enforced.csv"),
		   "Striatum" = list(
				     adjm = "Striatum_Adjm.csv",
				     netw = "Striatum_NE_Adjm.csv",
				     gmap = "Striatum_gene_map.RData",
				     data = "Striatum_norm_protein.csv",
				     part = "Striatum_NE_SurpriseVertexPartition.csv",
				     pres = "Striatum_partition_self_preservation_enforced.csv")
		   )[[analysis_type]]

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

# Global imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Functions.
TBmiscr::load_all()

# Directories.
root <- TBmiscr::load_all()
rdatdir <- file.path(root, "rdata")

# Load expression data:
# Load the data, subset, coerce to matrix, Log2 transform, and 
# finally transpose such that rows = samples and columns = proteins.
myfile <- file.path(rdatdir, input_data[['data']])
dm <- fread(myfile) %>%
	dcast(Accession ~ Sample,value.var="Intensity") %>%
	as.matrix(rownames="Accession") %>% log2() %>% t()

# Load adjmatrix--coerce to a matrix.
myfile <- file.path(rdatdir, input_data[['adjm']])
adjm <- fread(myfile) %>% as.matrix(rownames="Accession")

# Load network--coerce to a matrix.
myfile <- file.path(rdatdir, input_data[['netw']])
netw <- fread(myfile) %>% as.matrix(rownames="Accession")

# Load Leidenalg graph partition.
myfile <- file.path(rdatdir, input_data[['part']])
part_dt <- fread(myfile, drop=1)
resolutions <- nrow(part_dt)

# Load graph partition after enforcing module self-preservation.
myfile <- file.path(rdatdir, input_data[['pres']])
part_dt <- fread(myfile)
partition <- as.numeric(part_dt)

# Load gene identifier map.
gene_map <- readRDS(file.path(rdatdir, input_data[['gmap']]))

#---------------------------------------------------------------------------------
# Get proteins with any significant change.
idy <- lapply(c("Cortex", "Striatum"), function(x) {
  grep(x, colnames(glm_stats$FDR))
})
sigProts <- lapply(idy, function(x) {
  apply(glm_stats$FDR[, x], 1, function(pval) {
    any(pval < 0.05)
  })
})
names(sigProts) <- c("Cortex", "Striatum")
sigProts[["Combined"]] <- apply(cbind(sigProts[[1]], sigProts[[2]]), 1, any)

# Load expression data.
# Data should be transposed: rows, proteins.
# Insure that column names are not lost after transpose.
myfile <- file.path(rdatdir, input_files$data_file[[data_type]])
temp <- readRDS(myfile)
data <- t(temp)
colnames(data) <- rownames(temp)

# Load Sample info.
sampleTraits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load co-expression (adjacency) matrix.
myfile <- file.path(rdatdir, input_files$adjm_file[[data_type]])
adjm <- as.matrix(readRDS(myfile))
rownames(adjm) <- colnames(adjm)

# Load PPI adjacency matrix.
adjm_ppi <- fread(file.path(rdatdir, "3_PPI_Adjm.csv"), drop = 1)
adjm_ppi <- as.matrix(adjm_ppi)
rownames(adjm_ppi) <- colnames(adjm_ppi)

# Load network partitions--self-preservation enforced.
myfiles <- file.path(rdatdir, input_files$part_file[[part_type]])
partitions <- lapply(myfiles, function(x) unlist(readRDS(x)))
names(partitions) <- names(input_files$part_file[[part_type]])

# Load theme for plots.
ggtheme()

# Use arial font.
set_font("Arial")

# Load theme for tables.
# FIXME: ADD THIS! Include arial font.
# mytheme <- gtabtheme()

# Reset partition index for self-preserved modules.
partitions$self <- reset_index(partitions$self)

# Define opposite tissue type.
if (part_type == "Cortex") {
  other <- "Striatum"
} else if (part_type == "Striatum") {
  other <- "Cortex"
}

#---------------------------------------------------------------------
## SigProt annotations--which genotype is a given protein changing in?
#---------------------------------------------------------------------

# Create a df of protein-SigProt annotations.
df <- glm_stats$FDR
colnames(df) <- gsub(" FDR", "", colnames(df))
df$sigProt <- apply(df, 1, function(x) {
  paste(colnames(df)[x < 0.05], collapse = "; ")
})
df$sigProt[df$sigProt == ""] <- NA
sigProtAnno <- df$sigProt
names(sigProtAnno) <- rownames(df)

#---------------------------------------------------------------------
## Collect all modules in a list.
#---------------------------------------------------------------------

# Create list of modules.
module_list <- list()
partition <- partitions$self

# Module list with entrez ids.
idx <- match(names(partition), protmap$ids)
module_list[["Entrez"]] <- split(protmap$entrez[idx], partition)

# Module list with gene symbols.
module_list[["Symbols"]] <- split(protmap$gene[idx], partition)

# Module list with protein ids.
module_list[["IDs"]] <- split(partition, partition)

# Name modules with M prefix.
module_list <- lapply(module_list, function(x) {
  names(x) <- paste0("M", names(x))
  return(x)
})

# Drop M0 from lists.
module_list <- lapply(module_list, function(x) x[-which(names(x) == "M0")])

#---------------------------------------------------------------------
## Which modules are preserved in opposite tissue?
#---------------------------------------------------------------------

# Load partitions.
p1 <- partitions$self
p2 <- partitions$other

# Get preserved modules.
preserved_modules <- list()
modules <- split(p2, p1)
names(modules) <- paste0("M", names(modules))
preserved <- which(sapply(modules, function(x) unique(x) != 0))
preserved_modules[["other"]] <- names(modules)[preserved]
otherSig <- preserved_modules[["other"]]

# Fraction of modules that are preserved in the other tissue type.
nModules <- length(unique(p1[p1 != 0]))
nPres <- length(preserved_modules[["other"]])
pPres <- round(100 * (nPres / length(unique(p1))), 3)
message(paste0(
  nPres, " of ", nModules, " (", pPres, "%) ", part_type,
  " modules are preserved in the opposite tissue."
))

# Reformat other partition, so that module names are the same
# as in self-preserved partition (partitions$self)!
m <- split(p1, p2)
new_part <- sapply(seq_along(m), function(x) {
  p <- rep(as.numeric(names(m)[x]), length(m[[x]]))
  names(p) <- names(m[[x]])
  return(p)
})
partitions$other <- unlist(new_part)

#---------------------------------------------------------------------
## Which modules are preserved in the PPI graph?
#---------------------------------------------------------------------

# Load partitions.
p1 <- partitions$self
p2 <- partitions$ppi

# Get preserved modules.
modules <- split(p2, p1)
preserved <- which(sapply(modules, function(x) unique(x) != 0))
preserved_modules[["ppi"]] <- paste0("M", names(modules)[preserved])
PPIsig <- preserved_modules[["ppi"]]

# Fraction of modules that are preserved in the other tissue type.
nModules <- length(unique(p1[p1 != 0]))
nPres <- length(preserved_modules[["ppi"]])
pPres <- round(100 * (nPres / length(unique(p1))), 3)
message(paste0(
  nPres, " of ", nModules, " (", pPres, "%) ", part_type,
  " modules are preserved in the PPI network."
))

#---------------------------------------------------------------------
## Which modules are enriched for DBD-associated genes?
#---------------------------------------------------------------------

# Load Disease ontology.
DBDset <- "mouse_Combined_DBD_collection.RData"
myfile <- file.path(rdatdir, DBDset)
DBDcollection <- readRDS(myfile)

# Perform disease enrichment analysis.
gene_list <- module_list$Entrez
DBDenrichment <- gse(gene_list, DBDcollection)

# Collect modules with significant enrichment of DBD-genes.
method <- "FDR"
alpha <- 0.05
any_sig <- which(sapply(DBDenrichment, function(x) any(x[method] < alpha)))
DBDsig <- names(DBDenrichment)[any_sig]

# Add DBD genes and terms to list of preserved modules.
# Function to reformat pvalues and combine with term nam.e
fx <- function(x) {
  return(formatC(x, digits = 2, format = "e"))
}
term_pval <- sapply(DBDenrichment[DBDsig], function(x) {
  paste0(
    x$shortDataSetName[x$FDR < 0.05],
    " (FDR = ", fx(x$FDR[x$FDR < 0.05]), ")"
  )
})
preserved_modules$dbd <- term_pval

# Summary of Modules with DBD enrichment:
# term_pval[DBDsig]

# Status.
message(paste(
  "Total number of DBD-associated modules:",
  length(DBDsig)
))

# Write to file.
myfile <- file.path(tabsdir, "3_Module_DBD_Enrichment.xlsx")
write_excel(DBDenrichment[DBDsig], myfile)

# Collect all DBD genes.
DBDgenes <- lapply(DBDcollection$dataSets, function(x) {
  as.character(x$data$Entrez)
})
names(DBDgenes) <- sapply(DBDcollection$dataSets, function(x) x$name)

# All DBDprots.
DBDprots <- lapply(DBDgenes, function(x) {
  protmap$ids[which(x %in% protmap$entrez)]
})

# Create a vector of protein-DBD annotations.
DBDcols <- do.call(cbind, lapply(DBDprots, function(x) protmap$ids %in% x))
colnames(DBDcols) <- names(DBDprots)
DBDdf <- as.data.frame(DBDcols)
rownames(DBDdf) <- protmap$ids
DBDanno <- apply(DBDdf, 1, function(x) paste(colnames(DBDdf)[x], collapse = "; "))
DBDanno[DBDanno == ""] <- NA

# DBD colors:
DBD_colors <- c(
  "ASD" = "yellow",
  "ID" = "blue",
  "ADHD" = "orange",
  "Schizophrenia" = "red",
  "Bipolar disorder" = "green",
  "Epilepsy" = "purple"
)

# Function to mix colors.
mix_colors <- function(x) {
  idy <- as.logical(x)
  if (any(idy)) {
    # If any sig, mix colors.
    dm <- sapply(DBD_colors[idy], col2rgb)
    dm <- apply(dm, 1, mean)
    color <- rgb(dm[1], dm[2], dm[3], maxColorValue = 255)
    return(color)
  } else {
    # Else return gray.
    return("#808080")
  }
}

# DBD colors:
DBDcolors <- sapply(c(1:nrow(DBDdf)), function(x) mix_colors(DBDdf[x, ]))
names(DBDcolors) <- rownames(DBDdf)

#---------------------------------------------------------------------
## Is the Synaptosome enriched for DBD genes?
#---------------------------------------------------------------------


#---------------------------------------------------------------------
## Create a table summarizing disease genes.
#---------------------------------------------------------------------

# Summarize disease genes in final clustered dataset.
part <- partitions$self
prots <- names(part)[which(part != 0)]

# Filter DBDprots
myprots <- sapply(DBDprots, function(x) x[x %in% prots])

# Summarize all DBD genes identified in synaptosome.
df <- t(as.data.frame(sapply(myprots, length)))
rownames(df) <- NULL

# Theme for tables.
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params = list(cex = 0.75)),
  colhead = list(fg_params = list(cex = 0.75)),
  rowhead = list(fg_params = list(cex = 0.75))
)

# Create table and add borders.
mytable <- tableGrob(df, rows = NULL, theme = mytheme)

# Add Border around data rows.
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 1, b = nrow(mytable), l = 1, r = ncol(mytable)
)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 2, l = 1, r = ncol(mytable)
)

# Check the table.
fig <- plot_grid(mytable)

# Save.
# Note: ggssaveTable uses arial font.
myfile <- filename("DBD_Gene_Summary", fig_ext, figsdir)
ggsaveTable(mytable, myfile)

#--------------------------------------------------------------------
## Which modules are enriched for GO terms?
#--------------------------------------------------------------------

## Use the anRichment package.
# Build a GO collection.
message("Performing GO analysis with anRichment...")
GOcollection <- buildGOcollection(organism = "mouse")

# Perform gene set enrichment analysis.
GOresults <- gse(gene_list, GOcollection)

# Significant go terms for every module.
method <- "FDR"
alpha <- 0.05
topGO <- lapply(GOresults, function(x) {
  idx <- x[[method]] < alpha
  p <- x[[method]][idx]
  names(p) <- x$shortDataSetName[idx]
  return(p)
})

# Number of modules with any significant GO term enrichment.
GOsig <- names(which(sapply(topGO, function(x) length(x) > 0)))
message(paste(
  "Total number of modules with any significant GO",
  "enrichment:", length(GOsig)
))

# Add to list of preserved modules.
term_pval <- sapply(topGO, function(x) {
  paste0(names(x), " (FDR = ", fx(x), ")")
})
preserved_modules$go <- term_pval

#--------------------------------------------------------------------
## Which modules are enriched for ASD DEGs?
#--------------------------------------------------------------------

# Load ASD DEG collection.
GeneSet <- "mouse_ASD_DEG_collection.RData"
myfile <- file.path(rdatdir, GeneSet)
ASDcollection <- readRDS(myfile)

# All ASD DEGs
ASD_DEGs <- unique(ASDcollection$dataSets[[1]]$data$Entrez)
idx <- match(ASD_DEGs[which(ASD_DEGs %in% protmap$entrez)], protmap$entrez)
ASDprots <- protmap$ids[idx]

# Perform gene set enrichment analysis.
gene_list <- module_list$Entrez
ASDresults <- gse(gene_list, ASDcollection)

# Significant go terms for every module.
method <- "FDR"
alpha <- 0.05
topASD <- lapply(ASDresults, function(x) {
  idx <- x[[method]] < alpha & x$enrichmentRatio > 1
  p <- x[[method]][idx]
  names(p) <- x$shortDataSetName[idx]
  return(p)
})

# Number of modules with any significant GO term enrichment.
ASDsig <- names(which(sapply(topASD, function(x) length(x) > 0)))
message(paste(
  "Total number of modules with significant enrichment",
  "of ASD DEGs:", length(ASDsig)
))

# Add to list of preserved modules.
term_pval <- sapply(topASD[ASDsig], function(x) {
  paste0(names(x), " (FDR = ", fx(x), ")")
})
preserved_modules$asd <- term_pval

#---------------------------------------------------------------------
## Explore changes in module summary expression.
#---------------------------------------------------------------------

# Total Number of modules.
# P.values will be corrected for n comparisions.
all_modules <- module_list$IDs
nModules <- sum(names(all_modules) != "M0")
message(paste("Number of modules:", nModules))

# Get partition of interest.
# If analyzing the combined data, focus on modules that are preserved
# in opposite tissue type.
if (data_type == "Combined") {
  partition <- partitions$other
} else {
  # Otherwise focus on modules that are self-preserved.
  partition <- partitions$self
}

# Modules from the given partition.
modules <- split(partition, partition)
names(modules) <- paste0("M", names(modules))

# If not working with combined data,
# drop modules that are preserved in other tissue.
# These will be analyzed seperately.
if (data_type != "Combined") {
  out <- which(names(modules) %in% preserved_modules$other)
  modules <- modules[-out]
}

# Remove M0.
modules <- modules[-which(names(modules) == "M0")]

# Fix partition--only analyzing modules defined above.
out <- partition %notin% as.numeric(gsub("M", "", names(modules)))
partition[out] <- 0

# Calculate Module Eigengenes.
# Note: Soft power does not influence MEs.
# Note: Do not need to sort partition to be in the same order!
MEdata <- moduleEigengenes(data,
  colors = partition,
  excludeGrey = TRUE, # Ignore M0!
  softPower = 1,
  impute = FALSE
)
MEs <- as.matrix(MEdata$eigengenes)

# Create list of MEs.
# Do it this way to preserve names.
ME_list <- lapply(seq(ncol(MEs)), function(x) MEs[, x])
names(ME_list) <- names(modules)

# Module membership (KME).
KMEdata <- signedKME(data, MEdata$eigengenes,
  corFnc = "bicor",
  outputColumnName = "M"
)
KME_list <- lapply(seq(ncol(KMEdata)), function(x) {
  v <- vector("numeric", length = nrow(KMEdata))
  names(v) <- rownames(KMEdata)
  v[] <- KMEdata[[x]]
  v <- v[order(v, decreasing = TRUE)]
  return(v)
})
names(KME_list) <- colnames(KMEdata)

# Get Percent Variance explained (PVE).
PVE <- as.numeric(MEdata$varExplained)
names(PVE) <- names(modules)
medianPVE <- median(PVE)
message(paste(
  "Median module coherence (PVE):",
  round(100 * medianPVE, 2), "(%)."
))

# Sample to group mapping for statistical testing.
sampleTraits$Sample.Model.Tissue <- paste(sampleTraits$Sample.Model,
  sampleTraits$Tissue,
  sep = "."
)
idx <- match(rownames(MEs), sampleTraits$SampleID)
groups <- sampleTraits$Sample.Model.Tissue[idx]
names(groups) <- rownames(MEs)

# Group all WT samples from a tissue type together.
groups[grepl("WT.*.Cortex", groups)] <- "WT.Cortex"
groups[grepl("WT.*.Striatum", groups)] <- "WT.Striatum"

# If working with combined data, group all WT samples together...
if (data_type == "Combined") {
  # WT from both tissue will be annotated as WT.
  groups[grepl("WT", groups)] <- "WT"
}

# Perform Kruskal Wallis tests to identify modules whose summary
# expression profile is changing.
KWdata_list <- lapply(ME_list, function(x) {
  kruskal.test(x ~ groups[names(x)])
})
KWdata <- as.data.frame(do.call(rbind, KWdata_list))[-c(4, 5)]

# Correct p-values for n comparisons.
KWdata$p.adj <- as.numeric(KWdata$p.value) * nModules
KWdata$p.adj[KWdata$p.adj > 1.0] <- 1.0

# Significant modules.
alpha <- 0.1
sigModules <- rownames(KWdata)[KWdata$p.adj < alpha]
nSigModules <- length(sigModules)
message(paste0(
  "Number of modules with significant (p.adj < ", alpha, ")",
  " Kruskal-Wallis test: ", nSigModules, "."
))

# Perform Dunnetts test for post-hoc comparisons.
# Note: P-values returned by DunnettTest have already been adjusted for
# multiple comparisons!

# Define control group and levels (order) for DunnettTest.
group_order <- c("WT", "KO.Shank2", "KO.Shank3", "HET.Syngap1", "KO.Ube3a")
group_levels <- c(
  paste(group_order, "Cortex", sep = "."),
  paste(group_order, "Striatum", sep = ".")
)

# If combining data, set all WT samples to WT.
if (data_type == "Combined") {
  control_group <- "WT"
  group_levels[grepl("WT", group_levels)] <- "WT"
  group_levels <- unique(group_levels)
} else {
  control_group <- paste("WT", part_type, sep = ".")
}

# Loop to perform DTest. NOTE: This takes several seconds.
DTdata_list <- lapply(ME_list, function(x) {
  g <- factor(groups[names(x)], levels = group_levels)
  result <- DunnettTest(x ~ g, control = control_group)[[control_group]]
  return(as.data.frame(result))
})

# Number of modules with significant KW + DT changes.
DTalpha <- 0.05
nSigDT <- sapply(DTdata_list, function(x) sum(x$pval < DTalpha))
if (length(nSigDT[sigModules]) > 0) {
  message("Summary of Dunnett's test changes for significant modules:")
  print(nSigDT[sigModules])
}

# Numer of significant modules with any DBD association.
sigDBDmodules <- nSigDT[sigModules[which(sigModules %in% DBDsig)]]
nSigDisease <- length(sigDBDmodules)
if (nSigDisease > 0) {
  message(paste(
    "Number of significant modules with",
    "significant enrichment of DBD-associated genes:",
    nSigDisease
  ))
  message(paste(
    "Summary of Dunnett's test changes for significant,",
    "DBD-associated modules:"
  ))
  print(sigDBDmodules)
}

# Save significant modules.
myfile <- file.path(rdatdir, paste0(
  data_type, "_", part_type,
  "_sigModules.RData"
))
saveRDS(sigModules, myfile)

#---------------------------------------------------------------------
## Generate verbose boxplots.
#---------------------------------------------------------------------

# Reset sample to group mapping for generating boxplots.
sampleTraits$Sample.Model.Tissue <- paste(sampleTraits$Sample.Model,
  sampleTraits$Tissue,
  sep = "."
)
idx <- match(rownames(MEs), sampleTraits$SampleID)
groups <- sampleTraits$Sample.Model.Tissue[idx]
names(groups) <- rownames(MEs)

# Group all WT samples from a tissue type together.
groups[grepl("WT.*.Cortex", groups)] <- "WT.Cortex"
groups[grepl("WT.*.Striatum", groups)] <- "WT.Striatum"

# Reset group order.
group_order <- c("WT", "KO.Shank2", "KO.Shank3", "HET.Syngap1", "KO.Ube3a")
group_levels <- c(
  paste(group_order, "Cortex", sep = "."),
  paste(group_order, "Striatum", sep = ".")
)

# Generate boxplots summarizing module protein expression.
if (data_type != "Combined") {
  control_group <- NULL
}
plots <- lapply(ME_list, function(x) {
  ggplotVerboseBoxplot(
    x, groups,
    group_levels, control_group
  )
})
names(plots) <- names(ME_list)

## Loop to  clean-up plots.
# Loop to clean-up plots.
for (k in seq_along(plots)) {
  # Simplify x-axis labels.
  x_labels <- rep(c(
    "WT", "Shank2 KO", "Shank3 KO",
    "Syngap1 HET", "Ube3a KO"
  ), 2)
  # Add title and fix xlabels.
  plot <- plots[[k]]
  m <- names(plots)[k]
  # Add significance stars!
  DTdata_list[[m]]$pval
  dt <- DTdata_list[[m]]
  df <- data.table(plot$data)
  if (data_type != "Combined") {
    control_group <- paste("WT", part_type, sep = ".")
  }
  idx <- match(paste(df$g, control_group, sep = "-"), rownames(dt))
  df$pval <- dt$pval[idx]
  df$pval[is.na(df$pval)] <- 1
  df$xpos <- df$g
  df$ypos <- 1.1 * max(df$x)
  df$symbol <- ""
  df$symbol[df$pval < 0.05] <- "*"
  df$symbol[df$pval < 0.005] <- "**"
  df$symbol[df$pval < 0.0005] <- "***"
  df$tissue <- factor(df$tissue, levels = c("Cortex", "Striatum"))
  df$color <- "black"
  df$color[df$pval < 0.05] <- "red"
  # If KW Sig, then add stars.
  if (KWdata[m, "p.adj"] < 0.1) {
    plot <- plot +
      geom_text(data = df, aes(
        x = xpos,
        y = ypos, label = symbol, size = 7
      ))
  }
  # Fix title and xlabels.
  txt <- paste0(
    "P.adj = ", round(KWdata[m, "p.adj"], 3),
    "; ", "PVE = ", round(PVE[m], 3)
  )
  plot_title <- paste0(m, " (", txt, ")")
  plot$labels$title <- plot_title
  plot <- plot + scale_x_discrete(labels = x_labels)
  # Store results in list.
  plots[[k]] <- plot
} # Ends loop to fix plots.

#---------------------------------------------------------------------
## Save verbose box plots.
#---------------------------------------------------------------------

# Save all modules.
for (i in seq_along(plots)) {
  namen <- names(plots)[i]
  myfile <- filename(namen, fig_ext, modsdir)
  ggsave(myfile, plots[[i]],
    height = 7, width = 7)
}

# Save sig modules as single pdf.
# myfile <- filename("Sig_Module_Boxplots","pdf",figsdir)
# ggsavePDF(plots[sigModules],myfile)

#---------------------------------------------------------------------
## Save sigProt boxplots for sig modules.
#---------------------------------------------------------------------

# Load protein boxplots.
# NOTE: The plots are only of proteins with any sig change.
# If Cortex or striatum -- then only cortex or striatum are plotted.
# If Combined -- then data from both tissues are plotted.
myfiles <- c(
  Cortex = file.path(
    rdatdir,
    "All_Cortex_SigProt_Boxplots.RData"
  ),
  Striatum = file.path(
    rdatdir,
    "All_Striatum_SigProt_Boxplots.RData"
  ),
  Combined = file.path(
    rdatdir,
    "All_Faceted_SigProt_Boxplots.RData"
  )
)
myfile <- myfiles[data_type]
all_plots <- readRDS(myfile)

# Group by module.
plot_list <- split(all_plots, partition[names(all_plots)])
names(plot_list) <- paste0("M", names(plot_list))

# Save a single pdf for each sig module containing all of its
# sig proteins.
## TODO: suppress output from grid.arrange!
# for (i in 1:length(sigModules)){
# 	module_name <- sigModules[i]
# 	plots <- plot_list[[module_name]]
# 	n <- length(plots)
# 	groups <- rep(c(1:ceiling(n/4)),each=4)[c(1:n)]
# 	plot_groups <- split(plots,groups)
# 	figs <- lapply(plot_groups,function(x) {
# 			       fig <- grid.arrange(grobs=x,ncol=2,nrow=2)
# 			       return(fig)
# 		 })
# 	# Save.
# 	namen <- paste0(module_name,"_Module_SigProts")
# 	myfile <- filename(namen,"pdf",figsdir)
# 	ggsavePDF(figs,myfile)
# 	# Close the device.
# 	if (i == length(sigModules)) { dev.off() }
# }

# Save protein boxplots for sig modules.
for (module in sigModules) {
  plots <- plot_list[[module]]
  namen <- paste(module, names(plots), sep = "_")
  for (i in seq_along(plots)) {
    myfile <- filename(namen[i], fig_ext, protdir)
    ggsave(myfile, plots[[i]],
      height = 7, width = 7)
  }
}

#---------------------------------------------------------------------
## Examine overall structure of network.
#---------------------------------------------------------------------

# Examine relationships between modules by comparing their summary
# expression profiles (MEs). This is done for tissue specific MEs--
# i.e. do not use the Combined data.

# Recalculate MEs.
myfile <- file.path(rdatdir, input_files$data_file[[part_type]])
temp <- readRDS(myfile)
tempdat <- t(temp)
colnames(tempdat) <- rownames(temp)
MEdata <- moduleEigengenes(tempdat,
  colors = partitions$self,
  excludeGrey = TRUE,
  softPower = 1,
  impute = FALSE
)
MEs <- as.matrix(MEdata$eigengenes)
ME_list <- lapply(seq(ncol(MEs)), function(x) MEs[, x])
names(ME_list) <- names(all_modules)

# Get PVE.
pve <- MEdata$varExplained

# Calculate correlations between ME vectors.
adjm_me <- cor(do.call(cbind, ME_list))

## Perform Heirarchical clustering of the ME adjm.
# Convert similarity matrix to distance matrix, and cluster with hclust.
hc <- hclust(as.dist(1 - adjm_me), method = "ward.D2")

# Optimize cut of tree with modularity.
# First remove any negative edges from graph.
g <- graph_from_adjacency_matrix(adjm_me, mode = "undirected", weighted = TRUE)
g <- delete_edges(g, E(g)[E(g)$weight < 0])

# Generate a bunch of partitions by cutting the tree at various heights.
h <- seq(0, max(hc$height), by = 0.005)
hc_partitions <- lapply(h, function(x) cutree(hc, h = x))
# Number of groups, k.
k <- sapply(hc_partitions, function(x) length(unique(x)))
# Modularity, q.
q <- sapply(hc_partitions, function(x) {
  modularity(g, x, weights = abs(edge_attr(g, "weight")))
})

# Find the best cut height--the cut height that maximizes modularity.
best_q <- unique(q[seq(h)[q == max(q)]])
best_h <- median(h[seq(h)[q == max(q)]])
best_k <- unique(k[seq(k)[q == max(q)]])
message(paste0("Cut height that produces the best partition: ", best_h, "."))
message(paste0("Number of groups: ", best_k, " (Modularity = ", round(best_q, 3), ")."))

## With the optimal k, cut the graph.
hc_partition <- cutree(hc, k = best_k)
groups <- split(hc_partition, hc_partition)

# Get representative module from each group, its medoid.
# The medoid is the module which is closes (i.e. most similar)
# to all others in its group.
# rep_modules <- getMedoid(adjm_me,k=best_k)

# Try getting hub instead...
rep_modules <- sapply(groups, function(x) {
  idx <- idy <- names(x)
  col_sums <- apply(adjm_me[idx, idy], 2, sum)
  idmax <- names(which(col_sums == max(col_sums)))
  return(idmax)
})

# Get dendrogram data, update with group and rep_module.
dend_data <- ggdendro::dendro_data(as.dendrogram(hc))
dend_data <- dend_data$labels
dend_data$group <- hc_partition[as.character(dend_data$label)]
dend_data$rep_module <- dend_data$label %in% rep_modules

# Add actual position to dend data.
dend_data$x <- match(as.character(dend_data$label), hc$label[hc$order])
dend_data <- dend_data %>% arrange(x)
dend_data$label <- factor(dend_data$label, levels = dend_data$label)

# Generate dendrogram.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE, labels = FALSE) +
  geom_hline(yintercept = best_h, color = "red", size = 1)
dendro <- dendro +
  geom_text(
    data = dend_data,
    aes(x, y, label = label, color = rep_module),
    hjust = 1, angle = 90, size = 3
  ) +
  scale_colour_manual(values = c("black", "red")) +
  theme(legend.position = "none")

#---------------------------------------------------------------------
## Create a table summarizing key network stats.
#---------------------------------------------------------------------

# If working with the Combined data, get tissue specific sigModules,
# and add these to sigModules.
if (data_type == "Combined") {
  myfile <- file.path(rdatdir, paste0(
    part_type, "_", part_type,
    "_sigModules.RData"
  ))
  tissueSig <- readRDS(myfile)
  AllsigModules <- c(tissueSig, sigModules)
} else {
  AllsigModules <- sigModules
}

# Summarize key network stats.
nClust <- sum(partitions$self != 0)
nNodes <- length(partition)
pClust <- round(100 * nClust / nNodes, 2)
nodes <- paste0(nClust, "/", nNodes, " (", pClust, "%)")
medPVE <- round(median(as.numeric(pve)), 3)
nModules <- length(all_modules)
nSigModules <- length(AllsigModules)

# Collect in df.
df <- as.data.table(cbind(
  c(
    "Algorithm", "Modularity", "N Modules",
    "N Clustered", "Median PVE", "N SigModules"
  ),
  c(
    "Leiden", "Surprise", nModules, nodes, medPVE,
    nSigModules
  )
))

# Create table and add borders.
mytable <- tableGrob(df, cols = NULL, rows = NULL, theme = mytheme)
# Add Border around data rows.
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 1, b = nrow(mytable), l = 1, r = ncol(mytable)
)

# Check the table.
fig <- plot_grid(mytable)

# Save.
myfile <- filename("Network_Summary", fig_ext, figsdir)
ggsaveTable(mytable, myfile)

#---------------------------------------------------------------------
## Generate module colors based on their cor with rep modules.
#---------------------------------------------------------------------

# Assign modules a color based on similarity with rep modules.
df <- do.call(cbind, lapply(rep_modules, function(x) adjm_me[, x]))
colnames(df) <- paste0(rep_modules, "cor")
# Remove negative correlations.
# df <- df*as.numeric(df>0)

# Use RGB colors for three groups!
# Rescale the data to [0,1].
dm <- (df - min(df)) / (max(df) - min(df))
colnames(dm) <- c("R", "G", "B")
df <- data.table(cbind(df, dm))
df <- tibble::add_column(df, Module = colnames(adjm_me), .before = 1)
df <- df %>%
  arrange(match(df$Module, as.character(dend_data$label)))
# Convert RGB to hexadecimal color.
df$color <- rgb(255 * df$R, 255 * df$G, 255 * df$B,
  maxColorValue = 255
)

# Collect color assignments.
module_colors <- df$color
names(module_colors) <- names(all_modules)

# Add gray.
module_colors <- c(module_colors, c("M0" = col2hex("light gray")))

# Collect protein color assignments.
node_colors <- module_colors[paste0("M", partitions$self)]
names(node_colors) <- names(partitions$self)

# Generate colored bars for dendrogram.
dend_data$color <- module_colors[as.character(dend_data$label)]
dend_data$color <- factor(dend_data$color, levels = dend_data$color)
p2 <- ggplot(dend_data, aes(x = label, y = 1, fill = color)) +
  geom_tile() +
  scale_fill_manual(values = as.character(dend_data$color)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45))

# Save.
myfile <- filename("Modules_Dendro", fig_ext, figsdir)
ggsave(myfile,
  plot = dendro, height = 3, width = 3)

# Save.
myfile <- filename("Module_Colors", fig_ext, figsdir)
ggsave(myfile,
  plot = p2, height = 3, width = 3)

#---------------------------------------------------------------------
## Generate GO Scatter plots.
#---------------------------------------------------------------------

# Get subset of GO data.
GOresults <- GOresults[names(modules)]

# Loop to generate plots:
plots <- list()
alpha <- 0.05
for (i in 1:length(GOresults)) {
  # Get subset of data.
  df <- GOresults[[i]]
  namen <- names(GOresults)[i]
  if (nrow(df) <= 10) {
    topN <- nrow(df) - 1
  } else {
    topN <- 10
  }
  # Generate the plot.
  plot <- ggplotGOscatter(df, module_colors[namen], topN)
  # Add approximate significance threshold.
  if (any(df$FDR < alpha)) {
    threshold <- floor(-log(df$pValue[max(which(df$FDR < 0.05))]))
    plot <- plot +
      geom_hline(
        yintercept = threshold,
        linetype = "dashed", color = "red"
      )
  }
  # Add title.
  plot <- plot + ggtitle(namen) +
    theme(plot.title = element_text(color = "black", size = 14))
  plots[[i]] <- plot
}
names(plots) <- names(modules)

# Save a single pdf with all plots.
# myfile <- filename("Module_GOscatter","pdf",figsdir)
# ggsavePDF(plots,myfile)

# Loop to save plots.
for (i in seq_along(plots)) {
  plot <- plots[[i]]
  namen <- names(plots)[i]
  myfile <- filename(namen, fig_ext, scatdir)
  ggsave(myfile, plot,
    height = 7, width = 7)
}

#---------------------------------------------------------------------
## Generate PPI and co-expression graphs.
#---------------------------------------------------------------------

## NOTE: Coerce boolean attributes to integer to avoid warnings when
# loading into cytoscape.

# Create co-expression graph.
exp_graph <- graph_from_adjacency_matrix(adjm,
  mode = "undirected",
  weighted = TRUE, diag = FALSE
)

# Create PPI graph.
ppi_graph <- graph_from_adjacency_matrix(adjm_ppi,
  mode = "undirected",
  weighted = TRUE, diag = FALSE
)

# Remove NAs from PPI edges.
E(ppi_graph)$weight[which(is.na(E(ppi_graph)$weight))] <- 0

## Add attributes to igraph object.
# Add Gene symbols.
symbols <- protmap$gene[match(names(V(exp_graph)), protmap$ids)]
exp_graph <- set_vertex_attr(exp_graph, "symbol", value = symbols)

# Add sigProt vertex attribute.
anySig <- as.numeric(sigProts[[data_type]][names(V(exp_graph))])
exp_graph <- set_vertex_attr(exp_graph, "sigProt",
  value = anySig
)

# Add any DBDprot vertex attribute.
DBDnodes <- lapply(DBDprots, function(x) names(V(exp_graph)) %in% x)
for (DBD in names(DBDnodes)) {
  exp_graph <- set_vertex_attr(exp_graph,
    name = DBD,
    value = as.numeric(DBDnodes[[DBD]])
  )
}

# Collect PPI evidence and save to file.
myfile <- file.path(rdatdir, "3_All_PPIs.RData")
ppis <- readRDS(myfile)

# Map mouse entrez to protein ids.
ppis$ProteinA <- protmap$ids[match(ppis$osEntrezA, protmap$entrez)]
ppis$ProteinB <- protmap$ids[match(ppis$osEntrezB, protmap$entrez)]
out <- is.na(ppis$ProteinA) | is.na(ppis$ProteinB)
ppis <- ppis[!out, ]

# Get the relevant columns.
ppis <- ppis %>% select(
  ProteinA, ProteinB, osEntrezA, osEntrezB,
  Interactor_A_Taxonomy, Interactor_B_Taxonomy,
  Source_database, Confidence_score,
  Publications, Methods
)

# Save to file.
myfile <- file.path(tabsdir, paste0("3_All_PPIs.xlxs"))
write_excel(list(PPIs = ppis), myfile)

#---------------------------------------------------------------------
## Demonstrate that interacting proteins are highly co-expressed.
#---------------------------------------------------------------------

# Merge PPI and co-expression graphs into a single data.table.
df <- as_long_data_frame(igraph::union(ppi_graph, exp_graph))
df <- as.data.table(df)

# Get the relevant data columns.
df <- df %>% dplyr::select(c(from_name, to_name, weight_1, weight_2))

# Convert NA's to 0/FALSE (PPI edges).
df$weight_1[is.na(df$weight_1)] <- 0

# Seed seed for reproducibility.
set.seed(1207)

# Randomly sample edges.
n <- 1000

# Get random subset of interacting proteins.
idx <- sample(which(df$weight_1 == 1), n)

# Get random subset of non-interacting proteins.
idy <- sample(which(df$weight_1 == 0), n)

# Subset the data, coerce PPI column to factor.
subdat <- df[c(idx, idy), ]
subdat$weight_1 <- factor(subdat$weight_1, levels = c(0, 1))

# Calculate WRS p-value.
# Refactor, test that TRUE > FALSE.
WRS_test <- wilcox.test(subdat$weight_2 ~ subdat$weight_1, alternative = "less")
WRS_pval <- formatC(WRS_test$p.value, digits = 2, format = "e")

# Generate a plot.
plot <- ggplot(subdat, aes(x = weight_1, y = weight_2, fill = weight_1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 20, outlier.size = 1) +
  scale_x_discrete(labels = c("PPI = False", "PPI = True")) +
  ylab("Protein co-expression\n(bicor correlation)") +
  xlab(NULL) +
  scale_fill_manual(values = c("gray", "dark red")) +
  ggtitle(part_type) +
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold"),
    axis.text.x = element_text(color = "black", size = 11, face = "bold"),
    legend.position = "none"
  )

# Add Annotation.
plot <- plot +
  annotate("text",
    x = 1.5, y = 1.0,
    label = paste("p-value =", WRS_pval), size = 6, color = "black"
  )

# Save figure.
myfile <- filename("WRS_PPI_Bicor_Proteins", fig_ext, figsdir)
ggsave(myfile, plot,
  height = 4, width = 4)

#---------------------------------------------------------------------
## Generate cytoscape graph of all clustered proteins.
#---------------------------------------------------------------------

# Prompt user to open Cytoscape if it is not already open.

# If working with Combined data, append graphs to tissue specific
cytoscape_ping()
# Cytoscape file.
if (data_type == "Combined") {
  cysfile <- file.path(netsdir, paste0(part_type, ".cys"))
  if (file.exists(cysfile)) {
    message(paste("Adding graphs to", part_type, "file!"))
    winfile <- gsub("/mnt/d/", "D:/", cysfile)
    openSession(winfile)
  } else {
    message(paste(
      "Analyze", part_type,
      "data first. Combined graphs will be",
      "appended to this file."
    ))
  }
}

# Define threshold for Cortex and Striatum networks.
threshold <- c("Cortex" = 7.815, "Striatum" = 6.0)
network_layout <- "force-directed edgeAttribute=weight"

# Create a cytoscape graph of the co-expression network.
# Perform network enhancement and thresholding to improve layout.
net <- createCytoscapeCoExpressionGraph(partitions$self,
  node_colors, DBDcolors,
  threshold = threshold[part_type],
  background = 0,
  network_layout,
  title = part_type
)

## Generate subnetworks highlighting some modules of interest.
# Significant proteins. Need to remove any that were not clustered.
all_nodes <- getAllNodes(network = net)
nodes <- sigProts[[data_type]]
sigNodes <- names(nodes[nodes])
highlightNodes(
  nodes = sigNodes[sigNodes %in% all_nodes],
  main.network = net, subnetwork.name = paste(part_type, "SigProts")
)

# DBDprots
for (i in 1:length(DBDprots)) {
  namen <- names(DBDprots)[i]
  nodes <- DBDprots[[i]]
  highlightNodes(
    nodes = nodes[nodes %in% all_nodes],
    main.network = net, subnetwork.name = namen,
    new.values = DBD_colors[namen]
  )
}

#---------------------------------------------------------------------
## Generate cytoscape graph summarizing overall topolgy of the network.
#---------------------------------------------------------------------

# Create a cytoscape network showing all modules.
if (create_networks) {
  net <- createCytoscapeModuleGraph(partitions$self, ME_list,
    title = paste(data_type, "Modules")
  )
}

#---------------------------------------------------------------------
## Highlight some important modules.
#---------------------------------------------------------------------

## Highlight modules that are:
#     Meta modules.
#     Significant KW test.
#     Enriched for DBDs.
#     Enriched for ASD DEGs.
#     Underlying PPI toplogy is significant.
#     Significant GO enrichment.
#     Preserved in other tissue.

# Generate subnetworks highlighting some modules of intereest.
highlightNodes(
  nodes = AllsigModules, main.network = net,
  subnetwork.name = paste(part_type, "Sig")
)
highlightNodes(
  nodes = DBDsig[DBDsig %in% names(all_modules)],
  main.network = net, subnetwork.name = "DBD"
)
highlightNodes(
  nodes = ASDsig[ASDsig %in% names(all_modules)],
  main.network = net, subnetwork.name = "ASD DEG"
)
highlightNodes(
  nodes = GOsig[GOsig %in% names(all_modules)],
  main.network = net, subnetwork.name = "GO"
)
highlightNodes(
  nodes = PPIsig[PPIsig %in% names(all_modules)],
  main.network = net, subnetwork.name = "PPI"
)
highlightNodes(
  nodes = rep_modules[rep_modules %in% names(all_modules)],
  main.network = net, subnetwork.name = "MetaModule Hubs"
)

# Meta modules groups.
for (i in 1:length(groups)) {
  highlightNodes(
    nodes = names(groups[[i]]), main.network = net,
    subnetwork.name = paste("MetaModule", i)
  )
}

#---------------------------------------------------------------------
## Generate cytoscape graphs of all co-expression modules.
#---------------------------------------------------------------------

# Directory for output network images.
imgdir <- file.path(
  strsplit(figsdir, paste0("/", subdir))[[1]][1], subdir,
  "Networks", part_type
)

# Loop to create graphs.
if (create_networks) {
  for (i in c(1:length(modules))) {
    module_name <- names(modules)[i]
    message(paste("Working on module", module_name, "..."))
    nodes <- names(modules[[module_name]])
    module_kme <- KME_list[[module_name]]
    network_layout <- "force-directed edgeAttribute=weight"
    image_file <- file.path(imgdir, module_name)
    image_format <- "SVG"
    createCytoscapeGraph(
      exp_graph, ppi_graph, nodes,
      module_kme, module_name,
      module_colors, network_layout,
      output_file, image_file,
      image_format
    )
    # When done, save cytoscape session.
    if (i == length(modules)) {
      myfile <- file.path(netsdir, paste0(part_type, ".cys"))
      winfile <- gsub("/mnt/d/", "D:/", myfile)
      saveSession(winfile)
    }
  } # ENDS loop.
} # ENDS Network chunk.

#---------------------------------------------------------------------
## Save key results summarizing modules.
#---------------------------------------------------------------------

# Summarize modules: Name, N Nodes, PVE, Color, Prots.
module_summary <- data.frame(
  Module = names(PVE),
  Nodes = sapply(modules, length),
  PVE = PVE
)

# Combine with KWdata.
tempKW <- KWdata
colnames(tempKW) <- paste("KW", colnames(tempKW))
module_summary <- cbind(module_summary, tempKW)

# Remove parameter column.
module_summary$"KW parameter" <- NULL

# Add column for KW sig.
module_summary$"KW sig" <- module_summary$"KW p.adj" < 0.1

# Combine with DT results.
reformatDT <- function(x) {
  df <- as.data.table(x, keep.rownames = TRUE) %>% select(rn, diff, pval)
  df <- melt(df, id.vars = "rn")
  values <- df$value
  names(values) <- paste("DT", sapply(strsplit(df$rn, "-"), "[", 1), df$variable)
  return(values)
}
dm <- do.call(rbind, lapply(DTdata_list, reformatDT))
module_summary <- cbind(module_summary, dm)

# Number of sig changes.
module_summary$"N DT sig" <- nSigDT

# Any DT sig.
module_summary$"Any DT sig" <- module_summary$"N DT sig" > 0

# Other protein attributes: PPI pres, tissue pres, OMIM, DBD, GOsig.
module_summary$PPI <- as.character(module_summary$Module) %in% preserved_modules$ppi
module_summary$DBD <- sapply(preserved_modules$dbd, function(x) paste(x, collapse = " | "))[as.character(module_summary$Module)]
module_summary$OMIM <- sapply(preserved_modules$omim, function(x) paste(x, collapse = " | "))[as.character(module_summary$Module)]
module_summary$PFAM <- sapply(preserved_modules$pfam, function(x) paste(x, collapse = " | "))[as.character(module_summary$Module)]
module_summary$ASD_DEGs <- sapply(preserved_modules$asd, function(x) paste(x, collapse = " | "))[as.character(module_summary$Module)]
module_summary$"Preserved in opposite tissue" <- sapply(preserved_modules$other, function(x) paste(x, collapse = " | "))[as.character(module_summary$Module)]

## More detailed summary of every module.
dfs <- lapply(seq_along(KME_list), function(x) {
  df <- data.table(Protein = names(partition))
  df$Module <- paste0("M", partition)
  df$KME <- KME_list[[x]][df$Protein]
  df$sigProt <- df$Protein %in% sigProts
  df$"Sig Genotype(s)" <- sigProtAnno[df$Protein]
  df$DBDProt <- df$Protein %in% DBDprots
  df$"DBD Association(s)" <- DBDanno[df$Protein]
  df <- df %>% filter(Module == names(KME_list)[x])
  df <- df[order(df$KME, decreasing = TRUE), ]
})
names(dfs) <- names(modules)

# Write to file.
results <- list()
results[["Summary"]] <- module_summary
results <- c(results, dfs[sigModules])
# Insure that if working with combined data we know which
# data and partition type...
if (data_type == "Combined") {
  myfile <- file.path(
    tabsdir,
    paste0(
      "3_", data_type, "_", part_type,
      "_Module_Summary.xlsx"
    )
  )
} else {
  myfile <- file.path(
    tabsdir,
    paste0("3_", data_type, "_Module_Summary.xlsx")
  )
}
write_excel(results, myfile)

# Remove that pesky Rplots.
unlink("Rplots.pdf")
