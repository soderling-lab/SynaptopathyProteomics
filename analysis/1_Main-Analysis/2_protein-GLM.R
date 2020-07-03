#!/usr/bin/env Rscript

#' ---
#' title:
#' description: Preprocessing of TMT data.
#' authors: Tyler W Bradshaw
#' ---

## Optional Parameters:

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
  library(dplyr)
  library(edgeR)
  library(data.table)
  library(tibble)
})

# Load project specific data and functions.
suppressWarnings({ devtools::load_all() })

# Project directories.
datadir <- file.path(root,"data") # Input/Output data.
rdatdir <- file.path(root,"rdata") # Temporary data files.

# Load the data.
data(cortex) # tidy_prot

# Cast tp into data matrix for EdgeR. Don't log!
dm <- tidy_prot %>% 
	dcast(Accession ~ Sample, value.var="Intensity") %>% 
	as.matrix(rownames=TRUE)

#---------------------------------------------------------------------
## Create protein identifier map.
#---------------------------------------------------------------------

# Create gene map.
uniprot <- rownames(dm)
entrez <- mgi_batch_query(ids=uniprot)
symbols <- getPPIs::getIDs(entrez,from="entrez",to="symbol",species="mouse")
gene_map <- as.data.table(keep.rownames="uniprot",entrez)
gene_map$symbol <- symbols[as.character(gene_map$entrez)]

# Save as rda.
myfile <- file.path(datadir,"gene_map.rda")
save(gene_map,file=myfile,version=2)

#--------------------------------------------------------------------
## EdgeR glm
#--------------------------------------------------------------------

# Create dge object.
dge <- DGEList(counts=dm)

# Perform TMM normalization.
dge <- calcNormFactors(dge)

samples <- unique(rownames(dge$samples))
replicate <- as.numeric(as.factor(samples))
idx <- match(samples,tidy_prot$Sample)
genotype <- tidy_prot$Genotype[idx]
treatment <- tidy_prot$Treatment[idx]
group <- paste(genotype,treatment,sep = ".")
dge$samples$group <- as.factor(group)

# Basic design matrix for GLM -- all groups treated seperately.
design <- model.matrix(~ 0 + group, data = dge$samples)
colnames(design) <- levels(dge$samples$group)

# Estimate dispersion:
dge <- estimateDisp(dge, design, robust = TRUE)

# Fit a general linear model.
fit <- glmQLFit(dge, design, robust = TRUE)

# Generate contrasts.
g <- colnames(design)
g <- g[!grepl("WT",g)]
contrasts <- paste(g,paste0(sapply(strsplit(g,"\\."),"[",1),".WT"),sep="-")

# Make contrasts for EdgeR.
# For some reason loops or lapply dont work with the makeContrasts function.

contrast_list <- list(
  makeContrasts(contrasts[1], levels = design),
  makeContrasts(contrasts[2], levels = design),
  makeContrasts(contrasts[3], levels = design),
  makeContrasts(contrasts[4], levels = design))

names(contrasts) <- unlist({
  lapply(contrasts, function(x) sapply(strsplit(colnames(x), " "), "[", 1))
})

# Call glmQLFTest() to evaluate differences in contrasts.
contrasts<-colnames(design)
qlf <- lapply(contrasts, function(x) glmQLFTest(fit, coef = x))
names(qlf) <- contrasts

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


# Create sample groupings given contrasts of interest.
colnames(tidy_prot)
dt <- tidy_prot %>% dplyr::select(Sample,Genotype,Treatment) %>% unique() 
idx <- match(rownames(dge$samples),dt$Sample)
genotype <- dt$Genotype[idx]
treatment <- dt$Treatment[idx]
dge$samples$group <- as.character(interaction(genotype,treatment))[idx]

# Create a design matrix for GLM.
design <- model.matrix(~ 0 + group, data = dge$samples)

# Estimate dispersion.
dge <- estimateDisp(dge, design, robust = TRUE)

# Fit a general linear model.
fit <- glmQLFit(dge, design, robust = TRUE)

# Evaluate differences within genotypes.
genotypes = unique(genotype)
qlf_list <- list()
for (geno in genotypes){
	qlf <-  glmQLFTest(fit,coef=grep(geno,colnames(design)))
	qlf_list[[geno]] <- qlf
}

# Call topTags to add FDR. Gather tabularized results.
glm_results <- lapply(qlf_list, function(x) {
			      topTags(x, n = Inf, sort.by = "PValue")$table })

# Insure first column is Accession.
glm_results <- lapply(glm_results,function(x) {
			      as.data.table(x,keep.rownames="Accession") })

# Add percent WT and sort by pvalue.
glm_results <- lapply(glm_results,function(x) {
			      x[["Percent Control"]] <- 2^x$logFC
			      return(x) })


# Add gene annotations.
results_list <- list()
for (i in c(1:length(glm_results))) {
	df <- glm_results[[i]]
	idx <- match(df$Accession,gene_map$uniprot)
	Entrez <- gene_map$entrez[idx]
	Symbol <- gene_map$symbol[idx]
	df <- tibble::add_column(df,Entrez,.after=1)
	df <- tibble::add_column(df,Symbol,.after=2)
	results_list[[i]] <- df
}

names(results_list) <- gsub("\\(|\\)","",colnames(design))

# Save.
myfile <- file.path(root,"tables",paste0(tissue,"_GLM_Results.xlsx"))
write_excel(results_list,file=myfile)

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
# Save results to file.
myfile <- file.path(rdatdir,paste0(output_name,"_GLM_Results.RData"))
saveRDS(glm_results, myfile)

# Save results to file as spreadsheet.
myfile <- file.path(tabsdir,paste0(output_name,"_TMT_Results.xlsx"))
write_excel(glm_results, myfile)
