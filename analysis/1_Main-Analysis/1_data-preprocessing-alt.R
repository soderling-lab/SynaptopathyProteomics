#!/usr/bin/env Rscript

#' ---
#' title:
#' description: Preprocessing of TMT data.
#' authors: Tyler W Bradshaw
#' ---

## Optional Parameters:
oldham_threshold = -2.5 # Threshold for detecting sample level outliers.
FDR_threshold = 0.05 # Threshold for protein differential abundance.

#---------------------------------------------------------------------
## Overview of Data Preprocessing:
#---------------------------------------------------------------------

# [INTRA-Experiment processing]
#     Peptide-level operations:
#     | * SL normalization
#     | * Remove QC outliers
#     | * Impute missing values 
#     Protein-level proessing:
#     | * Summarize proteins
#     | * SL normalization
#     | * Remove intra-batch batch-effect.
# [INTER-Experiment operations]
#     | * Remove QC sample outliers
#     | * IRS normalization
#     | * Filter proteins
#     | * Protein imputing
# [edgeR STATISTICAL analysis]
#       * Final TMM normalization
#       * GLM for differential abundance.

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
parse_args <- function(default="Striatum", args=commandArgs(trailingOnly=TRUE)){
	# Input must be Cortex or Striatum.
	msg <- c("Please specify a tissue type to be analyzed:\n",
	 "Choose either 'Cortex' or 'Striatum'.")
	# If interactive, return default tissue.
	if (interactive()) { 
		return(default) 
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
  library(sva)
  library(grid)
  library(dplyr)
  library(WGCNA)
  library(edgeR)
  library(tibble)
  library(gtable)
  library(getPPIs)
  library(cowplot)
  library(ggplot2)
  library(gridExtra)
  library(flashClust)
  library(data.table)
})

# Load project specific data and functions.
suppressWarnings({ devtools::load_all() })

# Project directories.
datadir <- file.path(root,"data") # Input/Output data.
rdatdir <- file.path(root,"rdata") # Temporary data files.
tabsdir <- file.path(root,"tables") # Output excel tables.

# Number of threads for parallel processing.
nThreads <- parallel::detectCores() - 1

#---------------------------------------------------------------------
## Load the raw data and sample info (traits).
#---------------------------------------------------------------------
# The raw peptide intensity data were exported from ProteomeDiscover (PD)
# version 2.2. Note that the default export from PD2.x is a unitless signal to
# noise ratio, and it is not recommended to use ths for quantification.

# Load the TMT data.
data_files <- c(
	      "Cortex" = "Cortex_Raw_Peptides.csv",
	      "Striatum" = "Striatum_Raw_Peptides.csv"
	      )

# Load sample information.
meta_files <- c(
		"Cortex" = "Cortex_Samples.csv",
		"Striatum" = "Striatum_Samples.csv"
		)

# Load the TMT data and sample meta data.
raw_peptide <- fread(file = file.path(datadir, data_files[tissue]))
sample_info <- fread(file = file.path(datadir, meta_files[tissue]))

# Remove Iap immunoglobin protein.
drop <- c("P03975","")
raw_peptide <- raw_peptide %>% filter(Accession %notin% drop)

#---------------------------------------------------------------------
## Create gene identifier map.
#---------------------------------------------------------------------
message("\nCreating protein/gene identifier map.")

# Create gene map.
uniprot <- unique(raw_peptide$Accession)
entrez <- mgi_batch_query(ids=uniprot)

# Map missing ids by hand.
message("Mapping missing ids by hand.")
is_missing <- is.na(entrez)
missing <- c("P10853" = 319180) # H2bc7
entrez[names(missing)] <- missing
if (sum(is.na(entrez))>0) { stop("Missing identifiers.") }

# Map entrez to gene symbols.
symbols <- getPPIs::getIDs(entrez,from="entrez",to="symbol",species="mouse")
gene_map <- as.data.table(keep.rownames="uniprot",entrez)
gene_map$symbol <- symbols[as.character(gene_map$entrez)]

# Save gene_map as rda.
myfile <- file.path(datadir,paste0(tolower(tissue),"_gene_map.rda"))
save(gene_map,file=myfile,version=2)

#---------------------------------------------------------------------
## Sample loading normalization within experiments.
#---------------------------------------------------------------------
# The function __normalize_SL__ performs sample loading (SL) normalization to
# equalize the run level intensity (column) sums. The data in each column are
# multiplied by a factor such that the mean of the column sums are are equal.
# Sample loading normalization is performed within an experiment under the
# assumption that equal amounts of protein were used for each of the 11 TMT
# channels.
message("\nPerforming peptide-level sample loading normalization.")

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

# Filter peptides based on QC precision.
message("\nRemoving outlier QC peptides.")
filter_peptide <- filter_QC(SL_peptide, groups, nbins = 5, threshold = 4)

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
message("\nImputing missing peptide values with KNN algorithm.")

# Impute missing values using KNN algorithm for MNAR data.
# Rows with missing QC replicates are ingored (qc_threshold=0).
# Rows with more than 2 (50%) missing biological replicates are
# ignored (bio_threshold=2).
data_impute <- impute_peptides(filter_peptide, groups, method = "knn")
impute_peptide <- data_impute$data_imputed

#---------------------------------------------------------------------
## Protein level summarization and normalization across all batches.
#---------------------------------------------------------------------
# Summarize to protein level by summing peptide intensities. Note that the
# peptide column in the returned data frame reflects the total number of
# peptides identified for a given protein across all experiments.

# Summarize to protein level:
message("\nSummarizing proteins as the sum of their peptides.")
raw_protein <- summarize_protein(impute_peptide)

# Normalize across all columns.
message("\nPerforming protein-level sample loading normalization.")
SL_protein <- normalize_SL(raw_protein, colID="Abundance", group="Abundance")

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
message("\nUsing ComBat to remove intra-batch batch effect")

# Loop to perform ComBat on intraBatch batch effect (prep date).
# If there is no known batch effect (bicor(batch,PC1)<0.1) then
# the data is returned un-regressed.
# NOTE: The values of QC samples are not adjusted by ComBat.
# QC samples were prepared from a seperate batch of mice and
# represent a single batch.

# Loop:
data_in <- SL_protein
data_out <- list() # ComBat data.
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
  idx <- match(colnames(data), sample_info$Sample)
  traits_sub <- sample_info[idx, ]
  rownames(traits_sub) <- traits_sub$Sample
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
  if (!all(colnames(data) == CombatInfo$Sample)) {
    warning("Warning: Names of traits and expression data do not match.")
  }
  # Apply ComBat.
  message(paste("\nPerforming", groups[i], "ComBat..."))
  if (length(unique(CombatInfo$PrepDate)) > 1 & abs(r1) > 0.1) {
    # Create ComBat model.
    model <- model.matrix(~ as.factor(CombatInfo$Treatment),
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
  # Un-log.
  data_ComBat <- 2^data_ComBat
  # Recombine with QC data.
  data_out[[i]] <- cbind(info_cols[!rows_out, ], data_QC, data_ComBat)
  names(data_out)[[i]] <- group
} # Ends ComBat loop.

# Merge the data frames with purrr::reduce()
combat_protein <- data_out %>% 
	purrr::reduce(left_join, by = c(colnames(data_in)[c(1, 2)]))

#---------------------------------------------------------------------
## IRS Normalization.
#---------------------------------------------------------------------
# Internal reference sclaing (IRS) normalization equalizes the protein-wise means
# of reference (QC) samples across all batches. Thus, IRS normalization accounts
# for the random sampling of peptides at the MS2 level which results in the
# identification/quantificaiton of proteins by different peptides in each
# experiment. IRS normalization was first described by __Plubell et al., 2017__.
message("\nPerforming IRS normalization.")

# Perform IRS normalization.
IRS_protein <- normalize_IRS(combat_protein, "QC", groups, robust = TRUE)

#---------------------------------------------------------------------
## Identify and remove QC sample outliers.
#---------------------------------------------------------------------
# IRS normalization utilizes QC samples as reference samples. Outlier QC
# measurements (caused by interference or other artifact) would influence the
# create unwanted variability. Thus, outlier QC samples are removed, if
# identified. The method used by __Oldham et al., 2016__ is used to identify
# QC sample outliers. A threshold of -2.5 is used.
message("\nExamining QC samples for potential outliers using Oldham's method.")

# Illustrate Oldham's sample connectivity.
data_in <- IRS_protein
sample_connectivity <- ggplotSampleConnectivity(data_in,
  colID = "QC",
  threshold = oldham_threshold
)

# Check sample connectivity
tab <- sample_connectivity$table
df <- tibble::add_column(tab, "Name" = rownames(tab), .before = 1)
rownames(df) <- NULL
#knitr::kable(df)

# Loop to identify Sample outliers using Oldham's connectivity method.
n_iter <- 5
data_in <- IRS_protein
out_samples <- list()
plots <- list()
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
message(paste("\nTotal number of outlier QC samples identified:", 
	      sum(bad_samples != "none")))

# Remove outliers from data.
samples_out <- paste(bad_samples, collapse = "|")
out <- grepl(samples_out, colnames(SL_protein))

# Redo IRS after outlier removal.
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
message("\nRemoving irroproducibly quantified proteins.")

# Remove proteins that are identified by only 1 peptide as well as
# proteins identified in less than 50% of samples.
filter_protein <- filter_proteins(IRS_protein, "Abundance")

# Impute the remaining number of missing values with KNN.
message("\nImputing missing protein values.")
impute_protein <- impute_proteins(filter_protein, "Abundance", method = "knn")

#---------------------------------------------------------------------
## Collect the preprocessed data in a tidy format.
#---------------------------------------------------------------------

# Tidy the data.
tidy_prot <- reshape2::melt(impute_protein,id.var=c("Accession","Peptides"),
			    value.var="Intensity", variable.name="Sample",
			    value.name = "Intensity") %>% as.data.table()

# Merge data and meta data by sample name.
sample_info$Sample <- as.character(sample_info$Sample)
tidy_prot$Sample <- as.character(tidy_prot$Sample)
tidy_prot <- left_join(tidy_prot,sample_info,by="Sample")

# Drop QC.
tidy_prot <- tidy_prot %>% filter(Treatment != "QC")

#---------------------------------------------------------------------
## Create protein networks.
#---------------------------------------------------------------------

# Drop QC, log2 transform, and do bicor.
message("\nBuilding protein covariation network.")
dm <- tidy_prot %>% filter(!grepl("QC",Sample)) %>% 
	as.data.table() %>% 
	dcast(Sample ~ Accession, value.var="Intensity") %>% 
	as.matrix(rownames="Sample")
adjm <- WGCNA::bicor(log2(dm))

# Network enhancment of the bicor adjacency matrix.
message("\nPerforming network enhancment.")
ne_adjm <- neten::neten(adjm)

#--------------------------------------------------------------------
## Create PPI network.
#--------------------------------------------------------------------

# Load mouse PPIs.
message("\nBuilding PPI network.")
data(musInteractome) # from twesleyb/getPPIs.

# Collect all entrez cooresponding to proteins in our network.
proteins <- colnames(adjm)
entrez <- gene_map$entrez[match(proteins,gene_map$uniprot)]
names(proteins) <- entrez

# Collect PPIs among all proteins.
os_keep <- c(9606,10116,10090)
ppi_data <- musInteractome %>%
	filter(Interactor_B_Taxonomy %in% os_keep) %>%
	filter(Interactor_B_Taxonomy %in% os_keep) %>%
	filter(osEntrezA %in% entrez) %>% 
	filter(osEntrezB %in% entrez)

# Save to excel.
myfile <- file.path(root,"tables",paste0(tissue,"_Network_PPIs.xlsx"))
write_excel(list("Network PPIs" = ppi_data),file=myfile)

# Create simple edge list (sif) and matrix with node attributes (noa).
sif <- ppi_data %>% dplyr::select(osEntrezA, osEntrezB)
sif$uniprotA <- proteins[as.character(sif$osEntrezA)]
sif$uniprotB <- proteins[as.character(sif$osEntrezB)]

# Create igraph object from sif.
g <- graph_from_data_frame(sif[,c("uniprotA","uniprotB")], directed = FALSE)
g <- simplify(g)

# Extract as adjm.
ppi_adjm <- as.matrix(as_adjacency_matrix(g))

# Fill matrix.
all_proteins <- colnames(adjm)
missing <- all_proteins[all_proteins %notin% colnames(ppi_adjm)]
x <- matrix(nrow=dim(ppi_adjm)[1],ncol=length(missing))
colnames(x) <- missing
ppi_adjm <- cbind(ppi_adjm,x)
y <- matrix(nrow=length(missing),ncol=dim(ppi_adjm)[2])
rownames(y) <- missing
ppi_adjm <- rbind(ppi_adjm,y)
ppi_adjm <- ppi_adjm[all_proteins,all_proteins]
ppi_adjm[is.na(ppi_adjm)] <- 0

# Check, all should be the same.
c1 <- all(colnames(adjm) == colnames(ne_adjm)) 
c2 <- all(colnames(ne_adjm) == colnames(ppi_adjm))
if (!(c1 & c2)){ stop() }

# Number of edges and nodes.
n_edges <- sum(ppi_adjm[upper.tri(ppi_adjm,diag=FALSE)])
n_nodes <- ncol(ppi_adjm)
data.table("Nodes"=formatC(n_nodes,big.mark=",",format="d"),
	   "Edges"=formatC(n_edges,big.mark=",",format="d")) %>% knitr::kable()

# FIXME: examine topology of the ppi graph.

# Save as rda.
myfile <- file.path(datadir,paste0(tolower(tissue),"_ppi_adjm.rda"))
save(ppi_adjm,file=myfile,version=2)

#--------------------------------------------------------------------
## Intra-genotype contrats with EdgeR GLM.
#--------------------------------------------------------------------

message("\nAssessing intra-genotype differential abundance with EdgeR GLM.")

# Loop to perform intra-genotype WT v MUT comparisons:
results_list <- list()
for (geno in groups){
	# Cast data into matrix.
	dm <- tidy_prot %>% filter(Treatment != "QC") %>% 
		filter(Genotype == geno) %>% as.data.table() %>%
		dcast(Accession ~ Sample, value.var="Intensity") %>% 
		as.matrix(rownames=TRUE)
	# Create dge object.
	dge <- DGEList(counts=dm)
	# Perform TMM normalization.
	dge <- calcNormFactors(dge)
	# Sample to group mapping.
	samples <- rownames(dge$samples)
	idx <- match(samples,tidy_prot$Sample)
	genotype <- tidy_prot$Genotype[idx]
	treatment <- tidy_prot$Treatment[idx]
	#batch <- as.factor(tidy_prot$PrepDate[idx])
	dge$samples$group <- interaction(genotype,treatment)
	# Basic design matrix for GLM -- all groups treated seperately.
	design <- model.matrix(~ 0 + group, data = dge$samples)
	# Estimate dispersion:
	dge <- estimateDisp(dge, design, robust = TRUE)
	# Fit a general linear model.
	fit <- glmQLFit(dge, design, robust = TRUE)
	# Create contrast.
	wt <- colnames(design)[grepl("WT",colnames(design))]
	mut <- colnames(design)[grepl("HET|KO",colnames(design))]
	contr <- limma::makeContrasts(paste(mut,wt,sep="-"),levels=design)
	# Assess differences.
	qlf <- glmQLFTest(fit,contrast=contr)
	#qlf <- glmQLFTest(fit)
	# Call topTags to add FDR. Gather tabularized results.
	glm_results <- topTags(qlf, n = Inf, sort.by = "PValue")$table
	# Use gene map to annotate glm_results with entrez Ids and gene symbols.
	idx <- match(rownames(glm_results), gene_map$uniprot)
	glm_results <- add_column(glm_results, "Accession"=rownames(glm_results),
				  .before = 1)
	glm_results <- add_column(glm_results, "Entrez" = gene_map$entrez[idx], 
				  .after = 1)
	glm_results <- add_column(glm_results, "Symbol" = gene_map$symbol[idx], 
				  .after = 2)
	# Add expression data.
	wt_cols <- grep("WT",colnames(dm))
	mut_cols <- grep("HET|KO",colnames(dm))
	dt <- as.data.table(dm[,c(wt_cols,mut_cols)],keep.rownames="Accession")
	glm_results <- left_join(glm_results,dt,by="Accession")
	results_list[[geno]] <- glm_results
}

# Examine results.
message("\nSummary of DA proteins for intra-genotype contrasts:")
df <- sapply(results_list,function(x) sum(x$FDR<FDR_threshold))
knitr::kable(t(df))

intrabatch_results <- results_list

#---------------------------------------------------------------------
## Perform TAMPOR normalization.
#---------------------------------------------------------------------
message("\nPerforming TAMPOR normalization.")

# batch.channel annotations
sample_info$Batch <- paste0("b",as.numeric(as.factor(sample_info$Genotype)))
sample_info$Batch.Channel <- as.character(interaction(sample_info$Batch,sample_info$Channel))

# The data, with batch.channel names.
dm <- tidy_prot %>% filter(Treatment != "QC") %>% as.data.table() %>%
	dcast(Accession ~ Sample, value.var="Intensity") %>% 
	as.matrix(rownames=TRUE)
idx <- match(colnames(dm),sample_info$Sample)
colnames(dm) <- sample_info$Batch.Channel[idx]

# traits to be passed to TAMPOR.
# remove any samples that were dropped from traits.
traits <- sample_info[sample_info$Batch.Channel %in% colnames(dm),]
traits <- traits %>% dplyr::select(Sample,Batch,Channel,Batch.Channel,
				   Genotype,Treatment)
rownames(traits) <- traits$Batch.Channel

# GIS index is all WT samples.
controls <- traits$Batch.Channel[grepl("WT", traits$Treatment)]

# Perform TAMPOR.
results <- TAMPOR(
  dat = dm,
  traits = traits,
  batchPrefixInSampleNames = TRUE,
  samplesToIgnore = "None",
  GISchannels = controls,
  parallelThreads = nThreads
)

# Collect normalized abundance data.
TAMPOR_dm <- results$cleanRelAbun
idy <- match(colnames(TAMPOR_dm),traits$Batch.Channel)
colnames(TAMPOR_dm) <- traits$Sample[idy] # map back to sample name

# Merge with tidy_prot.
norm_prot <- TAMPOR_dm %>% as.data.table(keep.rownames="Accession") %>%
	reshape2::melt(id.var="Accession")
colnames(norm_prot) <- c("Accession", "Sample","TAMPOR.Intensity")
norm_prot$Sample <- as.character(norm_prot$Sample)
tidy_prot <- left_join(tidy_prot %>% filter(Treatment != "QC"),
		       norm_prot,by=c("Accession","Sample"))

#---------------------------------------------------------------------
## EdgeR statistical comparisons post-TAMPOR.
#---------------------------------------------------------------------
# Statistical comparisons are KO/HET versus all WT of a tissue type.
message("\nAssessing differential abundance using combined WT controls with EdgeR GLM.")

# Prepare data for EdgeR.
# Data should NOT be log2 transformed.
dm <- tidy_prot %>% as.data.table() %>%
	dcast(Accession ~ Sample, value.var="TAMPOR.Intensity") %>%
	as.matrix(rownames="Accession")

# Create dge object.
dge <- DGEList(counts=dm)

# Perform TMM normalization.
dge <- calcNormFactors(dge)

# Sample to group mapping.
samples <- rownames(dge$samples)
idx <- match(samples,tidy_prot$Sample)
genotype <- tidy_prot$Genotype[idx]
treatment <- tidy_prot$Treatment[idx]
group <- as.character(interaction(genotype,treatment))
group[grep("WT",group)] <- "WT"
dge$samples$group <- group

# Basic design matrix for GLM -- all groups treated seperately.
design <- model.matrix(~ 0 + group, data = dge$samples)

# Estimate dispersion:
dge <- estimateDisp(dge, design, robust = TRUE)

# Fit a general linear model.
fit <- glmQLFit(dge, design, robust = TRUE)

# Create contrast.
contr_list <- list(
		   "Shank2" = limma::makeContrasts('groupWT-groupShank2.KO',
						   levels=design),
		   "Shank3" = limma::makeContrasts('groupWT-groupShank3.KO',
						   levels=design),
		   "Syngap1" = limma::makeContrasts('groupWT-groupSyngap1.HET',
						    levels=design),
		   "Ube3a" = limma::makeContrasts('groupWT-groupUbe3a.KO',
						  levels=design)
		   )

# Assess differences.
qlf <- lapply(contr_list,function(x) glmQLFTest(fit,contrast=x))

# Call topTags to add FDR. Gather tabularized results.
glm_results <- lapply(qlf, function(x) {
			      topTags(x, n = Inf, sort.by = "PValue")$table})

# Use gene map to annotate glm_results with entrez Ids and gene symbols.
glm_results <- lapply(glm_results,function(x) {
       idx <- match(rownames(x), gene_map$uniprot)
       x <- add_column(x, "Accession"=rownames(x),.before = 1)
       x <- add_column(x, "Entrez" = gene_map$entrez[idx], .after = 1)
       x <- add_column(x, "Symbol" = gene_map$symbol[idx], .after = 2)
       return(x)
})
names(glm_results) <- groups

# Combine stats into single df and merge with tidy_prot.
glm_stats <- bind_rows(glm_results,.id="Genotype")
tidy_prot <- left_join(tidy_prot,glm_stats,by=c("Accession","Genotype"))

# Add expression data to glm_results.
glm_results <- lapply(names(glm_results),function(x) {
	wt_cols <- which(grepl("WT",colnames(dm)) & grepl(x,colnames(dm)))
	mut_cols <- which(grepl("HET|KO",colnames(dm)) & grepl(x,colnames(dm)))
	dt <- as.data.table(dm[,c(wt_cols,mut_cols)],keep.rownames="Accession")
	glm_results[[x]] <- left_join(glm_results[[x]],dt,by="Accession")
})
names(glm_results) <- groups

# Summary of DA proteins.
message("\nSummary of DA proteins:")
df <- sapply(glm_results,function(x) sum(x$FDR<FDR_threshold))
knitr::kable(t(df))

# Save results to file as excel spreadsheet.
myfile <- file.path(tabsdir,paste0(tissue,"_Protein_GLM_Results.xlsx"))
names(intrabatch_results) <- paste(names(intrabatch_results),"intrabatch")
names(glm_results) <- paste(names(glm_results),"Results")
results_list <- c("Raw peptide" = list(raw_peptide),
		  "Normalized protein" = list(as.data.table(dm,keep.rownames="Accession")),
		  intrabatch_results,
		  glm_results)
write_excel(results_list,myfile)

#---------------------------------------------------------------------
## Save output data.
#---------------------------------------------------------------------
message("\nSaving data for downstream analysis.")

## Output in root/rdata

# ne_adjm.csv - adjacency matrix saved as csv for Leidenalg.
myfile <- file.path(rdatdir,paste0(tolower(tissue),"_ne_adjm.csv"))
ne_adjm %>% as.data.table(keep.rownames="Accession") %>% fwrite(myfile)

## Output in root/data

# Save adjm as rda.
myfile <-file.path(datadir,paste0(tolower(tissue),"_adjm.rda"))
save(adjm,file=myfile,version=2)

# ne_adjm.rda
myfile <-file.path(datadir, paste0(tolower(tissue),"_ne_adjm.rda"))
save(ne_adjm,file=myfile,version=2)

# tidy_prot.rda
myfile <- file.path(datadir,paste0(tolower(tissue),".rda"))
tidy_prot <- tidy_prot %>% filter(!grepl("QC",Sample)) %>% as.data.table() # Drop QC!
save(tidy_prot,file=myfile,version=2)

message("\nDone!")
