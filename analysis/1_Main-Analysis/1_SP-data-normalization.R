#!/usr/bin/env Rscript

# title: Analysis of Synaptopathy TMT Proteomics
# description: preprocessing and statistical analysis of TMT proteomics
# author: Tyler W Bradshaw

## ---- INPUTs
# project root
root <- "~/projects/SynaptopathyProteomics"

# the TMT data in root/rdata:
# peptide level data exported from PD
data_files <- c(
  "Cortex" = "Cortex_Raw_Peptides.csv",
  "Striatum" = "Striatum_Raw_Peptides.csv"
)

# the sample metadata in root/rdata:
meta_files <- c(
  "Cortex" = "Cortex_Samples.csv",
  "Striatum" = "Striatum_Samples.csv"
)

# NOTE: all 4 TMT mixtures were processed together in proteome discover and saved as
# single csv file. Data is processed togther, and normalization is performed for each 
# mixture independently w/ the exception of IRS normalization which utilizes the
# QC sample analyzed in technical replication in all mixtures to normalize
# protein measurements between experiments. This step adjusts Biological groups
# equivalently and has no effect on intra-Mixture comparisons. 


## ---- R ENVIRONMENT 

# Prepare the R workspace for the analysis. Load custom functions and prepare
# the project directory for saving output files.

# load required packages
suppressPackageStartupMessages({
  library(dplyr) # for working with data
  library(tibble) # for add_column function
  library(data.table) # for working with tables
})

# load project specific data and functions
devtools::load_all(root)

# parse input
tissue <- parseArgs()

# gene mapping data.frame
data(gene_map)

# define project's directories:
datadir <- file.path(root, "data") # input/output data
rdatdir <- file.path(root, "rdata") # temporary data files
tabsdir <- file.path(root, "tables") # output excel tables


## ---- MAIN 
# Load the raw data and sample meta data (aka sample_info, traits)
# The raw peptide intensity data were exported from ProteomeDiscover (PD)
# version 2.2. Note that the default export from PD2.x is a unitless signal to
# noise ratio, and it is not recommendetraits = sample_traits/d to use ths for quantification.

# load the TMT data and sample metadata
raw_peptide <- fread(file = file.path(root,"inst","extdata", data_files[tissue]))
sample_info <- fread(file = file.path(root,"inst","extdata", meta_files[tissue]))


## ---- Summarize initial number of peptides and protein 
# The initial number of unique proteins and unique peptides quantified
# in all four experiments.

message("\nInitial number of peptides and proteins:")

nsamples <- sum(grepl("Abundance", colnames(raw_peptide)))
npep <- nrow(raw_peptide)
nprot <- raw_peptide %>%
  dplyr::select(Accession) %>%
  unlist() %>%
  unique() %>%
  length()
df <- data.table(
  "N Peptides" = formatC(npep, big.mark = ","),
  "N Proteins" = formatC(nprot, big.mark = ","),
  "N Samples" = nsamples
)
knitr::kable(df)


## ---- Sample loading normalization within experiments 
# The function normalize_SL performs sample loading (SL) normalization to
# equalize the sample-level intensity (column) sums. The data from each sample
# are multiplied by a factor such that the mean of the sample sums are are
# equal.  Sample loading normalization is performed within an experiment under
# the assumption that equal amounts of protein were used for each of the 11 TMT
# channels.

message("\nPerforming peptide-level sample loading normalization.")

# define data columns for SL within experiments:
colID <- "Abundance" # Column identifier of numeric data
groups <- c("Shank2", "Shank3", "Syngap1", "Ube3a")

# perform initial SL normalization:
sl_peptide <- normalize_SL(raw_peptide, colID, groups)


## ---- Peptide level filtering based on QC samples 
# Peptides that were not quantified in all three QC replicates are removed.
# The data are binned by intensity, and measurments that are more than 4
# standard deviations away from the mean ratio of the intensity bin are
# considered outliers and removed.

# filter peptides based on QC precision
message("\nRemoving outlier QC peptides.")
filt_peptide <- filter_QC(sl_peptide, groups, nbins = 5, threshold = 4)


## ---- Impute missing peptide values within an experiment 

# The function impute_peptides performs KNN imputing for values that are missing
# at random data. Impution is performed with an experiment, and rows with more
# than 50% missing values are censored and will not be imputed. Peptides with
# more than 2 missing biological replicates or any missing quality control (QC)
# replicates will be censored and are not imputed.  NOTE: there shouldn't really
# be any missing values within each TMT 11-plex experiment due to the
# simultanous quantification of species across all 11 channels. Thus, missing
# values may trully be missing at random.

# NOTE: Rows with missing QC replicates are ingored (qc_threshold=0).
# NOTE: Rows with more than 2 (50%) missing biological replicates are ignored
# (bio_threshold=2).
message("\nImputing missing peptide values with KNN(k=10) algorithm.")
data_impute <- impute_peptides(filt_peptide, groups, method = "knn")
impute_peptide <- data_impute$data_imputed


## ---- Protein level summarization
# Summarize to protein level by summing peptide intensities. Note that the
# peptide column in the returned data frame reflects the total number of
# peptides identified for a given protein across all experiments.

# summarize to protein level:
message("\nSummarizing proteins as the sum of their peptides.")
raw_protein <- summarize_protein(impute_peptide)


## ---- tidy raw protein

idy <- colnames(raw_protein)[grepl("Abundance",colnames(raw_protein))]
idx <- !apply(raw_protein[,idy],1,function(x) any(is.na(x)))

tidy_prot <- raw_protein %>% 
	filter(idx) %>%
	reshape2::melt(id.var=c("Accession","Peptides"),
					    value.var = idy,
					    value.name = "Intensity",
					    variable.name="Sample") %>% 
                             left_join(sample_info, by = "Sample") 

#myfile <- file.path(root,"rdata","tidy_prot.rda")
#save(tidy_prot,file=myfile,version=2)


## ---- IRS Normalization across all batches 
# Internal reference sclaing (IRS) normalization equalizes the protein-wise means
# of reference (QC) samples across all batches. Thus, IRS normalization accounts
# for the random sampling of peptides at the MS2 level which results in the
# identification/quantificaiton of proteins by different peptides in each
# experiment. IRS normalization was first described by Plubell et al., 2017.
message("\nPerforming IRS normalization.")

irs_protein <- normalize_IRS(raw_protein, "QC", groups, robust = TRUE)


## ---- Protein level filtering

# Proteins that are identified by only a single peptide are removed. 
message("\nRemoving irroproducibly quantified proteins.")

filt_protein <- filter_proteins(irs_protein, "Abundance")


## ---- Collect the preprocessed data in a tidy format 

# tidy the data
tidy_prot <- reshape2::melt(filt_protein,
  id.var = c("Accession", "Peptides"),
  value.var = "Intensity", variable.name = "Sample",
  value.name = "Intensity"
) %>% as.data.table()

# merge TMT data and sample metadata by sample name
sample_info$Sample <- as.character(sample_info$Sample)
tidy_prot$Sample <- as.character(tidy_prot$Sample)
tidy_prot <- left_join(tidy_prot, sample_info, by = "Sample") %>%
  #filter(Treatment != "SPQC") %>%
  mutate(Condition = interaction(Genotype,Treatment)) %>%
  mutate(Batch = as.numeric(as.factor(PrepDate))) %>%
  as.data.table()


## complete cases!
colnames(tidy_prot)
dm <- tidy_prot %>% 
	reshape2::dcast(Accession ~ Sample, value.var="Intensity") %>% 
	as.data.table() %>% as.matrix(rownames="Accession")
idx <- apply(dm, 1, function(x) any(is.na(x)))
drop <- rownames(dm)[idx]

tidy_prot <- tidy_prot %>% filter(!(Accession %in% drop))

## final munge
tidy_prot$Condition = as.character(tidy_prot$Condition)
tidy_prot$Genotype = as.character(tidy_prot$Genotype)

# Condition == WT|HET|SPQC
idx <- grepl("QC",tidy_prot$Condition)
tidy_prot$Condition[idx] <- "SPQC"
tidy_prot$Condition[!idx] <- sapply(strsplit(tidy_prot$Condition[!idx],"\\."),"[",2)

# Batch == 0-5 (which synaptosome batch prep)
idx <- is.na(tidy_prot$Batch)
tidy_prot$Batch[idx] <- 0


## ---- Summarize final number of peptides and proteins 

tidy_prot %>% summarize(nSamples=length(unique(Sample)),
			nProteins=length(unique(Accession))) %>% knitr::kable()


## ---- Save output data 

# tidy_prot.rda
myfile <- file.path(datadir, paste0(tolower(tissue), ".rda"))
save(tidy_prot, file = myfile, version = 2)
message("saved: ", myfile)
