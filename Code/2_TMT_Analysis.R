#' ---
#' title: TMT Analysis part 2. TAMPOR normalization. 
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

#-------------------------------------------------------------------------------
#' ## Prepare the workspace.
#-------------------------------------------------------------------------------
#+ eval = TRUE, echo = FALSE, error = FALSE

# Use ctl+alt+T to execute a code chunk.

# Run this chunk before doing anything!
rm(list = ls())
dev.off()
cat("\014") # alternative is cat("\f")
options(stringsAsFactors = FALSE)

# Sometimes, if you have not cleared the workspace of all loaded packages,
# you man incounter problems.
# To remove all packages, you can call the following:
library(magrittr)
library(JGmisc)
detachAllPackages(keep = NULL)

#  Load required packages.
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

# To install TBmiscr:
#library(devtools)
#devtools::install_github("twesleyb/TBmiscr")

# Define version of the code.
CodeVersion <- "TMT_Analysis_part2"

# Define tisue type: cortex = 1; striatum = 2; 3 = combined. 
type <- 3
tissue <- c("Cortex", "Striatum", "Combined")[type]

# Set the working directory.
rootdir <- "D:/Documents/R/Synaptopathy-Proteomics"
setwd(rootdir)

# Set any other directories.
functiondir <- paste(rootdir, "Functions", sep = "/")
datadir <- paste(rootdir, "Input", sep = "/")
Rdatadir <- paste(rootdir,"RData", sep = "/")

# Create code-version specific figure and tables folders if they do not already exist.
# Creat otuput direcotry for figures.
outputfigs <- paste(rootdir, "Figures", tissue, sep = "/")
outputfigsdir <- paste(outputfigs, CodeVersion, sep = "/")
if (!file.exists(outputfigsdir)) {
  dir.create(file.path(outputfigsdir))
} else {
  print("This directory already exists. Warning: Some files may be overwritten when running this script.")
}
# Create output directory for tables.
outputtabs <- paste(rootdir, "Tables", tissue, sep = "/")
outputtabsdir <- paste(outputtabs, CodeVersion, sep = "/")
if (!file.exists(outputtabsdir)) {
  dir.create(file.path(outputtabsdir))
} else {
  print("This directory already exists. Warning: Some files may be overwritten when running this script.")
}
# Create output directory for reports.
outputreports <- paste(rootdir, "Reports", tissue, sep = "/")
outputrepsdir <- paste(outputreports, CodeVersion, sep = "/")
if (!file.exists(outputrepsdir)) {
  dir.create(file.path(outputrepsdir))
} else {
  print("This directory already exists. Warning: Some files may be overwritten when running this script.")
}

# Load required custom functions.
functiondir <- paste(rootdir, "Functions", sep = "/")
my_functions <- paste(functiondir, "TMT_Preprocess_Functions.R", sep = "/")
source(my_functions)

# Define prefix for output figures and tables.
outputMatName <- paste(tissue, "_TMT_Analysis", sep = "")

# Globally set ggplots theme.
ggplot2::theme_set(theme_gray())

#-------------------------------------------------------------------------------
#' ## Merge cortex and striatum data.
#-------------------------------------------------------------------------------

################################################################################
## Note: Skip ahead to next chunk if you wish to load the data from file.     ##
################################################################################

## Merge Traits. 

# Load the cortex and striatum traits files.
inputTraitsCSV <- c(
  "4227_TMT_Cortex_Combined_PD_Protein_Intensity_EBD_traits.csv",
  "4227_TMT_Striatum_Combined_PD_Protein_Intensity_EBD_traits.csv"
)

# Load the sample info into a list, traits.
files <- paste(datadir,inputTraitsCSV, sep = "/")
traits <- lapply(files, function(x) read.csv(file = x, header = TRUE))

# Bind the elements of the list. 
traits <- do.call(rbind,traits)

# Rename SampleIDs
idx <- c(45:nrow(traits))
vars <- strsplit(traits$SampleID,"\\.")[idx]
new_batch <- c(rep("b5",11),rep("b6",11),rep("b7",11),rep("b8",11))
channel <- sapply(vars, "[", 2)
ids <- paste(new_batch,channel,sep=".")
traits$SampleID[idx] <- ids
traits[1:5,1:5]

# Insure that rownames are new sampleIDS
rownames(traits) <- traits$SampleID

# Add column for tissue type.
traits$Tissue <- c(rep("Cortex",44),rep("Striatum",44))
dim(traits)

## Merge expression data. 

# Load the Cortex and Striatum IRS + eBLM regressed data. 
type <- c("Cortex","Striatum")
files <- paste(Rdatadir,type,"CleanDat_IRS_eBLM_TAMPOR_format.Rds",sep="/")
data <- list(
  cortex_protein <- readRDS(files[1]),
  striatum_protein <- readRDS(files[2])
)

# Fortify and add accession column
data_fort <- lapply(data,
                    function(x) add_column(as.data.frame(x),ID=rownames(x),.before=1))

# Bind by Accession
# Warning that ID has different attributes is okay... I think.
data_merge <- data_fort %>% reduce(left_join, by = "ID")

# Remove rows with missing values.
data_clean <- na.omit(data_merge)

# SL across all columns.
data_norm <- normalize_SL(data_clean,"b","b")

## Clean-up formatting for TAMPOR.
rownames(data_norm) <- data_norm$info_cols
data_norm$info_cols <- NULL

# Batch b1-b4 are cortex. Batches b5-b8 are straitum. 
col_names <- colnames(data_norm)
# Cortex = group1...
group1 <- col_names[grepl(".x",col_names)]
group1 <- gsub(".x","",group1)
# Striatum = group2...
group2 <- col_names[!grepl(".x",col_names)]
group2 <- gsub(".y","",group2)
group2 <- gsub("b1","b5",group2)
group2 <- gsub("b2","b6",group2)
group2 <- gsub("b3","b7",group2)
group2 <- gsub("b4","b8",group2)

# Change column names to batch.channel. 
colnames(data_norm) <- c(group1,group2)

# Write as cleanDat
cleanDat <- data_norm
cleanDat[1:5,1:5]
dim(cleanDat)  # one sample has been removed. 

# GIS index is all WT samples.
controls <- colsplit(traits$SampleID[grepl("WT",traits$SampleType)],"\\.",c("batch","channel"))
controls

# Save merged data and traits to file.
datafile <- paste(Rdatadir,"Combined_Cortex_Striatum_cleanDat.Rds",sep="/")
saveRDS(cleanDat,datafile)
traitsfile <- paste(Rdatadir,"Combined_Cortex_Striatum_traits.Rds",sep="/")
saveRDS(traits,traitsfile)

#-------------------------------------------------------------------------------
#' ## User Parameters to change for TAMPOR Normalization strategy.
#-------------------------------------------------------------------------------

# Load traits and data from file.
#datafile <- paste(Rdatadir,tissue,"CleanDat_IRS_TMM_TAMPOR_format.Rds",sep="/")
#traitsfile <- paste(Rdatadir,tissue,"Sample_info.Rds",sep="/")
datafile <- paste(Rdatadir,"Combined_Cortex_Striatum_cleanDat.Rds",sep="/")
traitsfile <- paste(Rdatadir,"Combined_Cortex_Striatum_traits.Rds",sep="/")

cleanDat <- readRDS(datafile)
traits <- readRDS(traitsfile)
rownames(traits) <- traits$SampleID

cleanDat.original <- cleanDat
cleanDat[1:5,1:5]
dim(cleanDat)
dim(traits)

# Define batch and sampleIndex:
sampleIndex <- as.data.frame(do.call(rbind, strsplit(colnames(cleanDat), "\\.")))
colnames(sampleIndex) <- c("batch", "channel")
batchIndex <- unique(sampleIndex$batch)

# Parameters for normalization:
# if TRUE, all randomized samples in the batch are used.
# if FALSE, only GIS channels are used
useALLnonGIS <- FALSE 

# Define columns for used as GIS for TAMPOR normalization:
# GIS index is all WT samples.
controls <- colsplit(traits$SampleID[grepl("WT",traits$SampleType)],"\\.",c("batch","channel"))
dim(controls)                     
GISchannels <- controls$channel # WT control samples. 
iterations <- 250
samplesToIgnore <- c("NONE")
#samplesToIgnore <- c('b1.126','b1.129N','b1.131C',
#                     'b2.126','b2.129N','b2.131C',
#                     'b3.126','b3.129N','b3.131C',
#                     'b4.126','b4.129N','b4.131C') # All QC samples. 
inputMatName <- "TMT_Analysis" # prefix for output, should be similar to prefix of inputTraitsCSV.
meanORmedian <- "median" # must be a valid R function, e.g. mean or median --this text to call the function is substituted in eval(parse(text=... wrapped code at these positions ***
removeGISafter <- FALSE # Should samples designated as GIS be removed before network-related code runs?

# Set threshold for convergence:
FNthreshold <- 1e-8

#-------------------------------------------------------------------------------
#' ## Initiate loop for normalization.
#-------------------------------------------------------------------------------

################################################################################
## Note: Skip ahead to "Start WGCNA" if you have already run TAMPOR           ##
################################################################################

# Ignore samples in samplesToIgnore vector
if (length(na.omit(match(samplesToIgnore, colnames(cleanDat)))) == length(samplesToIgnore)) {
  for (i in 1:length(samplesToIgnore)) {
    cleanDat[, samplesToIgnore[i]] <- as.vector(rep(NA, nrow(cleanDat)))
  }
} else {
  cat("NOTE: one or more samplesToIgnore do not match sample names (colnames) in input abundance data.\nNot ignoring any samples.\n")
}

GISchannels <- as.character(GISchannels) #**** make sure channels is a vector of strings
sampleIndex <- as.data.frame(do.call(rbind, strsplit(colnames(cleanDat), "\\.")))[, 1:2]
colnames(sampleIndex) <- c("batch", "channel")
traits <- traits[match(colnames(cleanDat), rownames(traits)), ] # insure trait rows match abundance columns
# sampleIndex$channel <- as.character(traits$Group) #**** added as.character, since GISchannels is specified as a character vector
batchIndex <- unique(sampleIndex$batch)

## NOTE: Normalization 'channels' are specified for mean or median initial denominator in below code chunk:
## *Now allows setting GISindices even if your GIS is on different channels in different batches, essentially if channel order is not the same in each batch.
## This enables using TAMPOR for batched LFQ or LFQ of different cohorts.
GISindices <- list() # *GISindices as a list added THIS VERSION 2/5/2019
i.prev <- 0
offset <- 0
iter <- 0
for (i in batchIndex) {
  GISindices[[i]] <- vector()
  iter <- iter + 1
  for (denomChannel in GISchannels) {
    GISindices[[i]] <- c(GISindices[[i]], which(sampleIndex$batch == unique(sampleIndex$batch)[iter] & sampleIndex$channel == denomChannel))
  }
  if (!i.prev == 0) {
    offset <- offset + length(which(sampleIndex$batch == i.prev))
  }
  GISindices[[i]] <- unique(GISindices[[i]]) - offset # prior code assuming 'channel' names occur once per a batch: [match(unique(GISchannels), GISchannels)]
  i.prev <- i
}

iterationTrackingDF <- data.frame(Iteration = 1:iterations, FrobeniusNorm = NA, FrobenPrev = NA, FrobenDiff = NA, FrobenOverFirstFroben = NA)
cat(paste0("Starting # of Rows in data: ", nrow(cleanDat), ".\n"))

# START TAMPOR:
for (repeats in 1:iterations) {
  #-------------------------------------------------------------------------------
  # STEP 1a. Ratio data and prepare to row-normalize
  ratioedBatches <- batchGISavgs <- list()
  ratioCleanDatUnnorm <- data.frame(row.names = rownames(cleanDat))
  withinBatchGISgeomeans <- withinBatchRowGeomeans <- data.frame(row.names = rownames(cleanDat))
  
  for (batch in batchIndex) {
    tempForAvg <- matrix()
    tempForAvg <- as.data.frame(as.matrix(cleanDat[, which(sampleIndex$batch == batch)][, GISindices[[batch]] ], nrow = nrow(cleanDat), ncol = dim(cleanDat[, which(sampleIndex$batch == batch)][, GISindices[[batch]] ])[2])) # GISindices as list added THIS VERSION 2/5/2019
    batchGISavgs[[batch]] <- apply(tempForAvg, 1, function(x) eval(parse(text = paste0(meanORmedian, "(x,na.rm=TRUE)")))) # ADDED na.rm v04 ##MEAN/MEDIAN FUNCTION CHOICE***
    ratioedBatches[[batch]] <- cleanDat[, which(sampleIndex$batch == batch)] / batchGISavgs[[batch]]
    
    ## Below unnormed ratio data are only assembled for graphing purposes, for comparison to step 1b and final step 2 output
    ratioCleanDatUnnorm <- cbind(ratioCleanDatUnnorm, ratioedBatches[[batch]])
    ## If batches are randomized channels distributing cases and controls evenly across all batches, useALLnonGIS==TRUE
    if (useALLnonGIS) {
      withinBatchRowGeomeans <- withinBatchGISgeomeans <- as.data.frame(cbind(withinBatchRowGeomeans, as.data.frame(as.matrix(apply(ratioedBatches[[batch]][, -GISindices[[batch]] ], 1, function(x) eval(parse(text = paste0("2^", meanORmedian, "(log2(na.omit(x)))")))), ncol = dim(ratioedBatches[[batch]])[2], nrow = dim(ratioedBatches[[batch]])[1])))) ## MEAN/MEDIAN FUNCTION CHOICE***
      # as.matrix(), NOT matrix()
    } else {
      ## If we cannot rely on the robust assumption of batch-to-batch biological equivalence (with randomized sample order across all avalable channels in all batches), then use robust mean of GIS samples only
      withinBatchGISgeomeans <- withinBatchRowGeomeans <- as.data.frame(cbind(withinBatchRowGeomeans, as.data.frame(as.matrix(apply(ratioedBatches[[batch]], 1, function(x) eval(parse(text = paste0("2^", meanORmedian, "(log2(na.omit(x[GISindices[[batch]] ])))")))), ncol = dim(ratioedBatches[[batch]])[2], nrow = dim(ratioedBatches[[batch]])[1]))))
      # as.matrix(), NOT matrix()
    } ## MEAN/MEDIAN FUNCTION CHOICE***
  }
  
  #-------------------------------------------------------------------------------
  ## Step 1b. Complete row-normalization. This step normalizes rows within batch by batchCorrFactors; ratioCleanDatUnnorm does not go through this step
  meanBatchGeomeans <- apply(withinBatchRowGeomeans, 1, function(x) eval(parse(text = paste0(meanORmedian, "(na.omit(x))")))) ## MEAN/MEDIAN FUNCTION CHOICE***
  
  ## Rowwise (RW) relative abundances from GIS (or representative samples),
  # if we want to take the whole protein row (across batches) back to abundance after normalization step2 is complete**
  RW.relAbunFactors <- RW.GISavgs <- data.frame(row.names = rownames(cleanDat))
  RW.GISavgs <- cbind(apply(data.frame(column = batchIndex), 1, function(x) batchGISavgs[[x]]))
  colnames(RW.GISavgs) <- batchIndex
  RW.relAbunFactors <- apply(RW.GISavgs, 1, function(x) eval(parse(text = paste0("2^", meanORmedian, "(log2(na.omit(x)))")))) #** relative abundance multipliers for recovery of relative abundance (all rows from input cleanDat) ##MEAN/MEDIAN FUNCTION CHOICE***
  rownames(RW.GISavgs) <- names(RW.relAbunFactors) <- rownames(cleanDat)
  
  ## Calculate Step 1b multipliers to complete RW normalization
  batchCorrFactors <- meanBatchGeomeans / withinBatchRowGeomeans ## Tried FLIPPING, but modules changed somewhat.. ratio for LV analysis using PBS+Monomer for step 1a ratio***
  colnames(batchCorrFactors) <- colnames(withinBatchRowGeomeans) <- batchIndex
  
  normedBatches <- list()
  ratioCleanDatNorm <- as.matrix(data.frame(row.names = rownames(cleanDat)))
  for (batch in batchIndex) {
    normedBatches[[batch]] <- ratioedBatches[[batch]] * batchCorrFactors[, batch] # rowwise correct by multiplers (batchCorrFactors)
    ratioCleanDatNorm <- cbind(ratioCleanDatNorm, normedBatches[[batch]])
  }
  #-------------------------------------------------------------------------------
  
  ## log2 transform ratios output from step 1a, and 1b
  cleanDat.log2.ratioUnnorm <- log2(ratioCleanDatUnnorm)
  cleanDatNormNoColScaling <- log2(ratioCleanDatNorm)
  
  
  ## For cleanDat.log2.ratioUnnorm:
  ## Enforce <50% missingness (1 less than half of columns (or round down half if odd number of columns))
  LThalfSamples <- length(colnames(cleanDat.log2.ratioUnnorm)) / 2
  LThalfSamples <- LThalfSamples - if ((length(colnames(cleanDat.log2.ratioUnnorm)) %% 2) == 1) {
    0.5
  } else {
    1.0
  }
  
  removedRownames1 <- rownames(cleanDat.log2.ratioUnnorm[which(rowSums(as.matrix(is.na(cleanDat.log2.ratioUnnorm))) > LThalfSamples), ]) # list rows to be removed
  removedRownames1
  
  # remove rows with >=50% missing values (only if there are some rows to be removed)
  if (length(na.omit(match(removedRownames1, rownames(cleanDat.log2.ratioUnnorm)))) == length(removedRownames1) & length(removedRownames1) > 0) {
    cleanDat.log2.ratioUnnorm <- cleanDat.log2.ratioUnnorm[-match(removedRownames1, rownames(cleanDat.log2.ratioUnnorm)), ]
  } else {
    cat("")
  } # "no rows removed.\n"); }
  # can be error prone:  if ( !is.null(rownames(cleanDat.log2.ratioUnnorm[which(rowSums(as.matrix(is.na(cleanDat.log2.ratioUnnorm)))>LThalfSamples),])) ) { if( !length(rownames(cleanDat.log2.ratioUnnorm[which(rowSums(as.matrix(is.na(cleanDat.log2.ratioUnnorm)))>LThalfSamples),]))==0) { IndexHighMissing=which(rowSums(as.matrix(is.na(cleanDat.log2.ratioUnnorm)))>LThalfSamples); cleanDat.log2.ratioUnnorm<-cleanDat.log2.ratioUnnorm[-IndexHighMissing,]; } }
  dim(cleanDat.log2.ratioUnnorm) # 3166 rows.
  
  
  
  ## For cleanDatNormNoColScaling and companion relative abundance factors:
  ## Enforce <50% missingness (1 less than half of cleanDatNormNoColScaling columns (or round down half if odd number of columns))
  LThalfSamples <- length(colnames(cleanDatNormNoColScaling)) / 2
  LThalfSamples <- LThalfSamples - if ((length(colnames(cleanDatNormNoColScaling)) %% 2) == 1) {
    0.5
  } else {
    1.0
  }
  
  removedRownames <- rownames(cleanDatNormNoColScaling[which(rowSums(as.matrix(is.na(cleanDatNormNoColScaling))) > LThalfSamples), ]) # list rows to be removed
  #  removedRownames
  
  # remove rows with >=50% missing values (only if there are some rows to be removed)
  if (length(as.vector(na.omit(match(removedRownames, rownames(cleanDatNormNoColScaling))))) == length(removedRownames) & length(removedRownames) > 0) {
    cleanDatNormNoColScaling <- cleanDatNormNoColScaling[-match(removedRownames, rownames(cleanDatNormNoColScaling)), ]
    RW.relAbunFactors.HiMissRmvd <- RW.relAbunFactors[-match(removedRownames, names(RW.relAbunFactors))]
    cat(paste0("[iter_", repeats, "] Removed ", length(removedRownames), " high missingness rows. ", nrow(cleanDatNormNoColScaling), " rows remaining."))
  } else {
    cat(paste0("[iter_", repeats, "] No rows removed with >=50% missing values."))
    RW.relAbunFactors.HiMissRmvd <- RW.relAbunFactors
  }
  dim(cleanDatNormNoColScaling) # 2875 rows.
  
  
  
  prevIterCleanDatNorm2 <- data.frame(matrix(0, nrow = nrow(cleanDatNormNoColScaling), ncol = ncol(cleanDatNormNoColScaling)))
  if (exists("cleanDatNorm2")) {
    prevIterCleanDatNorm2 <- cleanDatNorm2
  }
  
  #-------------------------------------------------------------------------------
  ## Step 2, Enforce equal loading assumption on output of step 1b (all well-quantified proteins equally considered/weighted)
  colMeans(ratioCleanDatUnnorm, na.rm = TRUE)
  # should be zero, but they are usually not. E.g. here, we see b2 channels 130N, 130C, 131N are >5% away from even 1:1, and b4 channels 127N-128C are >10% away
  
  ## Set all column means to ~0 (10^-16 or less), essentially 0=log2(ratio/GIS) for all column means or an average ratio/GIS=1
  cleanDatNorm <- scale(cleanDatNormNoColScaling, scale = FALSE)
  colMeans(cleanDatNorm, na.rm = TRUE)
  
  ## alternative equivalent operation (if colAvg=mean(x,na.rm=TRUE))
  cleanDatNorm2 <- apply(cleanDatNormNoColScaling, 2, function(x) {
    colAvg <- eval(parse(text = paste0(meanORmedian, "(x,na.rm=TRUE)"))) ## MEAN/MEDIAN FUNCTION CHOICE***
    outputCol <- x - colAvg # rep(colAvg,length(x));
    outputCol
  })
  
  ## show columnwise normalization/scaling methods are equivalent (with rounding to a few decimals, at least)
  colMeans(cleanDatNorm, na.rm = TRUE) == colMeans(cleanDatNorm2, na.rm = TRUE) # TRUE if meanORmedian="mean"
  
  ## GIS-derived relative RW abundances can be applied back to the Step2 RW & CW-normalized data's rows.
  ## values of cleanDatNorm2 are log2(measurement/batch-stabilized GIS abundance).
  if (paste(names(RW.relAbunFactors.HiMissRmvd), collapse = ",") == paste(rownames(cleanDatNorm2), collapse = ",")) {
    relAbundanceNorm2 <- RW.relAbunFactors.HiMissRmvd * 2^cleanDatNorm2
  } else {
    cat(paste0("[iter_", repeats, "] ERROR: step 2 data with removed high missingness rows does not match relative abundance factors after trying to remove same rows.\n"))
  }
  
  
  #-------------------------------------------------------------------------------
  ## Step 2, revisited considering only top proteins by relative abundance rank and housekeeping (HKG) status
  
  ## Which among the top x genes are housekeeping? By use as loading controls...
  percentTopProteins <- 12.5 # 7.5-15 is reasonable for unsupervised; 15 percent ok if you choose your (small list of) HKG candidates manually
  topProteinCount <- as.integer(nrow(cleanDatNorm2) * percentTopProteins / 100)
  
  proteinsRankedByRelAbundance <- RW.relAbunFactors.HiMissRmvd[order(RW.relAbunFactors.HiMissRmvd, decreasing = TRUE)]
  proteinsRankedByRelAbundance[1:topProteinCount] # prior to v4, only top 60 were considered for manual check
  
  ## One could choose HKG candidates (for synaptosome preps) from above top 60 most abundant proteins list (fit 10/line)
  # HKGindex<-c(3,4,5,7,8,10,11,12,13,16,17,19,22,24,25,26,27,28,29,30,31:39,41,43,45,46,47,49,51,52,54,55,59,60) #promising known HKGs
  ## One could search for any HKG:
  # names(proteinsRankedByRelAbundance)[which(grepl("HPRT",names(proteinsRankedByRelAbundance)))] #not found
  # but we will use automatic HKG selection, so:
  HKGindex <- vector()
  
  ## Which are the top x genes by CV in the existing step 2 relative abundance normalized data?
  step2ProteinCVs <- apply(relAbundanceNorm2, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  normProteinsRankedByIncrCV <- step2ProteinCVs[order(step2ProteinCVs, decreasing = FALSE)]
  data.frame(CV.rankOrder = normProteinsRankedByIncrCV[1:topProteinCount], RelAbunRank = match(names(normProteinsRankedByIncrCV)[1:topProteinCount], names(proteinsRankedByRelAbundance)))
  
  # Must be in HKGindex from any hand-picked pre-specific or general HKGs AND in top [percentTopProteins]% by lowest CV
  if (length(HKGindex) > 0) {
    HKGindex <- intersect(HKGindex, na.omit(match(names(normProteinsRankedByIncrCV)[1:topProteinCount], names(proteinsRankedByRelAbundance)[1:topProteinCount])))
  } else {
    HKGindex <- as.vector(na.omit(match(names(normProteinsRankedByIncrCV)[1:topProteinCount], names(proteinsRankedByRelAbundance)[1:topProteinCount])))
  }
  
  ## Look at CV, CV rank, and Abundance rank of these candidate HKGs for redoing step 2 normalization using these robust rows only.
  HKGlist1 <- data.frame(Protein = names(proteinsRankedByRelAbundance)[HKGindex], CV.rank = match(names(proteinsRankedByRelAbundance)[HKGindex], names(normProteinsRankedByIncrCV)), CVpct = round(100 * normProteinsRankedByIncrCV[match(names(proteinsRankedByRelAbundance)[HKGindex], names(normProteinsRankedByIncrCV))], 1), RelAbunRank = HKGindex)
  HKGlist1 <- HKGlist1[order(HKGlist1$RelAbunRank, decreasing = FALSE), ]
  HKGlist1
  
  HKGlist2 <- rownames(HKGlist1)
  HKGindex2 <- match(HKGlist2, rownames(cleanDatNormNoColScaling))
  cleanDatNorm2.HKG <- apply(cleanDatNormNoColScaling, 2, function(x) {
    colAvg <- eval(parse(text = paste0(meanORmedian, "(x[HKGindex2],na.rm=TRUE)"))) ## MEAN/MEDIAN FUNCTION CHOICE***
    outputCol <- x - colAvg
    outputCol
  })
  
  cleanDat <- relAbundanceNorm2
  DFforFrobCurrent <- apply(cleanDatNorm2, 2, as.numeric)
  DFforFrobPrev <- apply(prevIterCleanDatNorm2, 2, as.numeric)
  removeColumnsCurrent <- which(apply(DFforFrobCurrent, 2, function(x) sum(is.na(x))) == nrow(DFforFrobCurrent))
  removeColumnsPrev <- which(apply(DFforFrobPrev, 2, function(x) sum(is.na(x))) == nrow(DFforFrobPrev))
  
  if (length(removeColumnsCurrent) > 0) {
    frobeniusNormCurrent <- norm(na.omit(DFforFrobCurrent[, -removeColumnsCurrent]), type = "F")
  } else {
    frobeniusNormCurrent <- norm(matrix((na.omit(DFforFrobCurrent))), type = "F")
  }
  if (length(removeColumnsPrev) > 0) {
    frobeniusNormPrev <- norm(na.omit(DFforFrobPrev[, -removeColumnsPrev]), type = "F")
  } else {
    frobeniusNormPrev <- norm(na.omit(DFforFrobPrev), type = "F")
  }
  initialFrobeniusNorm <- if (repeats == 1) {
    frobeniusNormCurrent
  } else {
    initialFrobeniusNorm
  }
  iterationTrackingDF[repeats, ] <- c(repeats, frobeniusNormCurrent, frobeniusNormPrev, frobeniusNormPrev - frobeniusNormCurrent, frobeniusNormCurrent / initialFrobeniusNorm)
  cat(paste0("       iteration convergence tracking (Frobenius Norm Difference):  ", signif(iterationTrackingDF[repeats, 4], 3), "\n"))
  if (abs(iterationTrackingDF[repeats, 4]) < FNthreshold) {
    cat("...Reached convergence criterion (Frobenius Norm Difference)<1e-8!\n")
    break
  }
} # closes "for (repeats in 1:iterations)"
## END NORMALIZATION

#-------------------------------------------------------------------------------
#' ## Plots to check the normalization progress.
#-------------------------------------------------------------------------------

iterations.intended <- iterations
iterations <- repeats

plot1 <- ggplot(data=iterationTrackingDF, aes(x=Iteration, y=abs(FrobenDiff), group=1)) +
  geom_line() + geom_point() + 
  geom_hline(yintercept = FNthreshold, linetype = "dashed", color = "red", size = 0.25) + 
  ggtitle("Iteration Tracking") + xlab("Iteration") + ylab("Frobenius Norm Diff from Previous") +
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold")
  )

plot2 <- ggplot(data=iterationTrackingDF, aes(x=Iteration, y=log10(abs(FrobenDiff)), group=1)) +
  geom_line() + geom_point() + 
  geom_hline(yintercept = log10(FNthreshold), linetype = "dashed", color = "red", size = 0.25) + 
  ggtitle("Iteration Tracking (log scale)") + xlab("Iteration") + ylab("Log10(Frobenius Norm Diff from Previous)") +
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold")
  )

# Figures
fig <- plot_grid(plot1,plot2, labels = "auto")
fig 

# Save figures as pdf. 
#file <- paste0(outputfigsdir,"/",tissue,"Iteration_Tracking.pdf")
#ggsavePDF(plots=list(plot1,plot2),file)

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "Iteration_Tracking.tiff")
ggsave(file,fig)

#-------------------------------------------------------------------------------
#' ## Remove any sample outliers.
#-------------------------------------------------------------------------------

# Data is...
data_in <- relAbundanceNorm2
data_in[1:5,1:5]

# Illustrate Oldham's sample connectivity.
sample_connectivity <- ggplotSampleConnectivityv2(data_in,log=TRUE,colID="b")
sample_connectivity$table
plot <- sample_connectivity$connectivityplot + ggtitle("Sample Connectivity post-TAMPOR")
plot

# Save as figure.
file <- paste0(outputfigsdir,"/",outputMatName,"TAMPOR_Outliers.tiff")
ggsave(file,plot, width = 3, height = 2.5, units = "in")

# Loop to identify Sample outliers using Oldham's connectivity method.
n_iter <- 5
threshold <- -3.0
out_samples <- list()

# Loop:
for (i in 1:n_iter){
  data_temp <- data_in
  oldham <- ggplotSampleConnectivityv2(data_temp,log=TRUE,colID="b",threshold = -3)
  bad_samples <- rownames(oldham$table)[oldham$table$Z.Ki < threshold]
  print(paste(length(bad_samples)," outlier sample(s) identified in iteration ",i,".",sep=""))
  if (length(bad_samples)==0) bad_samples <- "none"
  out_samples[[i]] <- bad_samples
  out <- grepl(paste(unlist(out_samples),collapse="|"),colnames(data_in))
  data_in <- data_in[,!out]
}

# Outlier samples.
bad_samples <- unlist(out_samples)
traits$Sample.Model[rownames(traits) %in% bad_samples]

# Save data to file. 
cleanDat <- data_in
dim(cleanDat)
datafile <- paste(Rdatadir,tissue,"TAMPOR_data_outliersRemoved.Rds",sep="/")
saveRDS(cleanDat,datafile)

#-------------------------------------------------------------------------------
#' ## Examine sample clustering with MDS and PCA post-TAMPOR Normalization. 
#-------------------------------------------------------------------------------

# Split into cortex and striatum plots. 

# Insure that any outlier samples have been removed.
traits <- traits[rownames(traits) %in% colnames(cleanDat),]

# Check, traits and cleanDat should match data.
all(rownames(traits) == colnames(cleanDat))

## PCA Plots.
# Relative abundance. 
colors <- traits$Color
traits$ColumnName <- rownames(traits)

# Cortex and striatum.
idx <- traits$Tissue=="Cortex"
idy <- traits$Tissue=="Striatum"
plot1 <- ggplotPCA(log2(cleanDat[,idx]), traits, colors[idx], title = "Cortex")
plot2 <- ggplotPCA(log2(cleanDat[,idy]), traits, colors[idy], title = "Striatum")

# plots
plot1
plot2

# Save as figure.
file <- paste0(outputfigsdir,"/",outputMatName,"Cortex_TAMPOR_PCA.tiff")
ggsave(file,plot1, width = 3, height = 3, units = "in")

file <- paste0(outputfigsdir,"/",outputMatName,"Striatum_TAMPOR_PCA.tiff")
ggsave(file,plot2, width = 3, height = 3, units = "in")

# Customize colors.
#colors <- c("#FFF200", "#FFF200",    #22B14C
#            "#00A2E8", "#00A2E8",    #A349A4
#            "#22B14C", "#22B14C",  #FFF200
#            "#A349A4", "#A349A4")    #00A2E8
#plot1 + scale_color_manual(values = colors) + theme(legend.position = "none")
#plot2 + scale_color_manual(values = colors) + theme(legend.position = "none")

# Figure
fig <- plot_grid(plot1,plot2)
fig 

## MDS Plots.
# Relative abundance. 
#plots <- ggplotMDSv2(log2(cleanDat), colID="b", traits, title = "Normalized Abundance")
#plot1 <- plots$plot + theme(legend.position = "none")
#plot2 <- plots$dendro
#plot1

# Save as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_PCA_Post_TAMPOR.tiff")
ggsave(file,fig)


#-------------------------------------------------------------------------------
#' ## Statistical comparisons post-TAMPOR.
#-------------------------------------------------------------------------------

# Data is...
file <- paste(Rdatadir,tissue,"TAMPOR_data_outliersRemoved.Rds",sep="/")
cleanDat <- readRDS(file)

# Insure traits are loaded.
traitsfile <- paste(Rdatadir,tissue,"Combined_Cortex_Striatum_traits.Rds",sep="/")
sample_info <- readRDS(traitsfile)

# Create DGEList object...
data <- cleanDat
data[1:5,1:5]
y_DGE <- DGEList(counts = data)

# Example, checking the the normalization with plotMD.
plotMD(cpm(y_DGE, log=TRUE), column=2)
abline(h=0, col="red", lty=2, lwd=2)

# Create sample mapping.
traits <- sample_info
traits$ColumnName <- traits$SampleID
traits <- subset(traits,rownames(traits) %in% colnames(data))
traits <- traits[match(colnames(data),rownames(traits)),]
all(traits$ColumnName==colnames(data))
group <- paste(traits$Tissue,traits$Sample.Model,sep=".")
unique(group)
group[grepl("Cortex.WT",group)] <- "Cortex.WT"
group[grepl("Striatum.WT",group)] <- "Striatum.WT"
unique(group)
y_DGE$samples$group <- as.factor(group)

# Basic design matrix for GLM.
design <- model.matrix(~0 + group, data=y_DGE$samples)
colnames(design) <- levels(y_DGE$samples$group)
design

# Estimate dispersion:
y_DGE <- estimateDisp(y_DGE,design,robust=TRUE)

# PlotBCV
plot <- ggplotBCV(y_DGE)
plot 

# Fit a general linear model.  
fit <- glmQLFit(y_DGE, design, robust = TRUE)

# Examine the QL fitted dispersion.
plotQLDisp(fit)

# Generate contrasts.
g1 <- colnames(design)[grepl("Cortex",colnames(design))][-5]
g2 <- colnames(design)[grepl("Striatum",colnames(design))][-5]

cont1 <- makePairwiseContrasts(list(g1),list("Cortex.WT"))
cont2 <- makePairwiseContrasts(list(g2),list("Striatum.WT"))

# For some reason loops or lapply dont work with the makeContrasts function.
contrasts <- list(
  makeContrasts(cont1[1], levels=design),
  makeContrasts(cont1[2], levels=design),
  makeContrasts(cont1[3], levels=design),
  makeContrasts(cont1[4], levels=design),
  makeContrasts(cont2[1], levels=design),
  makeContrasts(cont2[2], levels=design),
  makeContrasts(cont2[3], levels=design),
  makeContrasts(cont2[4], levels=design))

# Call glmQLFTest() to evaluate differences in contrasts.
qlf <- lapply(contrasts, function(x) glmQLFTest(fit,contrast=x))

## Determine number of significant results with decideTests().
summary_table <- lapply(qlf, function(x) summary(decideTests(x)))
overall <- t(matrix(unlist(summary_table),nrow=3,ncol=8))
rownames(overall) <- unlist(lapply(contrasts,function(x) colnames(x)))
colnames(overall) <- c("Down","NotSig","Up")
overall <- as.data.frame(overall)
overall <- add_column(overall,Contrast=rownames(overall),.before=1)
overall <- overall[,c(1,3,2,4)]
overall$TotalSig <- rowSums(overall[,c(3,4)])

# Table of DE candidates.
table <- tableGrob(overall, rows = NULL)
grid.arrange(table)

# Save table as tiff.
file <- paste0(outputfigsdir, "/", outputMatName, "_TAMPOR_DE_Table.tiff")
ggsave(file,table)

# Call topTags to add FDR. Gather tablurized results. 
results <- lapply(qlf, function(x) topTags(x, n = Inf, sort.by = "none")$table)

# Function to annotate DE candidates:
annotateTopTags <- function(y_TT){
  y_TT$logCPM <- 100*(2^y_TT$logFC)
  colnames(y_TT)[2] <- "%WT"
  colnames(y_TT)[3] <- "F Value"
  y_TT$candidate <- "no"
  y_TT[which(y_TT$FDR <= 0.10 & y_TT$FDR > 0.05), dim(y_TT)[2]] <- "low"
  y_TT[which(y_TT$FDR <= 0.05 & y_TT$FDR > 0.01), dim(y_TT)[2]] <- "med"
  y_TT[which(y_TT$FDR <= 0.01), dim(y_TT)[2]] <- "high"
  y_TT$candidate <- factor(y_TT$candidate, levels = c("high", "med", "low", "no"))
  return(y_TT)
}

# Rows should just be Uniprot for working with some of my old functions.
my_func <- function(x){
  rownames(x) <- sapply(strsplit(rownames(x),"\\|"),"[",2)
  return(x)
}
results <- lapply(results,function(x) my_func(x))

# Convert logCPM column to percent WT. 
# Annotate with candidate column.  
results <- lapply(results, function(x) annotateTopTags(x))

# Annotate with Gene names and Entrez IDS.
results <- lapply(results,function(x) annotate_Entrez(x))
names(results) <- c(
  sapply(strsplit(cont1," - "),"[",2),
  sapply(strsplit(cont2," - "),"[",2))

# Sort by pvalue.
results <- lapply(results,function(x) x[order(x$PValue),])
names(results) <- sapply(strsplit(overall$Contrast," "),"[",1)

# Write to excel.
file <- paste0(outputtabsdir,"/",outputMatName,"_TAMPOR_GLM_Results.xlsx")
write.excel(results,file)

# Write to RDS.
file <- paste0(Rdatadir,"/",outputMatName,"_TAMPOR_GLM_Results.RDS")
saveRDS(results,file)

#-------------------------------------------------------------------------------
#' ## Volcano plots for cortex and striatum.
#-------------------------------------------------------------------------------

# Add column for genotype and unique ID to results in list.
for (i in 1:length(results)){
  Experiment <- names(results)[i]
  df <- results[[i]]
  df <- add_column(df,Experiment, .before = 1)
  #ID <- paste(df$Experiment,df$Uniprot, sep = "_")
  #df <- add_column(df, ID, .after = 1)
  results[[i]] <- df
}

# Merge the results. 
df <- do.call(rbind,results)

# Add column for genotype specific colors.
colors <- as.list(c("#FFF200","#00A2E8","#22B14C","#A349A4"))
#colors <- as.list(c("yellow","blue","green","purple"))
names(colors) <- c("Shank2","Shank3","Syngap1","Ube3a")
df$Color <- unlist(colors[sapply(strsplit(df$Experiment,"\\."),"[",3)])

# Split into cortex and striatum datasets.
idx <- grepl("Cortex",df$Experiment)
results <- split(df,idx)
names(results) <- c("Striatum","Cortex")

# Function for producing volcano plots.
vp <- function(df){
  df$x <- df[,grep("FC",colnames(df))]
  df$y <- -log10(df[,grep("PValue",colnames(df))])
  logic <- df$FDR < 0.05
  df$Color[!logic] <- "gray"
  df$Color <- as.factor(df$Color)
  y_int <- -1*log10(max(df$PValue[df$FDR<0.05]))
  plot <- ggplot(data = df, aes(x = x, y = y, color = Color)) + 
    geom_point() + scale_color_manual(values=levels(df$Color)) +
    geom_hline(yintercept = y_int, linetype = "dashed", color = "black", size = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.6) +
    #geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 0.6) + 
    xlab("Log2(Fold Change) ASD vs Control") + ylab("-Log10PValue") + 
    theme(
      plot.title = element_text(hjust = 0.5, color="black", size=11, face="bold"),
      axis.title.x = element_text(color="black", size=11, face="bold"),
      axis.title.y = element_text(color="black", size=11, face="bold"),
      legend.position = "none")
  return(plot)
}

# Cortex plot.
plot1 <- vp(results$Cortex) + ggtitle("Cortex") + xlim(-3.5,3.5)
plot1

# Summarize up and down-regulated proteins.
updown <- function(results, tissue){
  res <- results[[tissue]]
  up <- sum(res$FDR < 0.05 & res$logFC > 0)
  down <- sum(res$FDR < 0.05 & res$logFC < 0)
  return(list(up=up,down=down))
}

updown(results,"Cortex")
updown(results,"Striatum")

# Save
file <- paste0(outputfigsdir,"/",outputMatName,"Cortex_VolcanoPlot.tiff")
ggsave(file,plot1, width = 3, height = 2.5, units = "in")

# Striatum plot.
plot2 <- vp(results$Striatum) + ggtitle("Striatum") + xlim(-3,3)
plot2

# Save
file <- paste0(outputfigsdir,"/",outputMatName,"Striatum_VolcanoPlot.tiff")
ggsave(file,plot2, width = 3, height = 2.5, units = "in")

#-------------------------------------------------------------------------------
#' ## Volcano plots for each genotype.
#-------------------------------------------------------------------------------

# Add column for genotype and unique ID to results in list.
for (i in 1:length(results)){
  Experiment <- names(results)[i]
  df <- results[[i]]
  df <- add_column(df,Experiment, .before = 1)
  #ID <- paste(df$Experiment,df$Uniprot, sep = "_")
  #df <- add_column(df, ID, .after = 1)
  results[[i]] <- df
}

# Merge the results. 
df <- do.call(rbind,results)

# Add column for genotype specific colors.
colors <- as.list(c("#FFF200","#00A2E8","#22B14C","#A349A4"))
names(colors) <- c("Shank2","Shank3","Syngap1","Ube3a")
df$Color <- unlist(colors[sapply(strsplit(df$Experiment,"\\."),"[",3)])

# Split by genotype (color).
results <- split(df,df$Color)
names(results) <- names(colors)[match(names(results),colors)]

df <- results$Shank3

# Function for producing volcano plots.
ggplotVolcanoPlot2 <- function(df){
  df$x <- df[,grep("FC",colnames(df))]
  df$y <- -log10(df[,grep("PValue",colnames(df))])
  logic <- df$FDR < 0.05
  df$Color[!logic] <- "gray"
  df$Color <- as.factor(df$Color)
  y_int <- -1*log10(max(df$PValue[df$FDR<0.05]))
  plot <- ggplot(data = df, aes(x = x, y = y, color = Color)) + 
    geom_point() + scale_color_manual(values=levels(df$Color)) +
    geom_hline(yintercept = y_int, linetype = "dashed", color = "black", size = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.6) +
    #geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 0.6) + 
    xlab("Log2(Fold Change) ASD vs Control") + ylab("-Log10PValue") + 
    theme(
      plot.title = element_text(hjust = 0.5, color="black", size=11, face="bold"),
      axis.title.x = element_text(color="black", size=11, face="bold"),
      axis.title.y = element_text(color="black", size=11, face="bold"),
      legend.position = "none")
  return(plot)
}

# Generate plots.
plots <- lapply(as.list(names(colors)),function(x) ggplotVolcanoPlot2(results[[x]]))
names(plots) <- names(colors)

# Add titles.
for (i in 1:length(plots)){
  plots[[i]] <- plots[[i]] + ggtitle(names(plots)[i])
}

# Make x-axis symmetrical.
for (i in 1:length(plots)){
  plot <- plots[[i]]
  xlim <- max(abs(layer_scales(plot)$x$range$range))
  plots[[i]] <- plot + xlim(-xlim,xlim)
}

plots$Shank2
plots$Shank3
plots$Syngap1
plots$Ube3a

# Save plots.
for (i in 1:length(plots)){
  plot <- plots[[i]]
  file <- paste0(outputfigsdir,"/",outputMatName,names(plots)[i],"_VolcanoPlot.tiff")
  ggsave(file,plot, width = 3, height = 2.5, units = "in")
}

# Summarize up and down-regulated proteins.
updown <- function(results, genotype){
  res <- results[[genotype]]
  res$UpDown <- "Up" 
  res$UpDown[res$logFC<0] <- "Down"
  res$Experiment <- paste(res$Experiment,res$UpDown)
  sub <- subset(res,res$FDR<0.05)
  return(  table(sub$Experiment))
}

rbind(
  updown(results,"Shank2"),
  updown(results,"Shank3"),
  updown(results,"Syngap1"),
  updown(results,"Ube3a")
)

#-------------------------------------------------------------------------------
#' ## GO Enrichment analysis of DEPs.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# GO testing with the EdgeR qlf object...
# Extract rownames and map to entrez ids.
uniprot <- sapply(strsplit(rownames(qlf[[1]]),"\\|"),"[",2)
symbols <- sapply(strsplit(rownames(qlf[[1]]),"\\|"),"[",1)

# Map Uniprot IDs to Entrez IDs:
entrez <- mapIds(org.Mm.eg.db, 
                 keys = uniprot, 
                 column = "ENTREZID", 
                 keytype = "UNIPROT", 
                 multiVals = "first")

# Number of missing entrez ids.
sum(is.na(entrez))

# Map unmapped gene symbols to Entrez IDs:
entrez[is.na(entrez)] <- mapIds(org.Mm.eg.db, 
                                keys = symbols[is.na(entrez)], 
                                column = "ENTREZID", 
                                keytype = "SYMBOL", 
                                multiVals = "first")

# Number of remaining missing entrez ids.
sum(is.na(entrez))

# Perform GO and KEGG enrichment testing (using proteins with FDR<0.1).
GO <- lapply(qlf, function(x) goana(x, geneid = entrez, species="Mm", FDR=0.1))
KEGG <- lapply(qlf, function(x) kegga(x, geneid = entrez, species="Mm", FDR=0.1))

# Add names.
contrasts <- unlist(lapply(qlf,function(x) x$comparison))
contrasts <- gsub(" -1", "",sapply(strsplit(contrasts,"\\*"),"[",2))
names(GO) <- contrasts
names(KEGG) <- contrasts

# Calculate adjusted p.values.
my_func <- function(x){
  x$P.adj.Up <- p.adjust(x$P.Up, method = "bonferroni")
  x$P.adj.Down <- p.adjust(x$P.Down, method = "bonferroni")
  x$FDR.Up <- p.adjust(x$P.Up, method = "BH")
  x$FDR.Down <- p.adjust(x$P.Down, method = "BH")
  return(x)
}

GO <- lapply(GO, function(x) my_func(x))
KEGG <- lapply(KEGG, function(x) my_func(x))

# Combined result
GO_result <- do.call(rbind,GO)
GO_result <- add_column(
  GO_result,
  class = sapply(strsplit(rownames(GO_result),".GO"),"[",1),
  .before = 1
)
GO_result <- add_column(
  GO_result,
  GOID = sapply(strsplit(rownames(GO_result),"\\."),"[",4),
  .after = 1
)
rownames(GO_result) <- NULL

# Subset sig results.
sigGO <- lapply(GO,function(x) subset(x,FDR.Up < 0.1 | x$FDR.Down < 0.1))

# Write to excel.
file <- paste0(outputtabsdir,"/",outputMatName,"_TAMPOR_GO_Results.xlsx")
write.excel(GO,file)

file <- paste0(outputtabsdir,"/",outputMatName,"_TAMPOR_sigGO_Results.xlsx")
write.excel(sigGO,file)

#-------------------------------------------------------------------------------
#' ## GO enrichment analaysis for DEPs from each genotype.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Load the GLM statistical results.
file <- paste(outputtabs,"Final_TAMPOR",
              "Combined_TMT_Analysis_TAMPOR_GLM_Results.xlsx", sep = "/")
results <- lapply(as.list(c(1:8)),function(x) read_excel(file,x))
names(results) <- excel_sheets(file)

# Build a df with statistical results.
stats <- lapply(results,function(x) 
  data.frame(Uniprot=x$Uniprot,
             Gene = x$Gene,
             FDR=x$FDR))
names(stats) <- names(results)
df <- stats %>% reduce(left_join, by = c("Uniprot","Gene"))
colnames(df)[c(3:ncol(df))] <- names(stats)

## Prepare a matrix of class labels (colors) to pass to enrichmentAnalysis().
labels <- data.frame(
  Shank2 = df$Cortex.KO.Shank2 < 0.05 | df$Striatum.KO.Shank2 < 0.05,
  Shank3 = df$Cortex.KO.Shank3 < 0.05 | df$Striatum.KO.Shank3 < 0.05,
  Syngap1 = df$Cortex.HET.Syngap1 < 0.05 | df$Striatum.HET.Syngap1 < 0.05,
  Ube3a3 = df$Cortex.KO.Ube3a < 0.05 | df$Striatum.KO.Ube3a < 0.05
)
rownames(labels) <- df$Uniprot

# Convert TRUE to column names. 
logic <- labels == TRUE # 1 will become TRUE, and 0 will become FALSE.
# Loop through each column to replace 1 with column header (color).
for (i in 1:ncol(logic)){
  col_header <- colnames(labels)[i]
  labels[logic[,i],i] <- col_header
}

# Map Uniprot IDs to Entrez.
# Get Uniprot IDs and gene symbols from rownames
#Uniprot_IDs <- as.character(rownames(colors))

# Go back and get Symbol|Uniprot rownames.
genes <- as.list(sapply(strsplit(rownames(cleanDat),"\\|"),"[",1))
names(genes) <- sapply(strsplit(rownames(cleanDat),"\\|"),"[",2)

uniprot <- rownames(labels)

# Map Uniprot IDs to Entrez IDs:
entrez <- mapIds(org.Mm.eg.db, 
                 keys = uniprot, 
                 column = "ENTREZID", 
                 keytype = "UNIPROT", 
                 multiVals = "first")

# Map unmapped IDs...
not_mapped <- unlist(genes[uniprot[is.na(entrez)]])
entrez[is.na(entrez)] <- mapIds(org.Mm.eg.db, 
                                keys = not_mapped, 
                                column = "ENTREZID", 
                                keytype = "SYMBOL", 
                                multiVals = "first")
sum(is.na(entrez))

# Insure that labels is a matrix.
labels <- as.matrix(labels)
head(labels)

# look at the number of genes assigned to each cluster. 
table(labels)

# The labels matrix and vector of cooresponding entrez IDs 
# will be passed to enrichmentAnalysis().

# Build a GO annotation collection:
if (!exists(deparse(substitute(musGOcollection)))){
  musGOcollection <- buildGOcollection(organism = "mouse")}
  
# Creates some space by clearing some memory.
collectGarbage()
  
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
  useBackground = "given", # options are: given, reference (all), intersection, and provided. 
  threshold = 0.05,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "FALSE")

# Create some space by clearing some memory.
collectGarbage()

# Collect the results. 
results_GOenrichment <- list()
for (i in 1:length(GOenrichment$setResults)){
  results_GOenrichment[[i]] <- GOenrichment$setResults[[i]]$enrichmentTable
}
length(results_GOenrichment)
names(results_GOenrichment) <- colnames(labels)

a <- results_GOenrichment[[1]]
b <- results_GOenrichment[[2]]
c <- results_GOenrichment[[3]]
d <- results_GOenrichment[[4]]

# Write results to file. 
file <- paste0(outputtabsdir,"/",tissue,"_WGCNA_Analysis_Module_GOenrichment_Results.xlsx")
write.excel(results_GOenrichment,file)

#-------------------------------------------------------------------------------
#' ## Generate GO scatter plots for each genotype.
#-------------------------------------------------------------------------------

# Change names of GOenrichment results to genotype specific colors.
data <- results_GOenrichment
names(data) <- c("yellow","blue","green","purple")

plots <- list(
  Shank2  = ggplotGOscatter(data, color = "yellow", topN = 1),
  Shank3  = ggplotGOscatter(data, color = "blue", topN = 1),
  Syngap1 = ggplotGOscatter(data, color = "green", topN = 1),
  Ube3a   = ggplotGOscatter(data, color = "purple", topN = 1)
)

# Save plots and legends.
for (i in 1:length(plots)){
  plot <- plots[[i]]
  file <- paste0(outputfigsdir,"/",outputMatName,"_",
                 names(plots)[[i]],"_GO_ScatterPlot.tiff")
  ggsave(file, plot, width = 3.0, height = 2.5, units = "in")
}

#-------------------------------------------------------------------------------
#' ## Condition overlap.
#-------------------------------------------------------------------------------

# Load statistical results..
file <- paste0(Rdatadir,"/",outputMatName,"_TAMPOR_GLM_Results.RDS")
results <- readRDS(file)

# Check.
lapply(results,function(x) sum(x$FDR<0.05))

# Combine by FDR. 
stats <- lapply(results,function(x) 
  data.frame(Uniprot=x$Uniprot, Gene = x$Gene, FDR=x$FDR))
names(stats) <- names(results)
df <- stats %>% reduce(left_join, by = c("Uniprot","Gene"))
colnames(df)[c(3:ncol(df))] <- paste("FDR",names(stats),sep=".")
rownames(df) <- paste(df$Gene,df$Uniprot,sep="|")
df$Uniprot <- NULL
df$Gene <- NULL

# Check.
apply(df,2, function(x) sum(x<0.05))

# Gather sigProts.
sigProts <- list()
for (i in 1:ncol(df)){
  idx <- df[,i]<0.05
  sigProts[[i]] <- rownames(df)[idx]
}
names(sigProts) <- names(stats)

# Build a matrix showing overlap.
col_names <- names(stats)
row_names <- names(stats)

# All possible combinations.
contrasts <- expand.grid(col_names,row_names)
contrasts$GenoA <- as.vector(contrasts$GenoA)
contrasts$GenoB <- as.vector(contrasts$GenoB)
colnames(contrasts) <- c("GenoA","GenoB")

# Loop
int <- list()
for (i in 1:dim(contrasts)[1]){
  a <- unlist(as.vector(sigProts[contrasts$GenoA[i]]))
  b <- unlist(as.vector(sigProts[contrasts$GenoB[i]]))
  int[[i]] <- intersect(a,b)
}

# Add to contrasts.
contrasts$Intersection <- unlist(lapply(int,function(x) length(x)))

# Make overlap matrix.
dm <- matrix(contrasts$Intersection,nrow=8,ncol=8)
rownames(dm) <- colnames(dm) <- row_names

# Remove upper tri and melt.
dm[lower.tri(dm)] <- NA

# Calculate percent overlap.
dm2 <- sweep(dm,1,apply(dm, 1, function(x) max(x, na.rm=TRUE)), FUN = "/")

# Melt
df <- melt(dm,na.rm = TRUE)
df$percent <- round(melt(dm2, na.rm = TRUE)$value,2)

# Generate plot.
plot <- ggplot(df, aes(Var2, Var1, fill = percent)) +
  geom_tile(color = "white") + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 2.5) +
  scale_fill_gradient2(name="Percent Overlap") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.7),
    legend.direction = "horizontal") + 
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  coord_fixed()

# Extract legend and save seperately.
legend <- get_legend(plot)
plot <- plot + theme(legend.position = "none")

# Save figure.
file <- paste0(outputfigsdir,"/",outputMatName,"Condition_Overlap.tiff")
ggsave(file,plot, width = 3, height = 3, units = "in")

# Save legend.
file <- paste0(outputfigsdir,"/",outputMatName,"Condition_Overlap_legend.tiff")
ggsave(file,legend, width = 3, height = 3, units = "in")

#-------------------------------------------------------------------------------
#' ## Condition overlap, combined tissue types. 
#-------------------------------------------------------------------------------

# Load statistical results..
file <- paste0(Rdatadir,"/",outputMatName,"_TAMPOR_GLM_Results.RDS")
results <- readRDS(file)

# Check.
lapply(results,function(x) sum(x$FDR<0.05))

# Combine by FDR. 
stats <- lapply(results,function(x) 
  data.frame(Uniprot=x$Uniprot, Gene = x$Gene, FDR=x$FDR))
names(stats) <- names(results)
df <- stats %>% reduce(left_join, by = c("Uniprot","Gene"))
colnames(df)[c(3:ncol(df))] <- paste("FDR",names(stats),sep=".")
rownames(df) <- paste(df$Gene,df$Uniprot,sep="|")
df$Uniprot <- NULL
df$Gene <- NULL

# Check.
apply(df,2, function(x) sum(x<0.05))

# Gather sigProts.
sigProts <- list()
for (i in 1:ncol(df)){
  idx <- df[,i]<0.05
  sigProts[[i]] <- rownames(df)[idx]
}
names(sigProts) <- names(stats)

# Combine sigprots for each genotype.
genos <- c("Shank2","Shank3","Syngap1","Ube3a")
u <- list()
for (geno in genos){
  print(geno)
  idx <- grep(geno,names(sigProts))
  u[[geno]] <- union(sigProts[[idx[1]]],sigProts[[idx[2]]])
}
sigProts <- u

# Build a matrix showing overlap.
col_names <- names(sigProts)
row_names <- names(sigProts)

# All possible combinations.
contrasts <- expand.grid(col_names,row_names)
colnames(contrasts) <- c("GenoA","GenoB")
contrasts$GenoA <- as.vector(contrasts$GenoA)
contrasts$GenoB <- as.vector(contrasts$GenoB)

# Loop
int <- list()
for (i in 1:dim(contrasts)[1]){
  a <- unlist(as.vector(sigProts[contrasts$GenoA[i]]))
  b <- unlist(as.vector(sigProts[contrasts$GenoB[i]]))
  int[[i]] <- intersect(a,b)
}

contrasts$Name <- paste(contrasts$GenoA,contrasts$GenoB,sep="_U_")
names(int) <- contrasts$Name
int$Syngap1_U_Ube3a

# Add to contrasts.
contrasts$Intersection <- unlist(lapply(int,function(x) length(x)))

# Make overlap matrix.
dm <- matrix(contrasts$Intersection,nrow=4,ncol=4)
rownames(dm) <- colnames(dm) <- row_names

# Remove upper tri and melt.
dm[lower.tri(dm)] <- NA

# Calculate percent overlap.
dm2 <- sweep(dm,1,apply(dm, 1, function(x) max(x, na.rm=TRUE)), FUN = "/")

# Melt
df <- melt(dm,na.rm = TRUE)
df$percent <- round(melt(dm2, na.rm = TRUE)$value,2)

# Generate plot.
plot <- ggplot(df, aes(Var2, Var1, fill = percent)) +
  geom_tile(color = "white") + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 2.5) +
  scale_fill_gradient2(name="Percent Overlap") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.7),
    legend.direction = "horizontal") + 
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  coord_fixed()

plot

# Extract legend and save seperately.
legend <- get_legend(plot)
plot <- plot + theme(legend.position = "none")

# Save figure.
file <- paste0(outputfigsdir,"/",outputMatName,"Condition_Overlap.tiff")
ggsave(file,plot, width = 3, height = 3, units = "in")

# Save legend.
file <- paste0(outputfigsdir,"/",outputMatName,"Condition_Overlap_legend.tiff")
ggsave(file,legend, width = 3, height = 3, units = "in")

#-------------------i------------------------------------------------------------
#' ## Generate protein boxplots for significantly DE proteins.
#-------------------------------------------------------------------------------
#+ eval = FALSE

# Load statistical results..
file <- paste0(Rdatadir,"/",outputMatName,"_TAMPOR_GLM_Results.RDS")
results <- readRDS(file)

# Expression data. 
exprDat <- log2(y_DGE$counts)

# Insure all WT samples are annotated as WT in traits.
traits_temp <- traits
traits_temp$Sample.Model <- paste(traits_temp$Tissue,traits_temp$Sample.Model,sep=".")
traits_temp$Sample.Model[grepl("Cortex.WT", traits_temp$Sample.Model)] <- "Cortex.WT"
traits_temp$Sample.Model[grepl("Striatum.WT", traits_temp$Sample.Model)] <- "Striatum.WT"
traits_temp <- subset(traits_temp, rownames(traits_temp) %in% colnames(exprDat))

# Generate plots.
plot_list <- ggplotProteinBoxes(
  data_in = exprDat,
  interesting.proteins = rownames(exprDat),
  dataType = "Relative Abundance",
  traits = traits_temp,
  order = c(2,7,4,6,5,8,1,9,3,10),
  scatter = TRUE
)

# Add custom colors.
colors <- c("gray","gray",
            "#FFF200", "#FFF200",
            "#00A2E8", "#00A2E8",
            "#22B14C", "#22B14C",
            "#A349A4","#A349A4")

plot_list <- lapply(plot_list, function(x) x + scale_fill_manual(values = colors))

# Example plot.
plot_list[[1]]

## Add significance stars.
# Build a df with statistical results.
stats <- lapply(results, function(x)
  as.data.frame(cbind(Uniprot = x$Uniprot, FDR = x$FDR)))
names(stats) <- names(results)
df <- stats %>% reduce(left_join, by = "Uniprot")
colnames(df)[c(2:ncol(df))] <- names(stats)
head(df)

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

# Modify x-axis labels-- significant bar axis labels are in red.
a <- c("black","black","red","red","black","black","black","black","black","black")
b <- c("black","black","black","black","red","red","black","black","black","black")
c <- c("black","black","black","black","black","black","red","red","black","black")
d <- c("black","black","black","black","black","black","black","red","red","red")

p1 <- p1 + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
p2 <- p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = b))
p3 <- p3 + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = c))
p4 <- p4 + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = d))

plots <- list(p1,p2,p3,p4)
names(plots) <- c("Shank2","Shank3","Syngap1","Ube3a")

# Save to pdf.
#file <- paste0(outputfigsdir, "/", tissue, "_WGCNA_Analysis_InterBatch_ProteinBoxPlots.pdf")
#ggsavePDF(plots = plot_list, file)

# Save as tiff.
for (i in 1:length(plots)){
  file <- paste0(outputfigsdir, "/", outputMatName, names(plots)[i],"_BoxPlot.tiff")
  ggsave(file,plots[[i]], width = 3, height = 3, units = "in")
}


## Inspecting other interesting proteins...

# Make function to change x-axis text to red if significant.
change_axis_color <- function(plot){
  build <- ggplot_build(plot)
  labs <- build$data[[3]]
  idx <- grep("\\*",labs$label[order(labs$x)]) + 2
  color_vect <- rep("black",10)
  color_vect[idx] <- "red"
  plot <- plot + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = color_vect))
  return(plot)
}


# Other interesting proteins.
plots <- plot_list[int$Syngap1_U_Ube3a]
names(plots) <- sapply(strsplit(names(plots),"\\|"),"[",1)

# Change axis color.
plots <- lapply(plots,function(x) change_axis_color(x))

# Save as tiff.
for (i in 1:length(plots)){
  file <- paste0(outputfigsdir, "/", outputMatName, names(plots)[i],"_BoxPlot.tiff")
  ggsave(file,plots[[i]], width = 3, height = 3, units = "in")
}

#-------------------------------------------------------------------------------
#' ## Render report.
#-------------------------------------------------------------------------------
#' This script is formatted for automated rendering of an RMarkdown report.
#'
#+ eval = FALSE

# Code directory. 
dir <- "D:/Documents/R/TMT-Analysis/Synaptosome-TMT-Analysis/Code/Bradshaw"
file <- paste(dir,"TMT_Analysis_v14.R",sep="/")

# Save and render.
rstudioapi::documentSave()
rmarkdown::render(file)
