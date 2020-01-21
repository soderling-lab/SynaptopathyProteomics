#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change:
net <- "Cortex" # Which network are we analyzing? One of: c("Cortex","Striatum")

# Global options and imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(getPPIs)
  library(purrr)
  library(WGCNA)
  library(org.Mm.eg.db)
  library(anRichment)
  library(getPPIs)
  library(DescTools)
  library(igraph)
  library(vegan)
})

# Directories.
if (rstudioapi::isAvailable()) {
	setwd("D:/projects/SynaptopathyProteomics/analysis/3_Network-Analysis")
}
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
figsdir <- file.path(root, "figs")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load GLM stats.
myfile <- file.path(rdatdir, "2_GLM_Stats.RData")
glm_stats <- readRDS(myfile)

# Proteins with any significant change.
sigProts <- apply(glm_stats$FDR,1,function(x) any(x<0.05))
sigProts <- names(sigProts)[sigProts]

# SigProts by genotype.
fdr_df <- glm_stats$FDR
fdr_df$Protein <- rownames(fdr_df)
sig_df <- melt(fdr_df,id="Protein") %>% 
	group_by(variable) %>% filter(value < 0.05) %>% group_split()
sigProts_geno <- sapply(sig_df,function(x) x$Protein)
names(sigProts_geno) <- gsub(" ","_", gsub(" FDR", "", colnames(fdr_df)[1:8]))

# Load expression data.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_cleanDat.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_cleanDat.RData")
)
data <- t(readRDS(myfiles[net])) # Data should be transposed: rows, proteins.

# Load Sample info.
traits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load co-expression (adjacency) matrices.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_Adjm.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_Adjm.RData")
)
adjm <- as.matrix(readRDS(myfiles[net]))
rownames(adjm) <- colnames(adjm)

# Load network partitions-- self-preservation enforced.
ids <- c("Cortex"="14942508","Striatum"="14940918")
myfile <- list.files(rdatdir, pattern = ids[net], full.names = TRUE)
partitions <- readRDS(myfile) 

# Reset partition index.
partitions <- lapply(partitions, reset_index)

#------------------------------------------------------------------------------
## Evaluate similarity between modules from select partitions.
#------------------------------------------------------------------------------

# Name partitions.
names(partitions) <- paste0("R",c(1:length(partitions)))

# Get all modules.
named_partitions <- lapply(partitions,function(x) {
			      m = split(x,x)
			      names(m) = paste0("M",names(m))
			      return(m)
})
all_modules <- unlist(named_partitions,recursive=FALSE)
modules <- names(all_modules)

# All comparisons (contrasts).
contrasts <- expand.grid("M1"=modules,"M2"=modules,stringsAsFactors=FALSE)

# Examine module jaacard similarity for all comparisons between 
# representative modules. This is overkill since we wont be using 
# comparisions. But the resulting matrix is easy to work with.

# Loop:
message("Calculating Module Jaacard Similarity...")
n <- nrow(contrasts)
modulejs <- vector("numeric",n)
for (i in 1:n){
	if (i==1) { pbar <- txtProgressBar(min=1,max=n,style=3) }
	setTxtProgressBar(pbar,i)
	x <- contrasts[i,]
	idm1 <- x[["M1"]]
	idm2 <- x[["M2"]]
	m1 <- names(all_modules[[idm1]])
	m2 <- names(all_modules[[idm2]])
	if (idm1 == idm2) {
		modulejs[i] <- 1
	} else {
		modulejs[i] <- js(m1,m2)
	}
	if (i == n) { close(pbar); message("\n") }
} # Ends loop.

# Save.
myfile <- file.path(rdatdir,paste0("3_",net,"_All_Module_JS.RData"))
saveRDS(modulejs,myfile)
