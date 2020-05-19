#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Parse command line arguments:
if (interactive()) {
	## If interactive:
	# User defined parameters (you only need to change these two):
	analysis_type = "Striatum" # Tissue type for analysis.
} else if (!interactive()) {
	## If not interactive, check that only 1 arg is passed.
	args <- commandArgs(trailingOnly=TRUE)
	if (length(args) == 1) { 
		analysis_type = commandArgs(trailingOnly=TRUE)[1]
	} else { 
		stop("Specify either 'Cortex' or 'Striatum'.",call.=FALSE) 
	}
}

## DEFAULT project root.
root = "/mnt/d/projects/SynaptopathyProteomics"

# Input data should be in root/rdata/:
input_data = list("Cortex" = "Cortex_norm_protein.csv",
		  "Striatum" = "Striatum_norm_protein.csv")[[analysis_type]]

## Output for downstream analysis:
output_name = analysis_type

# Outputs stored in root/rdata/
# * [output_name]_Adjm.csv
# * [output_name]_NE_Adjm.csv

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------

# Load renv.
renv::load(root)

# Imports.
suppressPackageStartupMessages({
  library(dplyr)
  library(WGCNA)
  library(neten)
  library(getPPIs)
  library(data.table)
})

# Additional functions.
TBmiscr::load_all()

# Directories.
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#---------------------------------------------------------------------
## Prepare the data.
#---------------------------------------------------------------------

# Load final normalized data.
myfile <- file.path(rdatdir,input_data)
prot_dt <- fread(myfile)

# Log transform data.
prot_dt$Abundance <- log2(prot_dt$Intensity)

# Coerce to data matrix.
dm <- prot_dt %>% 
	dcast(Accession ~ Sample, value.var="Abundance") %>%
	as.matrix(rownames="Accession")

#---------------------------------------------------------------------
## Create co-expression networks.
#---------------------------------------------------------------------

message(paste("\nCreating",analysis_type,
	      "protein co-variation networks."))

# Create signed adjacency (correlation) matrices.
adjm <- WGCNA::bicor(t(dm))

# Perform network enhancement.
message("\nPerforming network enhancement.")
adjm_ne <- neten(adjm)

# Write correlation matrices to file.
myfile <- file.path(rdatdir, paste0(output_name,"_Adjm.csv"))
adjm %>% as.data.table(keep.rownames="Accession") %>% fwrite(myfile)

# Write enhanced networks to file.
myfile <- file.path(rdatdir, paste0(output_name,"_NE_Adjm.csv"))
adjm_ne %>% as.data.table(keep.rownames="Accession") %>% fwrite(myfile)

#---------------------------------------------------------------------
## Generate the PPI graph.
#---------------------------------------------------------------------

# Load mouse interactome.
data("musInteractome")

# Subset mouse interactome, keep data from mouse, human, and rat.
orgs <- c(10090, 9606, 10116)
ppis <- musInteractome %>% filter(Interactor_A_Taxonomy %in% orgs)

# Get entrez IDs for all proteins in co-expression networks.
all_entrez <- unique(prot_dt$Entrez)

# Build a ppi graph.
message("\nBuilding PPI graph...")
g <- buildNetwork(ppis, all_entrez, taxid = 10090)

# Get ppi adjacency matrix.
PPI_adjm <- as.matrix(as_adjacency_matrix(g))

# Map entrez back to uniprot accession.
idx <- match(colnames(PPI_adjm),prot_dt$Entrez)
colnames(PPI_adjm) <- prot_dt$Accession[idx]

# Evaluate scale free fit.
node_degree <- apply(PPI_adjm, 2, sum)
fit <- WGCNA::scaleFreeFitIndex(node_degree, nBreaks = 10, removeFirst = FALSE)
r2 <- fit$Rsquared.SFT
message(paste("\nScale free fit of PPI graph:", round(r2, 3)))

#--------------------------------------------------------------------
## Demonstrate that interacting proteins are highly co-expressed.
#--------------------------------------------------------------------

# Randomly sample 10,000 edges drawn from interacting and 
# non-interacting proteins.
#n <- 10000

# Seed seed for reproducibility.
#set.seed(as.numeric(Sys.time()))

# Get random samples.
#adjm_edges <- reshape2::melt(adjm,value.name="bicor") %>% as.data.table()
#colnames(adjm_edges)[c(1,2)] <- c("ProtA","ProtB")
#adjm_edges$ProtA <- as.character(adjm_edges$ProtA)
#adjm_edges$ProtB <- as.character(adjm_edges$ProtB)

# Fix names.
#rownames(PPI_adjm) <- colnames(PPI_adjm)
#ppi_edges <- reshape2::melt(PPI_adjm,value.name="ppi") %>% as.data.table()
#colnames(ppi_edges)[c(1,2)] <- c("ProtA","ProtB")
#ppi_edges$ProtA <- as.character(ppi_edges$ProtA)
#ppi_edges$ProtB <- as.character(ppi_edges$ProtB)

# Combine, yeah its big.
#edges_dt <- left_join(adjm_edges,ppi_edges,by=c("ProtA","ProtB")) %>%
#	as.data.table()
#edges_dt$ppi[is.na(edges_dt$ppi)] <- 0 # Coerce NA to 0.

# Get random sample.
#df <- as.data.frame(edges_dt)
#idx <- sample(nrow(df),n)
#subdf <- df[idx,]
# Check the mean of each group.
#subdf %>% group_by(ppi) %>% summarize(mean(bicor))

# Calculate WRS p-value.
# FIXME: Permutation test would be more robust.
#WRS_test <- wilcox.test(subdf$bicor ~ subdf$ppi, alternative = "less")
#message(paste("\nWRS P-value:",round(WRS_test$p.value,4)))

# Save PPI adjm to file.
myfile <- file.path(rdatdir,paste0(output_name,"_PPI_Adjm.csv"))
PPI_adjm %>% as.data.table(keep.rownames="Accession") %>% 
	fwrite(myfile)

#--------------------------------------------------------------------
## Save the ppi data.
#--------------------------------------------------------------------

# Subset PPIs.
idx_a <- ppis$osEntrezA %in% all_entrez
idx_b <- ppis$osEntrezB %in% all_entrez
idx <- idx_a & idx_b
PPI_dt <- subset(ppis,idx) %>% as.data.table()

# Add Symbol annotation.
PPI_dt$SymbolA <- prot_dt$Symbol[match(PPI_dt$osEntrezA,prot_dt$Entrez)]
PPI_dt$SymbolB <- prot_dt$Symbol[match(PPI_dt$osEntrezB,prot_dt$Entrez)]

# Save to file.
myfile <- file.path(tabsdir,paste0(analysis_type,"_PPI_data.csv"))
fwrite(PPI_dt,myfile)

message("Done!\n")
