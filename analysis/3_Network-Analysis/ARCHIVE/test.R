#!/usr/bin/env Rscript

# Imports.
library(GOSemSim)

# Directories.
here = getwd()
root = dirname(dirname(here))
rdatdir = file.path(root,"rdata")

# Load Cortex partitions.
myfile = list.files(rdatdir,patter="10773682",full.names=TRUE)
cox_parts = readRDS(myfile)

# Load protein identifier map.
myfile = file.path(rdatdir,"2_Protein_ID_Map.RData")
prot_map = readRDS(myfile)

# Get Entrez IDs.
entrez <- prot_map$entrez[match(names(cox_parts[[1]]),prot_map$ids)]

# Build list of GO ontologies.
msGO <- lapply(c("BP","CC","MF"), function(x) {
		       godata('org.Mm.eg.db', ont=x)
})
names(msGO) <- c("BP","CC","MF")

# Calculate semantic similarity between all pairwise comparisons of genes.
gosemsim <- lapply(msGO,function(x) { 
	       mgeneSim(genes=head(entrez), semData=x, measure="Wang",verbose=TRUE)
})

# Save to file.
myfile <- file.path(rdatdir,"3_GO_Semantic_Similarity_Adjms.RData")
saveRDS(gosemsim,myfile)

# Combine BP + MF + CC results: RMS.

# Cluster.
