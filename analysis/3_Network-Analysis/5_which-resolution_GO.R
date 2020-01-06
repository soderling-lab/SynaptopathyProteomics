#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Set-up the workspace.
#-------------------------------------------------------------------------------

# User parameters to change.
net <- "Cortex" 
#mypart <- c(Cortex = "10360847",Striatum = "10342568")[net]
mypart <- c(Cortex = "10773682",Striatum = "10781799")[net] # relaxed criterion

# Global options and imports.
suppressPackageStartupMessages({
	library(org.Mm.eg.db)
	library(anRichment)
	library(data.table)
	library(getPPIs)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Functions.
myfun <- list.files(funcdir, full.names = TRUE)
invisible(sapply(myfun, source))

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load expression data.
myfile <- file.path(rdatdir, paste0("3_",net,"_cleanDat.RData"))
data <- readRDS(myfile)

# Load network partitions-- self-preservation enforced.
myfile <- list.files(rdatdir,pattern=mypart,full.names=TRUE)
partitions <- readRDS(myfile)

#------------------------------------------------------------------------------
## Build SynGO gene collection.
#------------------------------------------------------------------------------

# Load SynGO annotations.
# Data downloaded from: https://syngoportal.org/
myfile <- file.path(rdatdir,"SynGO_bulk_download_release_20180731",
		   "syngo_annotations.xlsx")
synGO <- readxl::read_excel(myfile)

# Load SynGO gene mapping table.
myfile <- file.path(rdatdir,"SynGO_bulk_download_release_20180731",
		   "syngo_genes.xlsx")
genes <- readxl::read_excel(myfile)

# Some rows contain multiple MGI ids, seperate these.
genes <- tidyr::separate_rows(gene_map, mgi_id,sep=",")

# Map MGI ids to mouse entrez.
mgi <- paste0("MGI:",genes$mgi_id)
entrez <- mapIDs(mgi,from="mgi",to="entrez",species="mouse")
names(entrez) <- genes$mgi_id
genes$mus_entrez <- entrez

# Map Human HGNC ids to mouse entrez.
idx <- match(synGO$"human ortholog gene hgnc_id",genes$hgnc_id)
synGO$mus_entrez <- genes$mus_entrez[idx]

# Remove rows with unmapped genes.
synGO <- subset(synGO,!is.na(synGO$mus_entrez))

# Collect as named list of genes.
mus_entrez <- synGO$mus_entrez
data_list <- split(mus_entrez,synGO$"GO term ID")

# Loop to build gene sets:
geneSets <- list()
for (i in 1:length(data_list)) {
	id <- names(data_list)[i]
	geneSets[[i]] <- newGeneSet(geneEntrez = data_list[[i]],
				    geneEvidence = "IEA", # Inferred from Electronic Annotation
				    geneSource = "SynGO",
				    ID = id, 
				    name = id,
				    description = "Synaptic gene ontology",
				    source = "https://syngoportal.org/data/download.php?file=SynGO_bulk_download_release_20180731.zip",
				    organism = "mouse",
				    internalClassification = "SynGO",
				    groups = "PL",
				    lastModified = "2020-01-03")
}
# Annotate with group name.
SynGOgroup = newGroup(name = "SynGO", 
		   description = "Currated synaptic gene ontology from SynGO database.",
		   source = "syngoportal.org")

# Combine as gene collection.
SynGOcollection <- newCollection(dataSets = geneSets, groups = list(SynGOgroup))

# Save as Rdata.
myfile <- file.path(rdatdir,"3_SynGOcollection.RData")
saveRDS(SynGOcollection,myfile)

## Combine SynGO with all other mouse GO data.

# Build mouse GO collection:
musGOcollection <- buildGOcollection(organism="mouse")

# Which GO groups would you like to use in your analysis?
keep <- c("GO","GO.BP","GO.MF","GO.CC")
musGOcollection <- subsetCollection(musGOcollection, tags = keep)

# Combine SynGO and GO datasets.
GOcollection <- newCollection()
GOcollection <- addToCollection(musGOcollection,SynGOcollection)

#------------------------------------------------------------------------------
## Perform GO analysis of modules at every resolution.
#------------------------------------------------------------------------------

# What about disease enrichment...
#myfile <- file.path(rdatdir,"mouse_DisGeneNETcollection.RData")
#GOcollection <- readRDS(myfile)

# Loop to perform GO enrichment for modules at every resolution.
message("Performing GO enrichment analysis...")
GOresults <- list()
#for (i in seq_along(partitions)) {
for (i in c(1,50,100)) {
  # Initialize progress bar.
  if (i == 1) {
    pb <- txtProgressBar(min = 0, max = length(partitions), style = 3)
  }
  # Perform GO analysis.
  GOresults[[i]] <- moduleGOenrichment(partitions,i, protmap,GOcollection)
  # Update progress bar.
  setTxtProgressBar(pb, i)
  # Close pb.
  if (i == length(partitions)) {
	  close(pb)
	  message("\n")
  }
} # Ends loop.

# Save results.
#myfile <- file.path(rdatdir,paste0("3_",net,"_Module_GO_Results.RData")) 
#saveRDS(GOresults, myfile)

data = GOresults[[50]]
write_excel(data,"temp.xlsx")


s = sapply(GOresults,function(x) sum(sapply(x,function(y) any(y$FDR <0.05))))

s == max(s)

sum(sapply(data,function(x) any(x$FDR <0.05)))
length(data)


sappl



#------------------------------------------------------------------------------
## Examine GO results in order to define ~best resolution.
#------------------------------------------------------------------------------

# Load GO results.
GOresults <- readRDS(myfile)

# Examine biological enrichment of modules at every resolution.
# Summarize the biological significance of a resolution as the sum of 
# -log(GO pvalues) for all modules.
modSig <- lapply(GOresults, function(x) 
		 sapply(x, function(y) sum(-log(y$pValue))))

modSig <- lapply(GOresults, function(x) 
		 sapply(x, function(y) mean(-log(y$pValue))))

# The code above is confusing, this is what it does:
#x = results[[1]] # list of go enrichment for all modules at res 1.
#y = x[[1]] # go enrichment df of module 1.
#y$pValue
#-log(y$pValue)
#sum(-log(y$pValue))
#sapply(x,function(y) sum(-log(y$pValue)))
#lapply(results, function(x) sapply(x,function(y) sum(-log(y$pValue))))
#out = lapply(results, function(x) sapply(x,function(y) sum(-log(y$pValue))))

# Summarize every resolution.
resSum <- sapply(modSig, sum)
best_res <- c(1:length(resSum))[resSum == max(resSum)]
names(best_res) <- net
# Status report. 
message(paste("Best resolution based on GO enrichment:",best_res))

# Save results.
myfile <- file.path(rdatdir,paste0("3_",net,"_GO_Best_Resolution.RData"))
saveRDS(best_res, myfile)

#------------------------------------------------------------------------------
# Define a function to compare partitions.

p1 = partitions[[90]]
p2 = partitions[[100]]
#module_assignment_similarity(p1,p2) # 

module_assignment_similarity <- function(p1,p2,ignore=0){

	# Check, partitions should have the same names!
	if (!identical(names(p1),names(p2))){
		stop()
	}
	# Set modules to be ignored to NA.
	p1[p1 == ignore] <- NA
	p2[p2 == ignore] <- NA
	# Create module assignment df.
	df <- data.table::CJ("ProtA"=names(p1),"ProtB"=names(p1), unique = TRUE)
	df <- df[!df$ProtA==df$ProtB,] # Remove self-interactions.
	df$ProtA.P1.Module <- p1[df$ProtA]
	df$ProtB.P1.Module <- p1[df$ProtB]
	df$ProtA.P2.Module <- p2[df$ProtA]
	df$ProtB.P2.Module <- p2[df$ProtB]
	# Check if ProtA and ProtB are in same modules in each partition (Pmk).
	df$Pmk1 <- p1[df$ProtA] ==  p1[df$ProtB]
	df$Pmk2 <- p2[df$ProtA] ==  p2[df$ProtB]
	# How to handle NA?
	df$Pmk1[is.na(df$Pmk1)] <- FALSE # 
	df$Pmk2[is.na(df$Pmk2)] <- FALSE
	# Calculate similarity statistic.
	# Euclidean inner (dot) product/Euclidean norm.

	dp <- pracma::dot(as.numeric(df$Pmk1),as.numeric(df$Pmk2))
	en/
		sqrt(sum(df$Pmk1^2)) * sqrt(sum(df$Pmk2^2))


	d <- dist(rbind(df$Pmk1,df$Pmk2))
	return(d)
}
