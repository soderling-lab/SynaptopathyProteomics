#!/usr/bin/env Rscript

# title: Analysis of Synaptopathy TMT Proteomics
# description: map Uniprot Acccession to Entrez IDs and gene Symbols
# author: Tyler W Bradshaw

## ---- INPUTs

# project root
root <- "~/projects/SynaptopathyProteomics"

# the TMT data in root/rdata:
data_files <- c(
  "Cortex" = "Cortex_Raw_Peptides.csv",
  "Striatum" = "Striatum_Raw_Peptides.csv"
)


## ---- R ENVIRONMENT 

# Prepare the R workspace for the analysis. Load custom functions and prepare
# the project directory for saving output files.

# load required packages
suppressPackageStartupMessages({
  library(dplyr) # for working with data
  library(geneLists) # twesleyb/geneLists
  library(data.table) # for working with tables
})

# load project specific data and functions
devtools::load_all()


## ---- MAIN 

myfiles <- file.path(root,"inst","extdata", data_files)
data_list <- lapply(myfiles, fread) 
proteins <- unique(unlist(lapply(data_list,function(x) x$Accession)))

## ---- Create gene identifier map 

# get Entrez ids with MGI database
all_entrez <- geneLists::queryMGI(proteins)

entrez <- all_entrez[!is.na(all_entrez)]

# map entrez to gene symbols with twesleyb/geneLists::getIDs()
symbols <- geneLists::getIDs(entrez, from = "entrez", to = "symbol", species = "mouse")

# create gene map
gene_map <- as.data.table(keep.rownames = "uniprot", entrez)
gene_map$symbol <- symbols[as.character(gene_map$entrez)]


## ---- Save output data 

# save gene_map as rda object in root/data
myfile <- file.path(root,"data", "gene_map.rda")
save(gene_map, file = myfile, version = 2)
message("saved: ", myfile)

# save some other genes
#shank2 <- gene_map$uniprot[match("Shank2",gene_map$symbol)]
#shank3 <- gene_map$uniprot[match("Shank3",gene_map$symbol)]
#syngap1 <- gene_map$uniprot[match("Syngap1",gene_map$symbol)]
#ube3a <- gene_map$uniprot[match("Ube3a",gene_map$symbol)]

saveGene <- function(gene, datadir, quiet=TRUE) {
  uniprot <- gene_map$uniprot[match(gene,gene_map$symbol)]
  myfile <- file.path(datadir,paste0(tolower(gene),".rda"))
  assign(tolower(gene),value=uniprot)
  save(list=tolower(gene), file=myfile, version=2)
  if (!quiet) { message("saved: ", myfile) }
}

# saves rogdi.rda (rogdi <- getUniprot(rogdi))
saveGene("Rogdi",datadir=file.path(root,"data"), quiet=FALSE)

genes <- c("Shank2","Shank3","Sygap1","Ube3a","Anks1b", "Gad1", "Insyn1")

invisible(sapply(genes, saveGene, datadir=file.path(root,"data"),quiet=FALSE))
