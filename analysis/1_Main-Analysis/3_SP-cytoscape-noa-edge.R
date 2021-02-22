#!/usr/bin/env Rscript

# author: twab
# description: create edge and noa files for Cytoscape
# title: SynaptopathyProteomics

## ---- prepare the R env

root <- "~/projects/SynaptopathyProteomics"
devtools::load_all(root, quiet=TRUE)


## ---- load the data

#library(SynaptopathyProteomics)
data(gene_map)
data(cortex)
data(cortex_results)
data(corum_prots)

library(dplyr)

library(getPPIs) # twesleyb/getPPIs
data(musInteractome) # mouse PPIs from HitPredict


## ---- create node annotation (noa) data frame

noa <- gene_map %>% 
	filter(uniprot %in% results$Protein) %>% 
	dplyr::mutate(Symbol = toupper(symbol))


#library(geneLists)

#data(corum)
#names(corum) <- namen

data(corum_prots)

# replace - and ' ' with _
# remove other special char ( []()./ , )
namen <- gsub("\\[|\\]|\\(|\\)|,|\\/|\\.","",gsub("-|\\ ", "_", names(corum_prots)))
names(corum_prots) <- namen

# munge to create df with rows == genes, cols = corum complex, 
# val = TRUE, FALSE (prot is in complex or not)
mylist <- lapply(noa$uniprot, function(prot) lapply(corum_prots, function(path) prot %in% path))
names(mylist) <- as.character(noa$uniprot)
df <- mylist %>% bind_rows(.id="uniprot")

# combine with noa
noa <- noa %>% left_join(df, by="uniprot")


sig_prots <- unique(results$Protein[results$FDR<0.05])

noa$sigProt <- noa$uniprot %in% sig_prots

# save
myfile <- file.path(root,"rdata","noa.csv")
data.table::fwrite(noa, myfile)
message("saved: ", myfile)


## ---- collect ppis

# keep ppis from human, mouse, and rat
os_keep <- c(9606, 10116, 10090) # taxonomix identifiers

# all entrez ids 
entrez <- noa$entrez

# collect interactions between swip and wash_interactome proteins
ppis <- musInteractome %>%
	filter(Interactor_A_Taxonomy %in% os_keep) %>%
	filter(Interactor_B_Taxonomy %in% os_keep) %>%
	subset(osEntrezA %in% entrez & osEntrezB %in% entrez)
edge_df <- ppis %>% select(osEntrezA, osEntrezB, Publications)

# entrez IDs to UniProt
idx <- match(edge_df$osEntrezA,gene_map$entrez)
protA <- gene_map$uniprot[idx]
idy <- match(edge_df$osEntrezB,gene_map$entrez)
protB <- gene_map$uniprot[idy]

# add to table
edge_df <- tibble::add_column(edge_df, protA, .before='osEntrezA')
edge_df <- tibble::add_column(edge_df, protB, .after='protA')

# save
myfile <- file.path(root,"rdata","edges.csv")
data.table::fwrite(edge_df, myfile)
message("saved: ", myfile)
