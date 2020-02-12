#!/usr/bin/env Rscript

here = getwd()
root = dirname(dirname(here))
rdatdir = file.path(root,"rdata")

# Load GO semantic similarity partitions.
myfile = file.path(rdatdir,"3_Combined_GOSemSim_partitions.csv")
gosim_parts = data.table::fread(myfile,drop=1,skip=1)

# Load GO semantic similarity adjm.
myfile = file.path(rdatdir,"3_All_GO_SemSim.RData")
gosemsim_adjm = readRDS(myfile)

# Load protein identifier map.
myfile = file.path(rdatdir,"2_Protein_ID_Map.RData")
prot_map = readRDS(myfile)

# Split go sem sim partition df into list of partitions.
# First, fix names.
colNames = prot_map$ids[match(colnames(gosemsim_adjm),prot_map$entrez)]
colnames(gosim_parts) = colNames
parts_list = lapply(c(1:100),function(x) {
	      part = as.numeric(gosim_parts[x,])
	      names(part) = colNames
	      return(part)
})

# Load Cortex partitions.
myfile = list.files(rdatdir,patter="10773682",full.names=TRUE)
cox_parts = readRDS(myfile)

# Compare partitions.
fmi = vector(mode = "numeric", length = 100)
for (i in 1:100) {
	p0 = cox_parts[[i]]
	p2 = parts_list[[i]]
	p1 = p0[c(1:length(p0))[names(p0) %in% names(p2)]]
	rho = sum(p1 == p2)/length(p1)
	#rho = cor(p1,p2,method = "spearman")
	# Check.
	if (!all(names(p1) == names(p2))) { stop() }
	# Calculate Folkes Mallow similarity index.
	fmi[i] = abs(rho)
	#fmi[i] = dendextend::FM_index_R(p1,p2)[1]
}

rbest = c(1:length(fmi))[fmi == max(fmi)]
table(cox_parts[[rbest]])

# vector Pmk of length Nk(Nk−1)/2
# where Nk is the number of genes in the network. 
# Each element of this vector corresponds to a pair of genes and equals
# one if the two genes are in the same module and zero otherwise. Accordingly, for
# any two module predictions (method m1 applied to network k1, and method m2
# applied to network k2), we calculated the distance as follows:


head(p1)
head(p2)

dendextend::FM_index_R(p1,p2)[1]


