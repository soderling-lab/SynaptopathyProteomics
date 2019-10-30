#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
## Set-up the workspace.
#-------------------------------------------------------------------------------

# Global options and imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(readxl)
	library(igraph)
	library(reshape2)
	library(ggplot2)
	library(anRichment)
	library(TBmiscr)
})

# Directories.
here <- getwd()
root <- dirname(dirname(here))
datadir <- file.path(root,"data")
tabsdir <- file.path(root,"tables")
funcdir <- file.path(root,"functions")

# Load functions.
source(file.path(funcdir,"0_Functions.R"))

# Load expression data
wtDat <- t(readRDS(file.path(datadir,"wtDat.Rds")))
koDat <- t(readRDS(file.path(datadir,"koDat.Rds")))

# Fix names.
colnames(wtDat) <- colnames(koDat) <- rownames(readRDS(file.path(datadir,"wtDat.Rds")))

# Compute adjmatrix:
wtAdjm <- quiet(WGCNA::bicor(wtDat))
koAdjm <- quiet(WGCNA::bicor(koDat))

# ~Best resolutions. 
# Most biological (GO) information: WT = 44; KO = 30.
r_wt <- 44
r_ko <- 30

# Load WT partitions.
wtParts <- data.table::fread(file.path(datadir,"WT_partitions.csv"),drop=1)
wtPartition <- as.integer(wtParts[r_wt,]) + 1
names(wtPartition) <- colnames(wtParts)
wtModules <- split(wtPartition,wtPartition)
message(paste("WT modules:", length(wtModules), "(before enforcing module self-preservation)."))

# Load KO partitions.
koParts <- data.table::fread(file.path(datadir,"KO_partitions.csv"),drop=1)
koPartition <- as.integer(koParts[r_ko,]) + 1
names(koPartition) <- colnames(koParts)
koModules <- split(koPartition,koPartition)
message(paste("KO modules:", length(koModules), "(before enforcing module self-preservation)."))

# Checks:
if (!all(colnames(wtDat) == colnames(koDat))) { stop("Input data don't match!") }
if (!all(colnames(wtAdjm) == colnames(koAdjm))) { stop("Input data don't match!") }
if (!all(names(wtPartition) %in% colnames(wtDat))) { stop("Input data don't match!") }
if (!all(names(koPartition) %in% colnames(koDat))) { stop("Input data don't match!") }

#-------------------------------------------------------------------------------
## Perform permutation testing for module self-preservation.
#-------------------------------------------------------------------------------

# Input for NetRep:
data_list        <- list(wt = wtDat,  ko = koDat)   
correlation_list <- list(wt = wtAdjm, ko = koAdjm) 
network_list     <- list(wt = wtAdjm, ko = koAdjm) 
module_list      <- list(wt = wtPartition, ko = koPartition)

# Perform permutation test for module self-preservation.
self <- as.list(c("wt","ko"))
selfPreservation <- lapply(self,function(x) {
			       NetRep::modulePreservation(
							  network = network_list,
							  data = data_list,
							  correlation = correlation_list,
							  moduleAssignments = module_list,
							  modules = NULL,
							  backgroundLabel = 0,
							  discovery = x,
							  test = x,
							  selfPreservation = TRUE,
							  nThreads = 8,
							  #nPerm = 100000, 
							  null = "overlap",
							  alternative = "greater",
							  simplify = TRUE,
							  verbose = TRUE)
})

## Remove modules that are not strongly preserved.
# Remove modules that are not strongly preserved--a module is not preserved if 
# any of its module preservation statistic adjusted p-values exceed 0.05.

# Get maximum p-value for each module's preservation stats. Corrected for 
# n module comparisons.
maxp <- function(preservation) {
	p <- apply(preservation$p.values,1,function(x) max(x,na.rm=TRUE))
	q <- p.adjust(p,"bonferroni")
	return(q)
}
q <- lapply(selfPreservation, maxp)

# Modules with NS preservation stats. 
out <- lapply(q,function(x)names(x)[x>0.05])

# For NS modules, set module membership to 0.
wtPartition[wtPartition %in% out[[1]]] <- 0
koPartition[koPartition %in% out[[2]]] <- 0

# Check module assignments.
wtModules <- split(wtPartition,wtPartition)
message(paste("WT modules:", length(wtModules), "(after enforcing module self-preservation)."))

koModules <- split(koPartition,koPartition)
message(paste("KO modules:", length(koModules), "(after enforcing module self-preservation)."))

# Save.
out <- list(wt = wtPartition, ko = koPartition)
saveRDS(out,file.path(datadir,"3_preserved_partitions.Rds"))

#-------------------------------------------------------------------------------
## Utilize permutation approach to identify divergent modules.
#-------------------------------------------------------------------------------

# Input for NetRep:
data_list        <- list(wt = wtDat,  ko = koDat) 
correlation_list <- list(wt = wtAdjm, ko = koAdjm) 
# Network list - THIS IS WHAT IS USED TO CALCULATE average edge weight not ADJM!
network_list     <- list(wt = wtAdjm,  ko = koAdjm)
module_list      <- list(wt = wtPartition, ko = koPartition)

# Generalize for discovery/test.
h0 = list(wt = c(discovery = "wt", test = "ko"), 
	  ko = c(discovery = "ko", test = "wt"))

# Perform permutation testing.
preservation <- lapply(h0, function(x) NetRep::modulePreservation(
					   network = network_list,
					   data = data_list,
					   correlation = correlation_list,
					   moduleAssignments = module_list,
					   modules = NULL,
					   backgroundLabel = 0,
					   discovery = x["discovery"],
					   test = x["test"],
					   selfPreservation = FALSE,
					   nThreads = 8,
					   #nPerm = 100000,  # determined by the function.
					   null = "overlap",
					   alternative = "two.sided", # c(greater,less,two.sided)
					   simplify = TRUE,
					   verbose = TRUE
					   )
)

## Identify divergent modules. 
# Divergent modules are those whose observed correlation structure is significantly
# less than the null model. 

# Just avg.weight.
q <- lapply(preservation, function(x) p.adjust(x$p.values[,1],"bonferroni"))

# Require all to be less than 0.05.
q <- lapply(preservation, function(x) p.adjust(apply(x$p.values,1,max)))

sigModules <- lapply(q,function(x) names(x)[x<0.05])

# Status.
sigWT <- sigModules$wt
sigKO <- sigModules$ko
message(paste("Number of WT modules that exhibit divergence:", length(sigWT)))
message(paste("Number of KO modules that exhibit divergence:", length(sigKO)))

# KO down: 5,6,13
# No WT Modules down...

#------------------------------------------------------------------------------
## Examine observed versus null distributions.
#------------------------------------------------------------------------------

plot_distributions <- function(x){
	module_results <- list()	
	for (module in 1:length(x$propVarsPresent)) {
		plots <- list()
		for (permstat in 1:dim(x$observed)[2]) {
			obs <- x$observed[module,permstat]
			p <- x$p.values[,permstat]
			statistic <- colnames(x$observed)[permstat]
			q <- round(p.adjust(p,"bonferroni")[module],3)
			mytitle <- paste0(statistic,"\n (p.adj = ",q,")")
			nulls <- data.frame("null"=x$nulls[module,permstat,])
			plot <- ggplot(data=nulls, aes(null)) + geom_histogram(bins=100,fill="gray") + 
				geom_vline(xintercept=obs, color = "red") + 
				ggtitle(mytitle) + 
				theme(plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
				      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
			              axis.title.y = element_text(color = "black", size = 11, face = "bold"))
			# annotate with observed value.
			ymaximum <- unlist(max(ggplot_build(plot)$layout$panel_params[[1]]$y.range))
			plot <- plot + annotate("text", x = obs, y = (ymaximum + 0.2*ymaximum), label = round(obs,3))
			plots[[statistic]] <- plot
		}
		module_results[[module]] <- plots
	}
	return(module_results)
}


# Generate plots.
plots <- list(wt = plot_distributions(preservation$wt),
	      ko = plot_distributions(preservation$ko))

# Save plots for average edge strength.
save_plots <- function(x) {
	for (i in 1:length(x)){
		namen <- paste0("M",i,"_avg_weight.tiff")
		ggsave(namen, plot = x[[i]][[1]])
	}
}

save_plots(plots$wt)
save_plots(plots$ko)

# WT modules...
p <- apply(preservation$wt$p.values,1,max)
q <- p.adjust(p,method="bonferroni")
length(wtModules)
sum(q<0.05)

# KO modules...
p <- apply(preservation$ko$p.values,1,max)
q <- p.adjust(p,method="bonferroni")
length(koModules)
sum(q<0.05)

#------------------------------------------------------------------------------
## Examine divergent modules. 
#------------------------------------------------------------------------------
# Under the null hypothesis that nothing is changing, modules with 
# average edge weight greater than or less than the null distribution are changing.

# Load statistical results.
myfile <- file.path(tabsdir, "2_Supplementary_TMT_GLM_Results.xlsx")
results <- lapply(as.list(c(1:8)),function(x) read_excel(myfile,x))
names(results) <- excel_sheets(myfile)

# Build a df with statistical results.
stats <- lapply(results, function(x)
		data.frame(Uniprot=x$Uniprot,
			   Symbol=x$Gene,
			   FDR=x$FDR))

names(stats) <- names(results)
statsdf <- stats %>% purrr::reduce(left_join, by = c("Uniprot","Symbol"))
colnames(statsdf)[c(3:ncol(statsdf))] <- names(stats)

# Proteins with any sig change.
statsdf$sigProt <- apply(statsdf,1,function(x) any(as.numeric(x[c(3:ncol(statsdf))])<0.05))

# Load protein identifier map for mapping protein names to entrez.
protmap <- data.table::fread(file.path(datadir,"ProtMap.csv"))

# Insure rownames are gene|uniprot.
rownames(statsdf) <- protmap$prots[match(as.character(statsdf$Uniprot),protmap$uniprot)]

# 1. Generate heat map.
generate_heatmaps <- function(modules) {
	for (i in 1:length(modules)) {
	     prots <- names(modules[[i]])
	     idx <- idy <- colnames(wtAdjm) %in% prots
	     subWT <- wtAdjm[idx,idy]
	     idx <- idy <- colnames(koAdjm) %in% prots
	     subKO <- koAdjm[idx,idy]
	     ## Generate Heatmap.
	     # Function to Reorder correlation matrix.
	     # Uses correlation between variables as distance.
	     reorder_cormat <- function(cormat){
		     dd <- as.dist((1-cormat)/2)
		     hc <- hclust(dd)
		     cormat <- cormat[hc$order, hc$order]}
	     # Reorder the correlation matrix based on WT values.
	     cormat <- reorder_cormat(subWT)
	     # Replace half of the correlation matrix with KO values.
	     cormat[upper.tri(cormat)] <- subKO[upper.tri(subKO)]
	     # Melt the correlation matrix
	     melted_cormat <- melt(cormat, na.rm = TRUE)
	     colnames(melted_cormat) <- c("WT","KO","value")
	     # Generate Heatmap
	     namen <- paste0("M",names(modules)[i])
	     plot <- ggplot(data = melted_cormat, aes(WT, KO, fill = value)) +
		     geom_tile(color = "white") +
		     scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
					  limit = c(min(melted_cormat$value),max(melted_cormat$value)), 
					  space = "Lab", name="Bicor") + theme_minimal() + 
	     theme(axis.text.x = element_blank(),
		   axis.text.y = element_blank(),
		   plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
	           axis.title.x = element_text(color = "black", size = 11, face = "bold"),
	           axis.title.y = element_text(color = "black", size = 11, face = "bold")) +
	     coord_fixed() + ggtitle(namen)
     # Save. 
     ggsave(paste0(namen,"heatmap.tiff"),plot) 
	}
}

# This will generate heatmaps and save them to file.
generate_heatmaps(wtModules)
generate_heatmaps(koModules)

# 2. Examine changes in edge weight...
generate_corplots <- function(modules) {
for (i in 1:length(modules)) {
	     prots <- names(modules[[i]])
	     idx <- idy <- colnames(wtAdjm) %in% prots
	     subWT <- wtAdjm[idx,idy]
	     idx <- idy <- colnames(koAdjm) %in% prots
	     subKO <- koAdjm[idx,idy]
	     df <- data.frame(protA = colnames(subWT),
		 protB = rep(colnames(subWT),each=ncol(subWT)),
		 wt = reshape2::melt(subWT,na.rm=TRUE)$value,
		 ko = reshape2::melt(subKO,na.rm=TRUE)$value)
	# Remove self-interactions.
	df <- df[!df$protA == df$protB,] 
	df$delta <- df$wt - df$ko
	# Sort.
	df <- df[order(df$delta,decreasing=TRUE),]
	# Get the top few prots.
	for (n in seq(2,10,by=2)) {
		prot1 <- as.character(df$protA[n])
	        prot2 <- as.character(df$protB[n])
	        p1 <- ggplotProteinScatterPlot(t(wtDat),prot1,prot2)
	        p2 <- ggplotProteinScatterPlot(t(koDat),prot1,prot2)
	        ggsave(paste0("M",i,"_","WT","_",gsub("\\|","_",prot1),"_",gsub("\\|","_",prot2),".tiff"),p1)
	        ggsave(paste0("M",i,"_","KO","_",gsub("\\|","_",prot1),"_",gsub("\\|","_",prot2),".tiff"),p2)
	}
}
}

# Generate corplots.
generate_corplots(wtModules)
generate_corplots(koModules)

# Any sig?
anySig <- rownames(statsdf)[statsdf$sigProt]
sum(prots %in% anySig)

# Sig enrichment???
obs = sum(prots %in% anySig)
exp = length(anySig)/3022 * length(prots)
obs/exp

# sig prots.
x = subset(statsdf,rownames(statsdf) %in% prots & statsdf$sigProt == TRUE)
head(x)

# Build  networks
generate_networks <- function(modules) {
	require(getPPIs)
	g <- list()
	for (i in 1:length(modules)) {
	     prots <- names(modules[[i]])
	     idx <- idy <- colnames(wtAdjm) %in% prots
	     subWT <- wtAdjm[idx,idy]
	     idx <- idy <- colnames(koAdjm) %in% prots
	     subKO <- koAdjm[idx,idy]
	     entrez <- protmap$entrez[match(colnames(subWT),protmap$prots)]
	     g[[i]] <- buildNetwork(musInteractome, mygenes=entrez, taxid=10090)
	}
	     return(g)
}

# Some modules are densely interconnected!
wtGraphs <- generate_networks(wtModules)
koGraphs <- generate_networks(koModules)
names(wtGraphs) <- names(wtModules)
names(koGraphs) <- names(koModules)

plot(koGraphs[[2]]$network)

#------------------------------------------------------------------------------
## GO Analysis.
#------------------------------------------------------------------------------

# Load previously compiled GO annotation collection:
musGOcollection <- readRDS(file.path(datadir,"musGOcollection.Rds"))

# Protein names (same for WT and KO).
prots <- colnames(wtAdjm)

# Function to perform GO analysis.
go_analysis <- function(partition){
	# import
	require(anRichment)
	# Get modules.
	modules <- split(partition,partition)
	# Build a matrix of labels.
	entrez <- protmap$entrez[match(names(partition),protmap$prots)]
	idx <- lapply(modules, function(x) names(partition) %in% names(x))
	labels_dm <- apply(as.matrix(do.call(cbind,idx)),2, function(x) as.numeric(x))
	# Perform GO Enrichment analysis with the anRichment library.
	GOenrichment <- enrichmentAnalysis(
						       classLabels = labels_dm,
						       identifiers = entrez,
						       refCollection = musGOcollection,
						       useBackground = "given",
						       threshold = 0.05,
						       thresholdType = "Bonferroni",
						       getOverlapEntrez = TRUE,
						       getOverlapSymbols = TRUE,
						       ignoreLabels = 0,
						       verbose = 1
						       )
	# Extract the results.
	GOdata <- lapply(GOenrichment$setResults, function(x) x[[2]])
	names(GOdata) <- paste0("M",names(modules))
	return(GOdata)
}

# Perform GO analysis for modules identified in WT and KO networks.
wtGO <- go_analysis(wtPartition)
koGO <- go_analysis(koPartition)

# Get top GO term associated with each module.
wtTopGO <- unlist(lapply(wtGO,function(x) x$shortDataSetName[1]))
koTopGO <- unlist(lapply(koGO,function(x) x$shortDataSetName[1]))

# Write GO results to file.
myfile <- file.path(tabsdir, "3_Supplementary_WT_Network_GO_Analysis.xlsx")
write.excel(wtGO, myfile)
myfile <- file.path(tabsdir, "3_Supplementary_KO_Network_GO_Analysis.xlsx")
write.excel(wtGO, myfile)

# ENDOFILE
#------------------------------------------------------------------------------

# Why were there more sig modules in WT before!!!?
# Under the null hypothesis that nothing is changing, modules with 
# average edge weight greater than or less than the null distribution are changing.

#-------------------------------------------------------------------------------
## Utilize permutation approach to identify divergent modules.
#-------------------------------------------------------------------------------



# Load expression data
wtDat <- t(readRDS(file.path(datadir,"wtDat.Rds")))
koDat <- t(readRDS(file.path(datadir,"koDat.Rds")))

# Fix rownames.
colnames(wtDat) <- colnames(koDat) <- rownames(readRDS(file.path(datadir,"wtDat.Rds")))

cleanDat <- rbind(wtDat,koDat)
colnames(cleanDat) <- colnames(wtDat)

cormat <- WGCNA::bicor(cleanDat)
rownames(cormat) <- colnames(cormat) <- colnames(cleanDat)

# RAise to power??
sft <-lapply(list(wt = wtDat,ko = koDat), function(x) {
		     WGCNA::pickSoftThreshold(x,
					powerVector=seq(4,20,by=1.0),
					corFnc="bicor",
					networkType="signed")})

plots <- lapply(sft,ggplotScaleFreeFit)

p1 <- plots$wt$ScaleFreeFit
p2 <- plots$ko$ScaleFreeFit

net <- wtAdjm - koAdjm

# Input for NetRep:
data_list        <- list(self = cleanDat)
correlation_list <- list(self = cormat) 
network_list     <- list(self = net) # THIS IS WHAT IS USED TO CALCULATE average edge weight!
module_list      <- list(wt = wtPartition, self = koPartition)

# Perform permutation testing.
preservation <- NetRep::modulePreservation(
					   network = network_list,
					   data = data_list,
					   correlation = correlation_list,
					   moduleAssignments = module_list,
					   modules = NULL,
					   backgroundLabel = 0,
					   discovery = "self",
					   test = "self",
					   selfPreservation = TRUE,
					   nThreads = 8,
					   #nPerm = 100000,  # determined by the function.
					   null = "overlap",
					   alternative = "two.sided", # c(greater,less,two.sided)
					   simplify = TRUE,
					   verbose = TRUE
					   )

q <- p.adjust(preservation$p.values[,1], method = "bonferroni")

#------------------------------------------------------------------------------
