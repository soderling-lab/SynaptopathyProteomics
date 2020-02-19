#!/usr/bin/env Rscript

#' ---
#' title: Network Analysis
#' description: Protein co-expression network analysis
#' authors: Tyler W Bradshaw
#' ---

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

## User parameters to change:
data_type <- "Combined" # Cortex, Striatum, or Combined...
part_type <- "Striatum" # Specify part type when working with comb data.

# Data files.
input_files <- list(adjm_files = list(Cortex="3_Cortex_Adjm.RData",
				      Striatum="3_Striatum_Adjm.RData",
				      Combined="3_Combined_Adjm.RData"),
		    data_files = list(Cortex="3_Cortex_cleanDat.RData",
				      Striatum="3_Striatum_cleanDat.RData",
				      Combined="3_Combined_cleanDat.RData"),
		    part_files = list(Cortex=list(self="2020-02-10_Cortex_Surprise_Module_Self_Preservation.RData",
						  ppi ="2020-02-13_Cortex_PPI_Module_Self_Preservation.RData",
						  other="2020-02-18_Cortex_Striatum_Module_Self_Preservation.RData"),
				      Striatum=list(self="2020-02-10_Striatum_Surprise_Module_Self_Preservation.RData",
						    ppi = "2020-02-13_Striatum_PPI_Module_Self_Preservation.RData",
						    other="2020-02-19_Striatum_Cortex_Module_Self_Preservation.RData"))
		    )

# Global imports.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(WGCNA)
  library(org.Mm.eg.db)
  library(anRichment)
  library(getPPIs)
  library(DescTools)
  library(igraph)
  library(ggplot2)
  library(enrichR)
  library(gtable)
  library(cowplot)
})

# Directories.
here <- getwd()
subdir <- basename(here)
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
netsdir <- file.path(root, "networks",data_type)
figsdir <- file.path(root, "figs",subdir,data_type)
tabsdir <- file.path(root, "tables", subdir,data_type)

# Functions.
suppressWarnings({ devtools::load_all() })

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load GLM stats.
myfile <- file.path(rdatdir, "2_GLM_Stats.RData")
glm_stats <- readRDS(myfile)

# Get proteins with any significant change.
idy <- lapply(c("Cortex","Striatum"),function(x) grep(x,colnames(glm_stats$FDR)))
sigProts <- lapply(idy, function(x) {
			   apply(glm_stats$FDR[,x],1,function(pval) any(pval<0.05)) })
names(sigProts) <- c("Cortex","Striatum")
sigProts[["Combined"]] <- apply(cbind(sigProts[[1]],sigProts[[2]]),1,any)

# Load expression data.
# Data should be transposed: rows, proteins.
# Insure that column names are not lost after transpose.
myfile <- file.path(rdatdir,input_files$data_file[[data_type]])
temp <- readRDS(myfile)
data <- t(temp)
colnames(data) <- rownames(temp)

# Load Sample info.
sampleTraits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load co-expression (adjacency) matrix.
myfile <- file.path(rdatdir,input_files$adjm_file[[data_type]])
adjm <- as.matrix(readRDS(myfile))
rownames(adjm) <- colnames(adjm)

# Load PPI adjacency matrix.
adjm_ppi <- fread(file.path(rdatdir,"3_PPI_Adjm.csv"),drop=1)
adjm_ppi <- as.matrix(adjm_ppi)
rownames(adjm_ppi) <- colnames(adjm_ppi)

# Load network partitions--self-preservation enforced.
myfiles <- file.path(rdatdir,input_files$part_file[[part_type]])
partitions <- lapply(myfiles,function(x) unlist(readRDS(x)))
names(partitions) <- names(input_files$part_file[[part_type]])

# Load theme for plots.
ggtheme()

# Reset partition index for self-preserved modules.
partitions$self <- reset_index(partitions$self)

#---------------------------------------------------------------------
## Check which modules are preserved in opposite tissue.
#---------------------------------------------------------------------

# Load partitions.
p1 <- partitions$self
p2 <- partitions$other

# Get preserved modules.
preserved_modules <- list()
modules <- split(p2,p1)
names(modules) <- paste0("M",names(modules))
preserved <- which(sapply(modules,function(x) unique(x)!=0))
preserved_modules[["other"]] <- names(modules)[preserved]

# Fraction of modules that are preserved in the other tissue type.
nModules <- length(unique(p1[p1!=0]))
nPres <- length(preserved_modules[["other"]])
pPres <- round(100*(nPres/length(unique(p1))),3)
message(paste0(nPres," of ",nModules," (",pPres,"%) ", part_type, 
	       " modules are preserved in the opposite tissue."))

# Reformat other partition, so that module names are the same 
# as in self-preserved partition (partitions$self)!
m <- split(p1,p2)
new_part <- sapply(seq_along(m),function(x) {
			   p <- rep(as.numeric(names(m)[x]),length(m[[x]]))
			   names(p) <- names(m[[x]])
			   return(p)
			   })
partitions$other <- unlist(new_part)

#---------------------------------------------------------------------
## Which modules are preserved in PPI graph?
#---------------------------------------------------------------------

# Load partitions.
p1 <- partitions$self
p2 <- partitions$ppi

# Get preserved modules.
modules <- split(p2,p1)
preserved <- which(sapply(modules,function(x) unique(x)!=0))
preserved_modules[["ppi"]] <- paste0("M",names(modules)[preserved])

# Fraction of modules that are preserved in the other tissue type.
nModules <- length(unique(p1[p1!=0]))
nPres <- length(preserved_modules[["ppi"]])
pPres <- round(100*(nPres/length(unique(p1))),3)
message(paste0(nPres," of ",nModules," (",pPres,"%) ", part_type, 
	       " modules are preserved in the PPI network."))

#---------------------------------------------------------------------
## SigProt annotations--which genotype is a given protein changing in?
#---------------------------------------------------------------------

# Create a df of protein-SigProt annotations.
df <- glm_stats$FDR
colnames(df) <- gsub(" FDR","",colnames(df))
df$sigProt <- apply(df,1, function(x){
			   paste(colnames(df)[x<0.05],collapse="; ")
})
df$sigProt[df$sigProt == ""] <- NA
sigProtAnno <- df$sigProt
names(sigProtAnno) <- rownames(df)

#---------------------------------------------------------------------
## Collect all modules in a list.
#---------------------------------------------------------------------

# Create list of modules.
module_list <- list()
partition <- partitions$self

# Module list with entrez ids.
idx <- match(names(partition),protmap$ids)
module_list[["Entrez"]] <- split(protmap$entrez[idx],partition)

# Module list with gene symbols.
module_list[["Symbols"]] <- split(protmap$gene[idx],partition)

# Module list with protein ids.
module_list[["IDs"]] <- split(partition,partition)

# Name modules with M prefix.
module_list <- lapply(module_list,function(x) {
			  names(x) <- paste0("M",names(x))
			  return(x)
})

# Drop M0 from lists.
module_list <- lapply(module_list,function(x) x[-which(names(x)=="M0")])

#---------------------------------------------------------------------
## Module enrichment for DBD-associated genes.
#---------------------------------------------------------------------

# Load Disease ontology.
DBDset <- "mouse_Combined_DBD_collection.RData"
myfile <- file.path(rdatdir,DBDset)
DBDcollection <- readRDS(myfile)

# Perform disease enrichment analysis.
gene_list <- module_list$Entrez
DBDenrichment <- gse(gene_list, DBDcollection)

# Collect modules with significant enrichment of DBD-genes.
method <- "FDR"
alpha <- 0.05
any_sig <- which(sapply(DBDenrichment,function(x) any(x[method] < alpha)))
DBDsig <- names(DBDenrichment)[any_sig]

# Status.
message(paste("Total number of DBD-associated modules:",
	      length(DBDsig)))

# Write to file.
myfile <- file.path(tabsdir,"3_Module_DBD_Enrichment.xlsx")
write_excel(DBDenrichment[DBDsig],myfile)

# Collect all DBD genes.
DBDgenes <- lapply(DBDcollection$dataSets,function(x) {
			   as.character(x$data$Entrez) 
})
names(DBDgenes) <- sapply(DBDcollection$dataSets,function(x) x$name)

# DBDprots.
DBDprots <- lapply(DBDgenes,function(x) {
			  protmap$ids[which(x %in% protmap$entrez)]
})

# Create a vector of protein-DBD annotations.
DBDcols <- do.call(cbind,lapply(DBDprots,function(x) protmap$ids %in% x))
colnames(DBDcols) <- names(DBDprots)
DBDdf <- as.data.frame(DBDcols)
#DBDdf <- tibble::add_column(DBDdf,"Protein"=protmap$ids,.before=1)
#DBDdf$anyDBD <- apply(DBDcols,1,any)
rownames(DBDdf) <- protmap$ids
DBDanno <- apply(DBDdf,1,function(x) paste(colnames(DBDdf)[x],collapse="; "))
DBDanno[DBDanno == ""] <- NA

#---------------------------------------------------------------------
## Create a table summarizing disease genes.
#---------------------------------------------------------------------

# Summarize all DBD genes identified in synaptosome.

# Summary of disease genes.
df <- t(as.data.frame(sapply(DBDprots,length)))
rownames(df) <- NULL

# Modify default table theme to change font size.
# Cex is a scaling factor relative to the defaults.
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params = list(cex = 0.75)),
  colhead = list(fg_params = list(cex = 0.75)),
  rowhead = list(fg_params = list(cex = 0.75))
)
# Create table and add borders.
# Border around data rows.
mytable <- tableGrob(df, rows = NULL, theme = mytheme)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 1, b = nrow(mytable), l = 1, r = ncol(mytable)
)
mytable <- gtable_add_grob(mytable,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
  t = 2, l = 1, r = ncol(mytable)
)

# Check the table.
fig <- plot_grid(mytable)

# Save.
myfile <- prefix_file(file.path(figsdir,"DBD_Gene_Summary.tiff"))
ggsaveTable(mytable,myfile)

#--------------------------------------------------------------------
## Module GO enrichment using anRichment.
#--------------------------------------------------------------------

# Build a GO collection.
GOcollection <- buildGOcollection(organism="mouse")

# Perform gene set enrichment analysis.
GOresults <- gse(gene_list,GOcollection)

# Top (1) go term for every module.
method <- "Bonferroni"
alpha <- 0.05
topGO <- lapply(GOresults,function(x) {
		       	p <- x[[method]][1] 
			names(p) <- x$shortDataSetName[1]
			return(p)
	    })

# Number of modules with any significant GO term enrichment.
message(paste("Total number of modules with any significant GO",
	      "enrichment:", sum(unlist(topGO)<alpha)))


# Modules with signifcant GO enrichment:
# Order by significance.
GOsig <- topGO[topGO < alpha]
GOsig <- GOsig[order(unlist(GOsig))]

# Write GO results to file.
myfile <- file.path(tabsdir,"Module_GO_Enrichment.xlsx")
write_excel(GOresults,myfile)

#--------------------------------------------------------------------
## Module enrichment using enrichR.
#--------------------------------------------------------------------

# EnrichR is an online platform for gene set enrichment.
# https://amp.pharm.mssm.edu/Enrichr/

# EnrichR can perform enrichment analysis querying a large number of 
# gene collections: https://amp.pharm.mssm.edu/Enrichr/#stats

# Collect list of module genes--input is gene symbols.
gene_list <- module_list$Symbols

# All available databases.
dbs <- listEnrichrDbs()
dbs <- dbs$libraryName[order(dbs$libraryName)]

# Enrichment analysis for OMIM disorders.
# FIXME: progress report would be nice.
#db <- "GO_Molecular_Function_2018"
db <- "OMIM_Disease"
OMIMresults <- enrichR(gene_list,db)
OMIMresults <- unlist(OMIMresults,recursive=FALSE) # Why nested list?
names(OMIMresults) <- sapply(strsplit(names(OMIMresults),"\\."),"[",1)

# Remove NA.
out <- which(sapply(OMIMresults,function(x) dim(x)[1]==0))
OMIMresults <- OMIMresults[-out]

# Which modules are enriched for OMIM disorders?
# Top OMIM term for every module.
alpha <- 0.05
topOMIM <- lapply(OMIMresults,function(x) {
		       	p <- x$Adjusted.P.value[1] 
			names(p) <- x$Term[1]
			return(p)
			})


# Collect any signficant terms.
OMIMsig <- topOMIM[which(topOMIM < alpha)]

# Save results to file...
myfile <- file.path(tabsdir,"Module_OMIM_Enrichment.xlsx")
write_excel(OMIMresults,myfile)

#---------------------------------------------------------------------
## Explore changes in module summary expression.
#---------------------------------------------------------------------

# Get partition of interest.
# If analyzing the combined data, focus on modules that are preserved
# in opposite tissue type.
if (data_type == "Combined") {
	partition <- partitions$other
} else {
	# Otherwise focus on modules that are self-preserved.
	partition <- partitions$self
}

# Modules from the given partition.
modules <- split(partition,partition)
names(modules) <- paste0("M",names(modules))
modules <- modules[-which(names(modules)=="M0")]

# Total Number of modules.
all_modules <- module_list$IDs
nModules <- sum(names(all_modules) != "M0")
message(paste("Number of modules:", nModules))

# Module size statistics.
mod_stats <- summary(sapply(modules, length))[-c(2, 5)]
message(paste("Minumum module size:",mod_stats["Min."]))
message(paste("Median module size:",mod_stats["Median"]))
message(paste("Maximum module size:",mod_stats["Max."]))

# Percent not clustered.
percentNC <- sum(partition == 0) / length(partition)
message(paste("Percent of proteins not clustered:", 
	      round(100 * percentNC, 2), "(%)"))

# Calculate Module Eigengenes.
# Note: Soft power does not influence MEs.
# Note: Do not need to sort partition to be in the same order!
MEdata <- moduleEigengenes(data,
  colors = partition, 
  excludeGrey = TRUE, # Ignore M0!
  softPower = 1, 
  impute = FALSE
)
MEs <- as.matrix(MEdata$eigengenes)

# Create list of MEs.
# Do it this way to preserve names.
ME_list <- lapply(seq(ncol(MEs)),function(x) MEs[,x]) 
names(ME_list) <- names(modules)

# Module membership (KME).
KMEdata <- signedKME(data,MEdata$eigengenes,corFnc="bicor",
		     outputColumnName = "M")
KME_list <- lapply(seq(ncol(KMEdata)),function(x) {
			   v <- vector("numeric",length=nrow(KMEdata))
			   names(v) <- rownames(KMEdata)
			   v[] <- KMEdata[[x]]
			   v <- v[order(v,decreasing=TRUE)]
			   return(v)
		     })
names(KME_list) <- colnames(KMEdata)

# Get Percent Variance explained (PVE).
PVE <- as.numeric(MEdata$varExplained)
names(PVE) <- names(modules)
medianPVE <- median(PVE)
message(paste("Median module coherence (PVE):", 
	      round(100 * medianPVE, 2), "(%)."))

# Sample to group mapping for statistical testing.
sampleTraits$Sample.Model.Tissue <- paste(sampleTraits$Sample.Model, 
					  sampleTraits$Tissue, sep = ".")
idx <- match(rownames(MEs), sampleTraits$SampleID)
groups <- sampleTraits$Sample.Model.Tissue[idx]
names(groups) <- rownames(MEs)

# Group all WT samples from a tissue type together.
groups[grepl("WT.*.Cortex", groups)] <- "WT.Cortex"
groups[grepl("WT.*.Striatum", groups)] <- "WT.Striatum"

# If grouping all WT samples together...
if (data_type == "Combined") {
	# WT from both tissue will be annotated as WT.
	groups[grepl("WT", groups)] <- "WT"
}

# Perform Kruskal Wallis tests to identify modules whose summary
# expression profile is changing.
KWdata_list <- lapply(ME_list, function(x) { 
			      kruskal.test(x ~ groups[names(x)]) 
					  })
KWdata <- as.data.frame(do.call(rbind,KWdata_list))[-c(4,5)]

# Correct p-values for n comparisons.
KWdata$p.adj <- as.numeric(KWdata$p.value) * nModules
KWdata$p.adj[KWdata$p.adj > 1.0] <- 1.0

# Significant modules.
alpha <- 0.1
sigModules <- rownames(KWdata)[KWdata$p.adj < alpha]
nSigModules <- length(sigModules)
message(paste0(
  "Number of modules with significant (p.adj < ", alpha, ")",
  " Kruskal-Wallis test: ", nSigModules,"."
))

# Perform Dunnetts test for post-hoc comparisons.
# Note: P-values returned by DunnettTest have already been adjusted for 
# multiple comparisons!

# Define control group and levels (order) for DunnettTest.
group_order <- c("WT","KO.Shank2","KO.Shank3", "HET.Syngap1","KO.Ube3a")
group_levels <- c(paste(group_order,"Cortex",sep="."),
		  paste(group_order,"Striatum",sep="."))

# If combining data, set all WT samples to WT.
if (data_type == "Combined") {
	control_group <- paste("WT", sep = ".")
	group_levels[grepl("WT", group_levels)] <- "WT"
	group_levels <- unique(group_levels)
} else {
	control_group <- paste("WT", part_type, sep = ".")
}

# Loop to perform DTest. NOTE: This takes several seconds.
DTdata_list <- lapply(ME_list, function(x) {
  g <- factor(groups[names(x)],levels=group_levels)
  result <- DunnettTest(x ~ g,control = control_group)[[control_group]] 
  return(as.data.frame(result))
})

# Number of modules with significant KW + DT changes.
DTalpha <- 0.05
nSigDT <- sapply(DTdata_list, function(x) sum(x$pval < DTalpha))
if (length(nSigDT[sigModules]) > 0) {
  message("Summary of Dunnett's test changes for significant modules:")
  print(nSigDT[sigModules])
}

# Numer of significant modules with disease association.
sigDBDmodules <- nSigDT[sigModules[which(sigModules %in% DBDsig)]]
nSigDisease <- length(sigDBDmodules)
if (nSigDisease > 0) {
  message(paste("Number of significant modules with",
		"significant enrichment of DBD-associated genes:", nSigDisease))
  message("Summary of Dunnett's test changes for significant, DBD-associated modules:")
  print(sigDBDmodules)
}

#---------------------------------------------------------------------
## Generate verbose boxplots.
#---------------------------------------------------------------------

# Reset sample to group mapping for generating boxplots.
sampleTraits$Sample.Model.Tissue <- paste(sampleTraits$Sample.Model, 
					  sampleTraits$Tissue, sep = ".")
idx <- match(rownames(MEs), sampleTraits$SampleID)
groups <- sampleTraits$Sample.Model.Tissue[idx]
names(groups) <- rownames(MEs)

# Group all WT samples from a tissue type together.
groups[grepl("WT.*.Cortex", groups)] <- "WT.Cortex"
groups[grepl("WT.*.Striatum", groups)] <- "WT.Striatum"

# Reset group order.
group_order <- c("WT","KO.Shank2","KO.Shank3", "HET.Syngap1","KO.Ube3a")
group_levels <- c(paste(group_order,"Cortex",sep="."),
		  paste(group_order,"Striatum",sep="."))

# Generate boxplots summarizing module protein expression.
plots <- lapply(ME_list,function(x) {
		 ggplotVerboseBoxplot(x,groups,
				      group_levels,control_group="WT")
  })
names(plots) <- names(ME_list)

## Loop to  clean-up plots.
# Simplify x-axis labels.
x_labels <- rep(c("WT","Shank2 KO","Shank3 KO",
		  "Syngap1 HET","Ube3a KO"),2)
# Loop to clean-up plots.
for (k in seq_along(plots)) {
	# Add title and fix xlabels.
	plot <- plots[[k]]
	m <- names(plots)[k]
	# Add significance stars!
	DTdata_list[[m]]$pval
	dt <- DTdata_list[[m]]
	df <- data.table(plot$data)
	idx <- match(paste(df$g,control_group,sep="-"),rownames(dt))
	df$pval <- dt$pval[idx]
	df$pval[is.na(df$pval)] <- 1
	df$xpos <- df$g
	df$ypos <- 1.1 * max(df$x)
	df$symbol <- ""
	df$symbol[df$pval<0.05] <- "*"
	df$symbol[df$pval<0.005] <- "**"
	df$symbol[df$pval<0.0005] <- "***"
	df$tissue <- factor(df$tissue,levels=c("Cortex","Striatum"))
	df$color <- "black"
	df$color[df$pval<0.05] <- "red"
	# If KW Sig, then add stars.
	if (KWdata[m,"p.adj"] < 0.1){
		plot <- plot + geom_text(data=df,aes(x=xpos, y=ypos,
						     label=symbol,size=7))
	}
	# Fix title and xlabels.
	txt <- paste0("P.adj = ", round(KWdata[m,"p.adj"],3),
		      "; ","PVE = ", round(PVE[m], 3))
	plot_title <- paste0(m, " (", txt, ")")
	plot$labels$title <- plot_title
	plot <- plot + scale_x_discrete(labels = x_labels)
	# Store results in list.
	plots[[k]] <- plot
} # Ends loop to fix plots.

# Create table summarizing partition statistics.
df <- data.table("Algorithm" = "Leiden",
		 "Quality Function" = "Surprise",
		 "N Modules"=nModules,
		 "Percent Not-clustered"=round(100*percentNC,3),
		 "Median PVE" = round(medianPVE,3),
		 "N Sig. Modules"=nSigModules,
		 "Min Size" = mod_stats[[1]],
		 "Median Size" = mod_stats[[2]],
		 "Max Size" = mod_stats[[4]])
mytable <- tableGrob(df, rows=NULL, theme = mytheme)

# Check the table.
fig <- plot_grid(mytable)

# Save.
myfile <- prefix_file(file.path(figsdir,"Module_Summary.tiff"))
ggsaveTable(mytable,myfile)

#---------------------------------------------------------------------
## Save verbose box plots.
#---------------------------------------------------------------------

# Save sig plots as single pdf.
myfile <- prefix_file(file.path(figsdir,"Sig_Module_Boxplots.pdf"))
ggsavePDF(plots[sigModules],myfile)

#---------------------------------------------------------------------
## Save sigProt boxplots for sig modules.
#---------------------------------------------------------------------

# Load plots.
myfiles <- c(Cortex=file.path(rdatdir,"All_Cortex_SigProt_Boxplots.RData"),
	     Striatum=file.path(rdatdir,"All_Striatum_SigProt_Boxplots.RData"),
	     Combined=file.path(rdatdir,"All_Faceted_SigProt_Boxplots.RData"))
myfile <- myfiles[data_type]
all_plots <- readRDS(myfile)

# Group by module.
plot_list <- split(all_plots,partition[names(all_plots)])
names(plot_list) <- paste0("M",names(plot_list))

# Save a single pdf containing all the sign proteins within a
# module for each module.
## FIXME: suppress output from grid.arrange!
for (i in 1:length(sigModules)){
	module_name <- sigModules[i]
	plots <- plot_list[[module_name]]
	groups <- rep(c(1:ceiling(length(plots)/4)),each=4)[c(1:length(plots))]
	plot_groups <- split(plots,groups)
	figs <- lapply(plot_groups,function(x) {
			       fig <- gridExtra::grid.arrange(grobs=x,ncol=2,nrow=2)
			       return(fig)
		 })
	myfile <- prefix_file(file.path(figsdir,
				paste0(module_name,"_Module_SigProts.pdf")))
	ggsavePDF(figs,myfile)
}

#---------------------------------------------------------------------
## Examine overall structure of network.
#---------------------------------------------------------------------

# Examine relationships between modules by comparing their summary
# expression profiles (MEs).

# Calculate correlations between ME vectors.
adjm_me <- cor(do.call(cbind,ME_list))

## Perform Heirarchical clustering of the ME adjm.
# Convert similarity matrix to distance matrix, and cluster with hclust.
hc <- hclust(as.dist(1 - adjm_me), method = "ward.D2")

# Optimize cut of tree with modularity.
# First remove any negative edges from graph.
g <- graph_from_adjacency_matrix(adjm_me,mode="undirected",weighted=TRUE)
g <- delete_edges(g,E(g)[E(g)$weight < 0])

# Generate a bunch of partitions by cutting the tree at various heights.
h <- seq(0,max(hc$height),by=0.005)
hc_partitions <- lapply(h,function(x) cutree(hc,h=x))
# Number of groups, k.
k <- sapply(hc_partitions, function(x) length(unique(x)))
# Modularity, q.
q <- sapply(hc_partitions,function(x) {
  modularity(g, x, weights = abs(edge_attr(g, "weight")))})

# Find the best cut height--the cut height that maximizes modularity.
best_q <- unique(q[seq(h)[q==max(q)]])
best_h <- median(h[seq(h)[q==max(q)]])
best_k <- unique(k[seq(k)[q==max(q)]])
message(paste0("Cut height that produces the best partition: ",best_h,"."))
message(paste0("Number of groups: ",best_k," (Modularity = ",round(best_q,3),")."))

## With the optimal k, cut the graph.
hc_partition <- cutree(hc, k=best_k)
groups <- split(hc_partition,hc_partition)

# Get representative module from each group, its medoid.
# The medoid is the module which is closes (i.e. most similar) 
# to all others in its group.
rep_modules <- getMedoid(adjm_me,h=best_k)

# Update dendro with cutheight and representative modules.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE, labels = FALSE) + 
  geom_hline(yintercept=best_h, color='red', size = 1)
dend_data <- ggdendro::dendro_data(as.dendrogram(hc))
dend_data <- dend_data$labels
dend_data$group <- as.factor(hc_partition[dend_data$label])
dend_data$rep_module <- dend_data$label %in% rep_modules
dendro <- dendro + 
  geom_text(data = dend_data, aes(x, y, label = label, color = rep_module),
            hjust = 1, angle = 90, size = 3) + 
  scale_colour_manual(values=c("black", "red")) +
  theme(legend.position="none")

#---------------------------------------------------------------------
## Generate module colors based on their dist to rep modules.
#---------------------------------------------------------------------

# Assign modules a color based on similarity with three rep modules.
df <- do.call(cbind,lapply(rep_modules,function(x) adjm_me[,x]))
colnames(df) <- paste0(rep_modules,"cor")

# Rescale the data to [0,1].
dm <- (df - min(df))/(max(df) - min(df))
colnames(dm) <- c("R","G","B")
df <- data.table(cbind(df,dm))
rownames(df) <- colnames(adjm_me)

# Convert RGB to hexadecimal color.
df$color <-  rgb(255*df$R, 255*df$G, 255*df$B, maxColorValue=255)

# Collect color assignments.
module_colors <- df$color
names(module_colors) <- names(modules)

# Add colored bars to dendro.
# FIXME: need to improve this.
#dend_data <- ggdendro::dendro_data(as.dendrogram(hc))
#dend_data <- dend_data$labels
#df <- data.frame("Color" = module_colors[as.character(dend_data$label)],
#		 "Module" = as.character(dend_data$label))

#p2<-ggplot(df,aes(x=Module,y=1,fill=Color)) + geom_tile() +
#	scale_fill_manual(values=as.character(df$Color)) + 
#	theme(legend.position="none")

#maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
#gp1$widths[2:5] <- as.list(maxWidth)
#gp2$widths[2:5] <- as.list(maxWidth)


# Save.
myfile <- prefix_file(file.path(figsdir,"Modules_Dendro.tiff"))
ggsave(myfile,plot=dendro, height=3, width = 3)

#---------------------------------------------------------------------
## Generate ppi graphs and co-expression graphs.
#---------------------------------------------------------------------

## NOTE: Coerce boolean attributes to integer to avoid warnings when
# loading into cytoscape.

# Create co-expression graph.
exp_graph <- graph_from_adjacency_matrix(adjm,mode="undirected",
				  weighted=TRUE, diag=FALSE)

# Create PPI graph.
ppi_graph <- graph_from_adjacency_matrix(adjm_ppi,mode="undirected",
				  weighted=TRUE, diag=FALSE)

# Remove NAs from PPI edges.
E(ppi_graph)$weight[which(is.na(E(ppi_graph)$weight))] <- 0

## Add attributes to igraph object.
# Add Gene symbols.
symbols <- protmap$gene[match(names(V(exp_graph)),protmap$ids)]
exp_graph <- set_vertex_attr(exp_graph,"symbol",value = symbols)

# Add sigProt vertex attribute.
anySig <- as.numeric(sigProts[[data_type]][names(V(exp_graph))])
exp_graph <- set_vertex_attr(exp_graph, "sigProt", 
			 value = anySig)

# Add any DBDprot vertex attribute.
DBDnodes <- lapply(DBDprots,function(x) names(V(exp_graph)) %in% x)
for (DBD in names(DBDnodes)){
	exp_graph <- set_vertex_attr(exp_graph, name=DBD, 
				 value = as.numeric(DBDnodes[[DBD]]))
}

# Collect PPI evidence.
myfile <- file.path(rdatdir,"3_All_PPIs.RData")
ppis <- readRDS(myfile)

# Map mouse entrez to protein ids.
ppis$ProteinA <- protmap$ids[match(ppis$osEntrezA,protmap$entrez)]
ppis$ProteinB <- protmap$ids[match(ppis$osEntrezB,protmap$entrez)]
out <- is.na(ppis$ProteinA) | is.na(ppis$ProteinB)
ppis <- ppis[!out,]

# Get the relevant columns.
ppis <- ppis %>% select(ProteinA,ProteinB,osEntrezA,osEntrezB,
			Interactor_A_Taxonomy,Interactor_B_Taxonomy,
			Source_database,Confidence_score,
			Publications,Methods)

# Save to file..
myfile <- file.path(tabsdir,paste0("3_All_PPIs.csv"))
fwrite(ppis,myfile)

#---------------------------------------------------------------------
## Generate cytoscape graphs.
#---------------------------------------------------------------------

# Create graphs.
for (i in c(1:length(modules))) {
	message(paste("Working on module",i,"..."))
	module_name = names(modules)[i]
	nodes = names(modules[[module_name]])
	module_kme = KME_list[[module_name]]
	network_layout = 'force-directed edgeAttribute=weight'
	image_file = file.path(figsdir,"NETWORKS",module_name)
	image_format = "SVG"
	createCytoscapeGraph(exp_graph,ppi_graph,nodes,
			     module_kme,module_name,
			     module_colors, network_layout,
			     output_file, image_file,
			     image_format)
	# When done, save cytoscape session.
	if (i == length(modules)) {
		if (data_type == "Combined") {
			myfile <- file.path(netsdir,paste0(data_type,"_",part_type,".cys"))
		} else {
			myfile <- file.path(netsdir,paste0(data_type,".cys"))
		}
		winfile <- gsub("/mnt/d/","D:/",myfile)
		saveSession(winfile)
	}
}

#---------------------------------------------------------------------
## Save key results summarizing modules.
#---------------------------------------------------------------------

# Summarize modules: Name, N Nodes, PVE, Color, Prots.
module_summary <- data.frame(Module = names(PVE),
			     Nodes = sapply(modules,length),
			     PVE = PVE,
			     Color = module_colors)

# Combine with KWdata.
tempKW <- KWdata
colnames(tempKW) <- paste("KW",colnames(tempKW))
module_summary <- cbind(module_summary,tempKW)

# Remove parameter column.
module_summary$"KW parameter" <- NULL

# Add column for KW sig.
module_summary$"KW sig" <- module_summary$"KW p.adj" < 0.05

# Combine with DT results.
reformatDT <- function(x){
	df <- as.data.table(x,keep.rownames=TRUE) %>% select(rn,diff,pval)
	df <- melt(df, id.vars="rn")
	values <- df$value
	names(values) <- paste("DT",sapply(strsplit(df$rn,"-"),"[",1),df$variable) 
	return(values)
}
dm <- do.call(rbind,lapply(DTdata_list,reformatDT))
module_summary <- cbind(module_summary,dm)

# Number of sig changes.
module_summary$"N DT sig" <- nSigDT

# Any DT sig.
module_summary$"Any DT sig" <- module_summary$"N DT sig" > 0

## More detailed summary of every module.
dfs <- lapply(seq_along(KME_list), function(x) {
	       df <- data.table(Protein = names(partition))
	       df$Module <- partition
	       df$KME <- KME_list[[x]][df$Protein]
	       df$sigProt <- df$Protein %in% sigProts
	       df$"Genotypes" <- sigProtAnno[df$Protein]
	       df$DBDProt <- df$Protein %in% DBDprots
	       df$"DBD Association(s)" <- DBDanno[df$Protein]
	       df <- df %>% filter(Module == x)
	       df <- df[order(df$KME,decreasing=TRUE),]
		 })
names(dfs) <- names(modules)

# Add expression data.
for (i in seq_along(dfs)){
	x <- dfs[[i]]
	y <- do.call(cbind,glm_stats[c(1,2,4,5)])
	colnames(y) <- sapply(strsplit(colnames(y),"\\."),"[",2)
	y$Protein <- rownames(y)
	y <- y %>% filter(Protein %in% x$Protein)
	z <- merge(x,y,by="Protein")
	dfs[[i]] <- z
}

# Write to file.
results <- list()
results[["Summary"]] <- module_summary
results = c(results,dfs[sigModules])
myfile <- file.path(tabsdir,"Module_Summary.xlsx")
write_excel(results,myfile)
