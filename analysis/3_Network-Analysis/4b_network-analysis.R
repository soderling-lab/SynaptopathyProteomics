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
data_type <- "Striatum" # Cortex, Striatum, or Combined...
part_type <- "Striatum" # Specify part type when working with comb data.
generate_cytoscape_graphs <- FALSE

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
  library(gtable)
  library(cowplot)
  library(RCy3)
})

# Directories.
here <- getwd()
subdir <- basename(here)
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
netsdir <- file.path(root, "networks",part_type)
figsdir <- file.path(root, "figs",subdir,data_type,part_type)
tabsdir <- file.path(root, "tables", subdir,data_type)

# Remove any existing figures and tables.
invisible(sapply(list.files(figsdir),unlink))
invisible(sapply(list.files(tabsdir),unlink))

# Functions.
suppressWarnings({ devtools::load_all() })

# Load protein identifier map.
protmap <- readRDS(file.path(rdatdir, "2_Protein_ID_Map.RData"))

# Load GLM stats.
myfile <- file.path(rdatdir, "2_GLM_Stats.RData")
glm_stats <- readRDS(myfile)

# Get proteins with any significant change.
idy <- lapply(c("Cortex","Striatum"),function(x) {
		      grep(x,colnames(glm_stats$FDR)) })
sigProts <- lapply(idy, function(x) {
			   apply(glm_stats$FDR[,x],1,function(pval) {
					any(pval<0.05) }) 
		      })
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
## Which modules are preserved in opposite tissue?
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
## Which modules are preserved in the PPI graph?
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
## Which modules are enriched for DBD-associated genes?
#---------------------------------------------------------------------

# Load Disease ontology.
DBDset <- "mouse_Combined_DBD_collection.RData"
DBDset <- "2020-02-21_mouse_Combined_DBD_collection.RData"
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

# Add DBD genes and terms to list of preserved modules.
# Function to reformat pvalues and combine with term nam.e
fx <- function(x) { return(formatC(x,digits=2,format="e")) }
term_pval <- sapply(DBDenrichment[DBDsig], function(x) {
			    paste0(x$shortDataSetName[x$FDR<0.05],
				   " (FDR = ",fx(x$FDR[x$FDR<0.05]),")") }) 
preserved_modules$dbd <- term_pval 

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
rownames(DBDdf) <- protmap$ids
DBDanno <- apply(DBDdf,1,function(x) paste(colnames(DBDdf)[x],collapse="; "))
DBDanno[DBDanno == ""] <- NA

#---------------------------------------------------------------------
## Create a table summarizing disease genes.
#---------------------------------------------------------------------

# Summarize all DBD genes identified in synaptosome.
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
## Which modules are enriched for GO terms?
#--------------------------------------------------------------------

## Use the anRichment package.
# Build a GO collection.
message("Performing GO analysis with anRichment...")
GOcollection <- buildGOcollection(organism="mouse")

# Perform gene set enrichment analysis.
GOresults <- gse(gene_list,GOcollection)

# Significant go terms for every module.
method <- "Bonferroni"
alpha <- 0.05
topGO <- lapply(GOresults,function(x) {
			idx <- x[[method]] < alpha
			p <- x[[method]][idx]
			names(p) <- x$shortDataSetName[idx]
			return(p)
	    })

# Number of modules with any significant GO term enrichment.
sigGO <- names(which(sapply(topGO,function(x) length(x) > 0)))
message(paste("Total number of modules with any significant GO",
	      "enrichment:", length(sigGO)))

# Add to list of preserved modules.
term_pval <- sapply(topGO,function(x) {
		      paste0(names(x)," (p.adj = ",fx(x),")") })
preserved_modules$go <- term_pval

#--------------------------------------------------------------------
## Which modules are enriched for ASD DEGs?
#--------------------------------------------------------------------

# Load ASD DEG collection.
GeneSet <- "mouse_ASD_DEG_collection.RData"
myfile <- file.path(rdatdir,GeneSet)
ASDcollection <- readRDS(myfile)

# All ASD DEGs
ASD_DEGs <- unique(ASDcollection$dataSets[[1]]$data$Entrez)
idx <- match(ASD_DEGs[which(ASD_DEGs %in% protmap$entrez)],protmap$entrez)
ASDprots <- protmap$ids[idx]

# Perform gene set enrichment analysis.
gene_list <- module_list$Entrez
ASDresults <- gse(gene_list,ASDcollection)

# Significant go terms for every module.
method <- "Bonferroni"
alpha <- 0.05
topASD <- lapply(ASDresults,function(x) {
			idx <- x[[method]] < alpha & x$enrichmentRatio > 1
			p <- x[[method]][idx]
			names(p) <- x$shortDataSetName[idx]
			return(p)
	    })

# Number of modules with any significant GO term enrichment.
sigASD <- names(which(sapply(topASD,function(x) length(x) > 0)))
message(paste("Total number of modules with significant enrichment",
	      "of ASD DEGs:",length(sigASD)))

# Add to list of preserved modules.
term_pval <- sapply(topASD[sigASD],function(x) {
		      paste0(names(x)," (p.adj = ",fx(x),")") })
preserved_modules$asd <- term_pval

#--------------------------------------------------------------------
## Module enrichment using enrichR.
#--------------------------------------------------------------------

# EnrichR is an online platform for gene set enrichment.
# https://amp.pharm.mssm.edu/Enrichr/

# EnrichR can perform enrichment analysis querying a large number of 
# gene collections: https://amp.pharm.mssm.edu/Enrichr/#stats
message("Performing enrichment analysis with enrichR...")

# Collect list of module genes--input is gene symbols.
gene_list <- module_list$Symbols

# Enrichment analysis for OMIM disorders.
db <- "OMIM_Disease"
OMIMresults <- enrichR(gene_list,db,quiet=TRUE)

# Remove NA.
out <- which(sapply(OMIMresults,function(x) dim(x)[1]==0))
OMIMresults <- OMIMresults[-out]

# Check, all enriched?
if (sum(sapply(OMIMresults,function(x) any(x$Odds.Ratio<1)))>0){
	message("Warning, some terms are depleted.")
}

# Which modules are enriched for OMIM disorders?
# Top OMIM term for every module.
alpha <- 0.05
topOMIM <- lapply(OMIMresults,function(x) {
			idx <- x$Adjusted.P.value < alpha
			p <- x$Adjusted.P.value[idx]
			names(p) <- x$Term[idx]
			return(p)
	    })

# Number of modules with any significant GO term enrichment.
sigOMIM <- names(which(sapply(topOMIM,function(x) length(x) > 0)))
message(paste("Total number of modules with any significant OMIM",
	      "enrichment:", length(sigOMIM)))

# Add to list of preserved modules.
term_pval <- sapply(topOMIM[sigOMIM],function(x) {
		      paste0(names(x)," (p.adj = ",fx(x),")") })
preserved_modules$omim <- term_pval

#---------------------------------------------------------------------
## Other EnrichR enrichment.
#---------------------------------------------------------------------

# Enrichment analysis for PFAM domains.
message("Performing enrichment analysis with enrichR...")
db <- "Pfam_Domains_2019"
PFAMresults <- enrichR(gene_list,db,quiet=TRUE)

# Remove NA.
out <- which(sapply(PFAMresults,function(x) dim(x)[1]==0))
PFAMresults <- PFAMresults[-out]

# Sig PFAM terms from every module.
alpha <- 0.05
topPFAM <- lapply(PFAMresults,function(x) {
			idx <- x$Adjusted.P.value < alpha & x$Odds.Ratio > 1
			p <- x$Adjusted.P.value[idx]
			names(p) <- x$Term[idx]
			return(p)
	    })

# Number of modules with any significant enrichment.
sigPFAM <- names(which(sapply(topPFAM,function(x) length(x) > 0)))
message(paste("Total number of modules with any significant PFAM",
	      "domain enrichment:", length(sigPFAM)))

# Add to list of preserved modules.
term_pval <- sapply(topPFAM[sigPFAM],function(x) {
		      paste0(names(x)," (p.adj = ",fx(x),")") })
preserved_modules$pfam <- term_pval

#---------------------------------------------------------------------
## Explore changes in module summary expression.
#---------------------------------------------------------------------

# Total Number of modules.
# P.values will be corrected for n comparisions.
all_modules <- module_list$IDs
nModules <- sum(names(all_modules) != "M0")
message(paste("Number of modules:", nModules))

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

# If not working with combined data,
# drop modules that are preserved in other tissue.
# These will be analyzed seperately.
if (data_type != "Combined") {
	out <- which(names(modules) %in% preserved_modules$other)
	modules <- modules[-out]
}

# Remove M0.
modules <- modules[-which(names(modules)=="M0")]

# Fix partition--only analyzing modules defined above.
out <- partition %notin% as.numeric(gsub("M","",names(modules)))
partition[out] <- 0

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

# Numer of significant modules with any DBD association.
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
if (data_type != "Combined") { control_group <- NULL }
plots <- lapply(ME_list,function(x) {
		 ggplotVerboseBoxplot(x,groups,
				      group_levels,control_group)
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

#---------------------------------------------------------------------
## Save verbose box plots.
#---------------------------------------------------------------------

# Save all modules.
myfile <- prefix_file(file.path(figsdir,"All_Module_Boxplots.pdf"))
ggsavePDF(plots,myfile)

# Save sig modules as single pdf.
myfile <- prefix_file(file.path(figsdir,"Sig_Module_Boxplots.pdf"))
ggsavePDF(plots[sigModules],myfile)

#---------------------------------------------------------------------
## Save sigProt boxplots for sig modules.
#---------------------------------------------------------------------

# Load plots.
# If Cortex or striatum -- then only cortex or striatum are plotted.
# If Combined -- then data from both tissues are plotted.
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
	n <- length(plots)
	groups <- rep(c(1:ceiling(n/4)),each=4)[c(1:n)]
	plot_groups <- split(plots,groups)
	figs <- lapply(plot_groups,function(x) {
			       fig <- grid.arrange(grobs=x,ncol=2,nrow=2)
			       return(fig)
		 })
	# Save.
	myfile <- prefix_file(file.path(figsdir,
				paste0(module_name,"_Module_SigProts.pdf")))
	ggsavePDF(figs,myfile)
	# Close the device.
	if (i == length(sigModules)) { dev.off() }
}

#---------------------------------------------------------------------
## Examine overall structure of network.
#---------------------------------------------------------------------

# Examine relationships between modules by comparing their summary
# expression profiles (MEs).

# Load the data and partition.
# Don't use the combined dataset.
myfile <- file.path(rdatdir,input_files$data[[part_type]])
tempdat <- readRDS(myfile)
tempdat <- t(tempdat)
part <- partitions$self

# Calculate module eigengenes.
MEdata <- moduleEigengenes(tempdat,
  colors = part, 
  excludeGrey = TRUE, # Ignore M0!
  softPower = 1, 
  impute = FALSE
)
MEs <- as.matrix(MEdata$eigengenes)

# List of MEs.
ME_list <- lapply(seq(ncol(MEs)),function(x) MEs[,x]) 
names(ME_list) <- names(all_modules)

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

# Get module which is most different!
#rep_modules <- sapply(groups,function(x) {
#	       colSums <- apply(adjm_me[names(x),names(x)],2,sum)
#	       rep_module <- names(which(colSums == max(colSums)))
#	       return(rep_module)
#  })

# Get dendrogram data, update with group and rep_module.
dend_data <- ggdendro::dendro_data(as.dendrogram(hc))
dend_data <- dend_data$labels
dend_data$group <- hc_partition[as.character(dend_data$label)]
dend_data$rep_module <- dend_data$label %in% rep_modules

# Add actual position to dend data.
dend_data$x <- match(as.character(dend_data$label),hc$label[hc$order])
dend_data <- dend_data %>% arrange(x)
dend_data$label <- factor(dend_data$label,levels=dend_data$label)

# Generate dendrogram.
dendro <- ggdendro::ggdendrogram(hc, rotate = FALSE, labels = TRUE) + 
  geom_hline(yintercept=best_h, color='red', size = 1)
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
df <- tibble::add_column(df,Module=colnames(adjm_me),.before=1)
df <- df %>% arrange(match(df$Module,as.character(dend_data$label)))

# Convert RGB to hexadecimal color.
df$color <-  rgb(255*df$R, 255*df$G, 255*df$B, maxColorValue=255)

# Collect color assignments.
module_colors <- df$color
names(module_colors) <- names(all_modules)

# Generate colored bars.
dend_data$color <- module_colors[as.character(dend_data$label)]
dend_data$color <- factor(dend_data$color,levels=dend_data$color)
p2 <- ggplot(dend_data,aes(x=label,y=1,fill=color)) + geom_tile() +
	scale_fill_manual(values=as.character(dend_data$color)) + 
	theme(legend.position="none", axis.text.x = element_text(angle=45))

# Save.
myfile <- prefix_file(file.path(figsdir,"Modules_Dendro.tiff"))
ggsave(myfile,plot=dendro, height=3, width = 3)

# Save.
myfile <- prefix_file(file.path(figsdir,"Module_Colors.tiff"))
ggsave(myfile,plot=p2, height=3, width = 3)

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
write_excel(list(PPIs=ppis),myfile)

#---------------------------------------------------------------------
## Generate cytoscape graph summarizing overall topolgy of the network.
#---------------------------------------------------------------------

#---------------------------------------------------------------------
## Generate cytoscape graphs of all modules.
#---------------------------------------------------------------------

# If working with Combined data, append graphs to tissue specific 
# Cytoscape file.
if (generate_cytoscape_graphs) {
# Prompt the user to open Cytoscape if it is not open.
	cytoscape_ping()
	if (data_type == "Combined") {
		cysfile <- file.path(netsdir,paste0(part_type,".cys"))
		if (file.exists(cysfile)){
			message(paste("Adding graphs to",part_type,"file!"))
			winfile <- gsub("/mnt/d/","D:/",cysfile)
			openSession(winfile)
		} else {
			message(paste("Analyze",part_type,"data first.",
			              "Combined graphs will be appended to",
			              "this file."))
	}
}

# Create graphs.
	for (i in c(1:length(modules))) {
		module_name = names(modules)[i]
		message(paste("Working on module", module_name,"..."))
		nodes = names(modules[[module_name]])
		module_kme = KME_list[[module_name]]
		network_layout = 'force-directed edgeAttribute=weight'
		image_file = file.path(dirname(figsdir),"Networks",module_name)
		image_format = "SVG"
		createCytoscapeGraph(exp_graph,ppi_graph,nodes,
			     module_kme,module_name,
			     module_colors, network_layout,
			     output_file, image_file,
			     image_format)
		# When done, save cytoscape session.
		if (i == length(modules)) {
			myfile <- file.path(netsdir,paste0(part_type,".cys"))
			winfile <- gsub("/mnt/d/","D:/",myfile)
			saveSession(winfile)
		}
	}
} # ENDS IF CHUNK

#---------------------------------------------------------------------
## Save key results summarizing modules.
#---------------------------------------------------------------------

# Summarize modules: Name, N Nodes, PVE, Color, Prots.
module_summary <- data.frame(Module = names(PVE),
			     Nodes = sapply(modules,length),
			     PVE = PVE)

# Combine with KWdata.
tempKW <- KWdata
colnames(tempKW) <- paste("KW",colnames(tempKW))
module_summary <- cbind(module_summary,tempKW)

# Remove parameter column.
module_summary$"KW parameter" <- NULL

# Add column for KW sig.
module_summary$"KW sig" <- module_summary$"KW p.adj" < 0.1

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

# Other protein attributes: PPI pres, tissue pres, OMIM, DBD, GOsig.
module_summary$PPI <- as.character(module_summary$Module) %in% preserved_modules$ppi
module_summary$DBD <- sapply(preserved_modules$dbd,function(x) paste(x,collapse=" | "))[as.character(module_summary$Module)]
module_summary$OMIM <- sapply(preserved_modules$omim,function(x) paste(x,collapse=" | "))[as.character(module_summary$Module)]
module_summary$PFAM <- sapply(preserved_modules$pfam,function(x) paste(x,collapse=" | "))[as.character(module_summary$Module)]
module_summary$ASD_DEGs <- sapply(preserved_modules$asd,function(x) paste(x,collapse=" | "))[as.character(module_summary$Module)]
module_summary$"Preserved in opposite tissue" <- sapply(preserved_modules$other,function(x) paste(x,collapse=" | "))[as.character(module_summary$Module)]

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
results <- c(results,dfs[sigModules])
# Insure that if working with combined data we know which
# data and partition type...
if (data_type == "Combined") {
	myfile <- file.path(tabsdir,
			    paste0("3_",data_type,"_",part_type,
				   "_Module_Summary.xlsx"))
} else {
	myfile <- file.path(tabsdir,
			    paste0("3_",data_type,"_Module_Summary.xlsx"))
}
write_excel(results,myfile)

# Remove that pesky Rplots.
unlink("Rplots.pdf")
