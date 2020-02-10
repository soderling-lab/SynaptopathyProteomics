#!/usr/bin/env Rscript

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

## User parameters to change:
pfile = "2020-02-10_Cortex_Surprise_Module_Self_Preservation.RData"
#pfile = "2020-02-10_Striatum_Surprise_Module_Self_Preservation.RData"
net = "Cortex" # Which network are we analyzing? 
do_DBD_enrichment = TRUE
do_GO_enrichment = FALSE
do_module_analysis = FALSE

# Global options and imports.
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
})

# Directories.
if (rstudioapi::isAvailable()) {
	setwd("D:/projects/SynaptopathyProteomics/analysis/3_Network-Analysis")
}
here <- getwd()
root <- dirname(dirname(here))
funcdir <- file.path(root, "R")
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
netsdir <- file.path(root, "networks")
figsdir <- file.path(root,"figs", Sys.Date())

# Create directory for figure output.
if (dir.exists(figsdir)) { 
	message(paste("Warning, overwriting files in:\n",figsdir))
} else { 
	dir.create(figsdir,recursive = TRUE) 
}

# Functions.
source_funcdir <- function(funcdir){
  myfun <- list.files(funcdir, pattern="*.R",full.names = TRUE)
  invisible(sapply(myfun, source))
}
source_funcdir(funcdir)

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
sampleTraits <- readRDS(file.path(rdatdir, "2_Combined_traits.RData"))

# Load co-expression (adjacency) matrices.
myfiles <- c(
  Cortex = file.path(rdatdir, "3_Cortex_Adjm.RData"),
  Striatum = file.path(rdatdir, "3_Striatum_Adjm.RData")
)
adjm <- as.matrix(readRDS(myfiles[net]))
rownames(adjm) <- colnames(adjm)

# Load network partitions-- self-preservation enforced.
myfile <- file.path(rdatdir,pfile)
partition <- readRDS(myfile)[[1]]

# Reset index.
partition <- reset_index(partition)

# Load theme for plots.
ggtheme()

#---------------------------------------------------------------------
## Module enrichment for DBD-associated genes.
#---------------------------------------------------------------------

# Load Disease ontology.
geneSet <- "mouse_Combined_SFARI_geneSets.RData"
#geneSet <- "mouse_Combined_DBD_geneSets.RData"
myfile <- file.path(rdatdir,geneSet)
DBDcollection <- readRDS(myfile)

# Perform disease enrichment analysis.
DBDresults <- moduleGOenrichment(list(partition),1, protmap, DBDcollection)
names(DBDresults) <- gsub("-",".",names(DBDresults))

# Save.
#myfile <- file.path(rdatdir,paste0("3_",net, "_", ptype, "_Module_DBD_Enrichment.RData"))
#saveRDS(DBDresults,myfile)

# Collect modules with significant enrichment of DBD-genes.
method <- "Bonferroni" 
alpha <- 0.05
is_sig <- sapply(DBDresults,function(x) any(x[[method]] < alpha))
DBDsig <- names(DBDresults)[is_sig]
nDBDsig <- length(DBDsig)

# Status.
message(paste("Total number of disease associated modules:",nDBDsig))


# Get Modules.
modules <- split(partition, partition)
names(modules) <- paste0("M", names(modules))

# Number of modules.
nModules <- sum(names(modules) != "M0")
message(paste("Number of modules:", nModules))

# Module size statistics.
mod_stats <- summary(sapply(modules, length)[!names(modules) == "M0"])[-c(2, 5)]
message(paste("Minumum module size:",mod_stats["Min."]))
message(paste("Median module size:",mod_stats["Median"]))
message(paste("Maximum module size:",mod_stats["Max."]))

# Percent not clustered.
percentNC <- sum(partition == 0) / length(partition)
message(paste("Percent of proteins not clustered:", 
	      round(100 * percentNC, 2), "(%)"))

# Calculate Module Eigengenes.
# Note: Soft power does not influence MEs.
MEdata <- moduleEigengenes(data,
  colors = partition, # Do not need to sort in same order!
  softPower = 1, impute = FALSE
)
MEs <- as.matrix(MEdata$eigengenes)

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
medianPVE <- median(PVE[names(PVE) != "M0"])
message(paste("Median module coherence (PVE):", 
	      round(100 * medianPVE, 2), "(%)."))

# Create list of MEs.
ME_list <- lapply(seq(ncol(MEs)),function(x) MEs[,x]) # Do it this way to preserve names.
names(ME_list) <- names(modules) # Same as colnames MEs.

# Remove M0. Do this before p-value adjustment.
ME_list <- ME_list[which(names(ME_list)!="M0")]

# Sample to group mapping.
sampleTraits$Sample.Model.Tissue <- paste(sampleTraits$Sample.Model, sampleTraits$Tissue, sep = ".")
groups <- sampleTraits$Sample.Model.Tissue[match(rownames(MEs), sampleTraits$SampleID)]
names(groups) <- rownames(MEs)

# Group all WT samples from a tissue type together.
groups[grepl("WT.*.Cortex", groups)] <- "WT.Cortex"
groups[grepl("WT.*.Striatum", groups)] <- "WT.Striatum"

# Fix levels (order).
group_order <- c("WT","KO.Shank2","KO.Shank3", "HET.Syngap1","KO.Ube3a")
group_levels <- paste(group_order,net,sep=".")

# Perform Kruskal Wallis tests to identify modules whose summary
# expression profile is changing.
KWdata_list <- lapply(ME_list, function(x) { kruskal.test(x ~ groups[names(x)]) })
KWdata <- as.data.frame(do.call(rbind,KWdata_list))[-c(4,5)]

# Correct p-values for n comparisons.
method <- "bonferroni"
KWdata$p.adj <- p.adjust(as.numeric(KWdata$p.value), method)

# Significant modules.
alpha <- 0.05
sigModules <- rownames(KWdata)[KWdata$p.adj < alpha]
nSigModules <- length(sigModules)
message(paste0(
  "Number of modules with significant (p.adj < ", alpha, ")",
  " Kruskal-Wallis test: ", nSigModules,"."
))

# Dunnetts test for post-hoc comparisons.
# Note: P-values returned by DunnettTest have already been adjusted for 
# multiple comparisons!
control_group <- paste("WT", net, sep = ".")
DTdata_list <- lapply(ME_list, function(x) {
  g <- factor(groups[names(x)],levels=group_levels)
	df <- as.data.frame({
			    DunnettTest(x ~ g,control = control_group)[[control_group]] 
			  })
	return(df)
})

# Number of significant changes.
alpha <- 0.05
nSigDT <- sapply(DTdata_list, function(x) sum(x$pval < alpha))
if (length(nSigDT[sigModules]) > 0) {
  message("Summary of Dunnett's test changes for significant modules:")
  print(nSigDT[sigModules])
}

# Numer of significant modules with disease association.
idx <- which(paste0("R",r,".",names(nSigDT[sigModules])) %in% DBDsig) 
diseaseSig <- names(nSigDT[sigModules])[idx]
nSigDisease <- length(diseaseSig)
if (nSigDisease > 0) {
  message(paste("Number of significant modules with",
		"significant enrichment of DBD-associated genes:", nSigDisease))
  message("Summary of Dunnett's test changes for significant, DBD-associated modules:")
  print(nSigDT[diseaseSig])
}

# Generate boxplots.
bplots <- lapply(ME_list,function(x) {
		 ggplotVerboseBoxplot(x,groups,group_levels)
  })
names(bplots) <- names(ME_list)
# Add R#.M# + PVE + pvalue to plot titles. Simplify x-axis labels.
x_labels <- rep(c("WT","Shank2 KO","Shank3 KO",
		  "Syngap1 HET","Ube3a KO"),2)
# Loop to clean-up plots.
for (k in seq_along(bplots)) {
	# Add title and fix xlabels.
	plot <- bplots[[k]]
	m <- names(bplots)[k]
	namen <- paste0("R",r,".",m)
	txt <- paste0("P.adj = ", round(KWdata[m,"p.adj"],3),
		      "; ","PVE = ", round(PVE[m], 3))
	plot_title <- paste0(namen, " (", txt, ")")
	plot$labels$title <- plot_title
	plot <- plot + scale_x_discrete(labels = x_labels)
	# Add significance stars!
	df <- data.table(xpos=c(2:5),
			 ypos = 1.01 * max(plot$data$x),
			 p=DTdata_list[[m]]$pval,
			 symbol="")
	df$symbol[df$p<0.05] <- "*"
	df$symbol[df$p<0.005] <- "**"
	df$symbol[df$p<0.0005] <- "***"
	if (any(df$p<0.05)) {
		plot <- plot + annotate("text",x=df$xpos,y=df$ypos,label=df$symbol,size=7) }
	# Store results in list.
	bplots[[k]] <- plot
} # Ends loop to fix plots.

