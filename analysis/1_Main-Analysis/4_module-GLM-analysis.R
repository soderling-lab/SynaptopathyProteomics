#!/usr/bin/env Rscript

#' ---
#' title: 
#' description: 
#' authors: Tyler W A Bradshaw
#' ---

## OPTIONS:
BF_alpha = 0.05

#---------------------------------------------------------------------
## Misc function - getrd().
#---------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
	in_root <- function(h=here, dir=dpat) { 
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE))) 
		return(check)
	}
	# Loop to find root.
	while (!in_root(here)) { 
		here <- dirname(here) 
	}
	root <- here
	return(root)
}

# Parse the command line arguments.
parse_args <- function(default="Cortex", args=commandArgs(trailingOnly=TRUE)){
	# Input must be Cortex or Striatum.
	msg <- c("Please specify a tissue type to be analyzed:\n",
	 "Choose either 'Cortex' or 'Striatum'.")
	# If interactive, return default tissue.
	if (interactive()) { 
		return("Cortex") 
	} else {
		# Check arguments.
		check <- !is.na(match(args[1], c("Cortex", "Striatum")))
		if (length(args == 1) & check) { 
			tissue  <- args[1]
			start <- Sys.time()
			message(paste("Starting analysis at:", start))
			message(paste0("Analyzing ", tissue,"..."))
		} else {
			stop(msg) 
		}
		return(tissue)
	}
}

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. Load custom functions and prepare
# the project directory for saving output files.

# Parse input.
tissue <- parse_args()

# Load renv -- use renv::load NOT activate!
root <- getrd()
renv::load(root,quiet=TRUE)

# Load required packages and functions.
suppressPackageStartupMessages({
	library(dplyr) # For manipulating data.
	library(edgeR) # For statistical analysis.
	library(data.table) # For working with tables.
})

# Load additional functions.
suppressWarnings({ devtools::load_all() })

# Project directories:
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Load the data.
data(list=tolower(tissue)) # tidy_prot

# Load the gene map.
data(gene_map) # gene_map

# Load the graph partition.
data(list=paste0(tolower(tissue),"_partition")) # partition

# Load the network.
data(list=paste0(tolower(tissue),"_ne_adjm")) # ne_adjm

#---------------------------------------------------------------------
## Perform module level statistical analysis.
#---------------------------------------------------------------------

# Calculate the number of proteins per module.
module_sizes <- sapply(split(partition,partition), length)

# Annotate data with module membership.
tidy_prot$Module <- partition[tidy_prot$Accession]

# Cast data into a dm, summarize a module as the sum of its proteins.
dm <- tidy_prot %>% 
	group_by(Module, Genotype, Sample, Treatment) %>%
	dplyr::summarize(Sum.Intensity=sum(Intensity),.groups="drop") %>% 
	as.data.table() %>%
	dcast(Module ~ Sample,value.var="Sum.Intensity") %>%
	as.matrix(rownames=TRUE)

# Create dge object.
dge <- DGEList(counts=dm)

# Sample to group mapping.
groups <- rownames(dge$samples)
idx <- match(groups,tidy_prot$Sample)
treatment <- tidy_prot$Treatment[idx]
treatment[grepl("KO|HET",treatment)] <- "MUT"
dge$samples$Batch <- tidy_prot$PrepDate[idx]
dge$samples$Genotype <- tidy_prot$Genotype[idx]
dge$samples$Treatment <- treatment

# Create design matrix.
design <- model.matrix( ~ Batch + Genotype + Treatment, data=dge$samples)

# Create contrast.
#contr <- limma::makeContrasts("GenotypeKO - GenotypeWT",levels=design)

# Estimate dispersion.
dge <- estimateDisp(dge, design, robust=TRUE)

# Fit a model.
fit <- glmQLFit(dge, design, robust=TRUE)

# Default comparison is last contrast: genotype.
#qlf <- glmQLFTest(fit,contrast=contr)
qlf <- glmQLFTest(fit)

# Collect results.
glm_results <- topTags(qlf,n=Inf,sort.by="p.value")$table %>%
	as.data.table(keep.rownames="Module")

# Drop M0.
glm_results <- glm_results %>% filter(Module != 0)

# Adjust p-values for n module comparisons.
glm_results$PAdjust <- p.adjust(glm_results$PValue,method="bonferroni")

# Fix logCPM column -- > convert to percentWT.
idy <- which(colnames(glm_results)=="logCPM")
colnames(glm_results)[idy] <- "PercentWT"
glm_results$PercentWT <- 2^glm_results$logFC

# Number of nodes per module:
n <- module_sizes[as.character(glm_results$Module)]
glm_results <- tibble::add_column(glm_results,Nodes=n,.after="Module")

# Number of significant modules.
message(paste("\nTotal number of modules:",length(unique(partition)) -1))
nsig <- sum(glm_results$PAdjust < BF_alpha)
message(paste0("\nNumber of significant ",
	      "(p-adjust < ", BF_alpha,") ",
	      "modules: ", nsig,"."))

# Pretty print sig results:
glm_results %>% filter(PAdjust < BF_alpha) %>%
	knitr::kable()

# Save sig modules.
myfile <- file.path(root,"data","sig_modules.rda")
sig_modules <- paste0("M",glm_results$Module[glm_results$PAdjust < BF_alpha])
save(sig_modules,file=myfile,version=2)

#--------------------------------------------------------------------
## Calculate Module PVE
#--------------------------------------------------------------------

# Calculate module eigengenes.
dm <- tidy_prot %>% filter(!grepl("QC",Sample)) %>% as.data.table() %>% 
	dcast(Sample ~ Accession, value.var = "Intensity") %>% 
	as.matrix(rownames="Sample") %>% log2()

ME_data <- WGCNA::moduleEigengenes(dm, colors = partition, 
				   excludeGrey = TRUE, softPower = 1 ,
				   impute = FALSE)

# Extract PVE.
pve <- as.numeric(ME_data$varExplained)
names(pve) <- gsub("X","M",names(ME_data$varExplained))

# Add to glm_results.
glm_results <- tibble::add_column(glm_results,
				  PVE=pve[paste0("M",glm_results$Module)],
				  .after="Nodes")

# Status.
message(paste("\nMedian module PVE:",
	      round(100*median(glm_results$PVE),3),"%."))

#--------------------------------------------------------------------
# Determine module hubs.
#--------------------------------------------------------------------

# All modules.
module_list <- split(partition,partition)
names(module_list) <- paste0("M",names(module_list))

# Get top three proteins.
# Sorted by node weighted degree.
module_hubs <- lapply(module_list, function(x) {
	       prots <- names(x)
	       subadjm <- ne_adjm[prots,prots]
	       node_degree <- apply(subadjm,2,sum)
	       node_degree <- node_degree[order(node_degree,decreasing=TRUE)]
	       hubs <- names(head(node_degree,3))
	       symbols <- gene_map$gene[match(hubs,gene_map$uniprot)]
	       ids <- paste(symbols,hubs,sep="|")
	       return(ids)
				  } )

# We will add hubs to the data in the code below.
				  
#--------------------------------------------------------------------
## Add module level data.
#--------------------------------------------------------------------

df <- tmt_protein %>% group_by(Module,Genotype,Treatment) %>% 
	dplyr::summarize(Intensity = sum(Adjusted.Intensity),.groups="drop")
df$.groups <- NULL
dm <- df %>% as.data.table() %>%
	dcast(Module ~ Fraction + Genotype, value.var = "Intensity")
dm$Module <- as.character(dm$Module)

# Sort columns.
idy <- grepl("F[0-9]{1,2}",colnames(dm))
lvls <- c("F4","F5","F6","F7","F8","F9","F10")
f <- factor(sapply(strsplit(colnames(dm)[idy],"_"),"[",1),levels=lvls)
g <- factor(sapply(strsplit(colnames(dm)[idy],"_"),"[",2),levels=c("WT","MUT"))
idy <- c("Module",gsub("\\.","_",as.character(levels(interaction(f,g)))))
dm <- dm %>% select(all_of(idy))

# Add to the data.
glm_results <- left_join(glm_results,dm,by="Module")

# Annotate with hubs.
hubs_list <- module_hubs[paste0("M",glm_results$Module)]
glm_results$Hubs <- sapply(hubs_list,paste,collapse="; ")

# Annotate with module proteins.

#--------------------------------------------------------------------
## Calculate mean and SEM of groups.
#--------------------------------------------------------------------

df <- glm_results
cols <- c(grep("F[0-9]{1,2}_WT",colnames(df)),grep("MUT",colnames(df)))
dm <- df %>% dplyr::select(Module,all_of(cols)) %>% 
	as.data.table() %>% as.matrix(rownames="Module")
idy <- grepl("WT",colnames(dm))
WT_means <- apply(dm,1,function(x) log2(mean(x[idy])))
WT_SEM <- apply(dm,1,function(x) sd(log2(x[idy])))/WT_means
idy <- grepl("MUT",colnames(dm))
MUT_means <- apply(dm,1,function(x) log2(mean(x[idy])))
MUT_SEM <- apply(dm,1,function(x) sd(log2(x[idy])))/MUT_means
df <- tibble::add_column(df,"WT Mean" = WT_means, .after="PVE")
df <- tibble::add_column(df,"WT SEM" = WT_SEM, .after="WT Mean")
df <- tibble::add_column(df,"MUT Mean" = MUT_means, .after="WT SEM")
df <- tibble::add_column(df,"MUT SEM" = MUT_SEM, .after="MUT Mean")
glm_results <- df

#--------------------------------------------------------------------
# Save results.
#--------------------------------------------------------------------

# Save sig proteins.
sig_proteins <- unique(tmt_protein$Accession[tmt_protein$PAdjust < 0.05])
myfile <- file.path(root,"data","sig_proteins.rda")
save(sig_proteins,file=myfile,version=2)
print(length(sig_proteins))

# Data.table describing graph partition.
Uniprot <- names(partition)
idx <- match(Uniprot,gene_map$uniprot)
Entrez <- gene_map$entrez[idx]
Symbol <- gene_map$gene[idx]
part_dt <- data.table(Uniprot,Entrez,Symbol,Module=partition)

# Save as excel table.
myfile <- file.path(tabsdir,"S3_Swip_TMT_Module_GLM_Results.xlsx")
results <- list("Network Partition" = part_dt,
		"Module GLM Results" = glm_results)
write_excel(results,file=myfile)

# Save as rda object.
module_stats <- glm_results
myfile <- file.path(datadir,"module_stats.rda")
save(module_stats,file=myfile,version=2)
