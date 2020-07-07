#!/usr/bin/env Rscript

#' ---
#' title: 
#' description: 
#' authors: Tyler W A Bradshaw
#' ---

## OPTIONS:
BF_alpha = 0.05 # bonferonni significance threshold for module DA

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
tabsdir <- file.path(root, "tables")

# Load the data.
data(list=tolower(tissue)) # tidy_prot

# Load the gene map.
data(list=paste0(tolower(tissue),"_gene_map")) # gene_map

# Load the graph partition.
data(list=paste0(tolower(tissue),"_partition")) # partition

# Load the protein covariation network.
data(list=paste0(tolower(tissue),"_ne_adjm")) # ne_adjm

# Load the ppi network.
data(list=paste0(tolower(tissue),"_ppi_adjm")) # ppi_adjm

#---------------------------------------------------------------------
## Perform module level statistical analysis.
#---------------------------------------------------------------------
message("\nPerforming module-level analysis of differential abundance with EdgeR GLM.")


# Annotate data with module membership.
tidy_prot$Module <- partition[tidy_prot$Accession]

# Cast data into a dm, summarize a module as the sum of its proteins.
dm <- tidy_prot %>%
	group_by(Module, Genotype, Sample, Treatment) %>%
	dplyr::summarize(Sum.Intensity=sum(TAMPOR.Intensity),.groups="drop") %>% 
	as.data.table() %>%
	dcast(Module ~ Sample,value.var="Sum.Intensity") %>%
	as.matrix(rownames=TRUE)

# Create dge object.
dge <- DGEList(counts=dm)

# Sample to group mapping.
idx <- match(rownames(dge$samples),tidy_prot$Sample)
treatment <- tidy_prot$Treatment[idx]
genotype <- tidy_prot$Genotype[idx]
groups <- as.character(interaction(genotype,treatment))
groups[grepl("WT",groups)] <- "WT"
dge$samples$group <- groups

# Create design matrix.
design <- model.matrix( ~ 0 + group, data=dge$samples)

# Estimate dispersion.
dge <- estimateDisp(dge, design, robust=TRUE)

# Fit a model.
fit <- glmQLFit(dge, design, robust=TRUE)

# Create contrasts.
contr_list <- list(
		   "Shank2" = limma::makeContrasts('groupWT-groupShank2.KO',
						   levels=design),
		   "Shank3" = limma::makeContrasts('groupWT-groupShank3.KO',
						   levels=design),
		   "Syngap1" = limma::makeContrasts('groupWT-groupSyngap1.HET',
						    levels=design),
		   "Ube3a" = limma::makeContrasts('groupWT-groupUbe3a.KO',
						  levels=design)
		   )

# Assess differences.
qlf <- lapply(contr_list,function(x) glmQLFTest(fit,contrast=x))

# Collect results.
glm_results <- lapply(qlf, function(x) {
			      topTags(x,n=Inf,sort.by="p.value")$table
		   })
glm_results <- lapply(glm_results, as.data.table, keep.rownames="Module")

# Drop M0.
glm_results <- lapply(glm_results, function(x) x %>% filter(Module != 0))

# Adjust p-values for n module comparisons.
glm_results <- lapply(glm_results, function(x) {
			      x$PAdjust <- p.adjust(x$PValue,method="bonferroni")
			      return(x) })

# Calculate the number of proteins per module.
module_sizes <- sapply(split(partition,partition), length)

# Annotate results with the number of nodes per module:
glm_results <- lapply(glm_results, function(x) {
			      n <- module_sizes[as.character(x$Module)]
			      x <- tibble::add_column(x,Nodes=n,.after="Module")
			      })

# Annotate results with module proteins.
modules <- split(names(partition),partition)
module_ids <- lapply(modules,function(x) {
	       paste(gene_map$symbol[match(x,gene_map$uniprot)],x,sep="|")
			      })
glm_results <- lapply(glm_results, function(x){
			      x$Proteins <- sapply(module_ids[x$Module],
						   paste,collapse=";")
			      return(x)
			      })

#--------------------------------------------------------------------
## Calculate Module PVE
#--------------------------------------------------------------------
# NOTE: the data pre-TAMPOR was used to create the networks!

# Calculate module eigengenes using  WGCNA functions.
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
glm_results <- lapply(glm_results, function(x) {
			      PVE <- pve[paste0("M",x$Module)]
			      x <- tibble::add_column(x,PVE,.after="Nodes")
			      return(x)
				   })

# Summary of network:
df <- data.table("Tissue"= tissue,
		 "Number of Nodes" = formatC(length(partition),big.mark=","),
		 "Number of Modules" = length(module_sizes[-1]),
		 "Median Module Size" = median(module_sizes[-1]),
		 "Median Module PVE" = round(median(pve),3),
		 "Percent Clustered" = round(100*(sum(partition != 0)/length(partition)),3))
knitr::kable(df)

# Summary of sig modules.
lapply(names(glm_results),function(x) {
	       df <- glm_results[[x]] %>% filter(PAdjust < BF_alpha) %>% 
		       select(-one_of("Proteins"))
	       df <- tibble::add_column(df,"Genotype"=x,.before=1)
	       return(knitr::kable(df))
})

#--------------------------------------------------------------------
## Add module level data to glm results.
#--------------------------------------------------------------------

# Summarize a module as the sum of its proteins.
df <- tidy_prot %>%
	group_by(Module, Genotype, Sample, Treatment) %>%
	dplyr::summarize(Sum.Intensity=sum(TAMPOR.Intensity),.groups="drop") %>% 
	as.data.table() %>%
	dcast(Module ~ Sample,value.var="Sum.Intensity")
df$Module <- as.character(df$Module)

# Sort columns.
column_order <- c(1, # 'Module' 
		  # WT columns:
		  which(grepl("Shank2",colnames(df)) & grepl("WT",colnames(df))),
		  which(grepl("Shank3",colnames(df)) & grepl("WT",colnames(df))),
		  which(grepl("Syngap1",colnames(df)) & grepl("WT",colnames(df))),
		  which(grepl("Ube3a",colnames(df)) & grepl("WT",colnames(df))),
		  # KO columns
		  which(grepl("Shank2",colnames(df)) & grepl("KO|HET",colnames(df))),
		  which(grepl("Shank3",colnames(df)) & grepl("KO|HET",colnames(df))),
		  which(grepl("Syngap1",colnames(df)) & grepl("KO|HET",colnames(df))),
		  which(grepl("Ube3a",colnames(df)) & grepl("KO|HET",colnames(df)))
		  )
df <- df %>% select(all_of(column_order))
module_df <- df

# Add to the data.
for (i in c(1:length(glm_results))){
	glm_results[[i]] <- left_join(glm_results[[i]],module_df,by="Module")
}

#--------------------------------------------------------------------
## Calculate mean and SEM of groups.
#--------------------------------------------------------------------

#df <- glm_results
#cols <- c(grep("F[0-9]{1,2}_WT",colnames(df)),grep("MUT",colnames(df)))
#dm <- df %>% dplyr::select(Module,all_of(cols)) %>% 
#	as.data.table() %>% as.matrix(rownames="Module")
#idy <- grepl("WT",colnames(dm))
#WT_means <- apply(dm,1,function(x) log2(mean(x[idy])))
#WT_SEM <- apply(dm,1,function(x) sd(log2(x[idy])))/WT_means
#idy <- grepl("MUT",colnames(dm))
#MUT_means <- apply(dm,1,function(x) log2(mean(x[idy])))
#MUT_SEM <- apply(dm,1,function(x) sd(log2(x[idy])))/MUT_means
#df <- tibble::add_column(df,"WT Mean" = WT_means, .after="PVE")
#df <- tibble::add_column(df,"WT SEM" = WT_SEM, .after="WT Mean")
#df <- tibble::add_column(df,"MUT Mean" = MUT_means, .after="WT SEM")
#df <- tibble::add_column(df,"MUT SEM" = MUT_SEM, .after="MUT Mean")
#glm_results <- df

#--------------------------------------------------------------------
# Save results.
#--------------------------------------------------------------------

# Data.table describing graph partition.
Uniprot <- names(partition)
idx <- match(Uniprot,gene_map$uniprot)
Entrez <- gene_map$entrez[idx]
Symbol <- gene_map$symbol[idx]
part_dt <- data.table(Uniprot,Entrez,Symbol,Module=partition)

# Save as excel table.
myfile <- file.path(tabsdir,paste0(tissue,"_Module_GLM_Results.xlsx"))
results <- c(list("Network Partition" = part_dt), glm_results)
write_excel(results,file=myfile)
