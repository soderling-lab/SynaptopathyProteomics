#!/usr/bin/env Rscript

# title: Analysis of Synaptopathy TMT Proteomics
# description: preprocessing and statistical analysis of TMT proteomics
# author: Tyler W Bradshaw

## ---- INPUTs

# project root
root <- "~/projects/SynaptopathyProteomics"


## ---- functions 

fx2str <- function(fx) {
	# return a formula (function, fx) as a single character string 
	fx_char <- as.character(fx)
	fx_str <- paste(fx_char[2],fx_char[1],fx_char[3])
	return(fx_str)
} #EOF


testProtContrast <- function(tidy_prot, fx, LT) {
  # performs statistical comparisons for a given model and contrast, LT
  # all proteins
  proteins <- unique(tidy_prot$Accession)
  # loop to fit linear model to each protein
  fit_list <- list()
  pbar <- txtProgressBar(max=length(proteins),style=3)
  for (prot in proteins) {
  	fm <- lm(fx, data = tidy_prot %>% filter(!grepl("QC",Treatment)) %>% subset(Accession==prot))
  	fit_list[[prot]] <- fm
  	setTxtProgressBar(pbar,value=match(prot,proteins))
  }
  close(pbar)
  # get s2 and df
  s2 <- sapply(fit_list, function(x) sigma(x)^2)
  df <- sapply(fit_list, function(x) x$df.residual)
  # moderated stats
  eb_var <- limma::squeezeVar(s2, df)
  df_prior <- eb_var$df.prior
  s2_prior <- eb_var$s2.prior
  if (is.null(s2_prior)) { s2_prior <- 0 ; warning("S2 prior = 0")}
  # perform tests with moderation
  results <- lapply(fit_list, function(x) lmTestContrast(x, LT, s2_prior, df_prior)) %>% 
  	bind_rows(.id="Protein") %>% 
  	mutate(FDR = p.adjust(Pvalue,method="fdr")) %>% 
  	arrange(Pvalue)
  # return the results
  return(results)
} #EOF


annotProts <- function(df, gene_map) {
  # annot with gene symbols
  idx <- match(df$Protein, gene_map$uniprot)
  df <- tibble::add_column(df, Symbol=gene_map$symbol[idx], .after="Protein")
  df <- tibble::add_column(df, Entrez=gene_map$entrez[idx], .after="Symbol")
  return(df)
} # EOF


## ---- R environment 

# Prepare the R workspace for the analysis. Load custom functions and prepare
# the project directory for saving output files.

# load required packages
suppressPackageStartupMessages({
  library(dplyr) # for working with data
  library(data.table) # for working with tables
  library(tidyProt) # for lmTestContrast
})

# load project specific data and functions
devtools::load_all(root)

# parse input
tissue <- parseArgs()


# the TMT data in root/data:
data(gene_map)
data(list=tolower(tissue)) # e.g. data(cortex); tidy_prot


## ---- get uniprot IDs for genes of interest

# saved in root/data
data(shank2)
data(shank3)
data(syngap1)
data(ube3a)


## ---- fit individual protein models for each Genotype and mutant Gene

# fit the same model for each experiment (mixture)
fx <- log2(Intensity) ~ 0 + Condition

message("fit: ", fx2str(fx))

# fit a model: Shank2
df_shank2 <- tidy_prot %>% 
	filter(Genotype == "Shank2") %>% 
	filter(Accession == shank2) %>% 
	filter(Treatment != "QC")
fm1 <- lm(fx, data = df_shank2)
LT1 <- getContrast(fm1, "KO","WT")
res1 <- lmTestContrast(fm1, LT1) %>% 
	mutate(Contrast=gsub("Condition","",Contrast))
res1 %>% knitr::kable()

# fit a model: Shank3
df_shank3 <- tidy_prot %>% 
	filter(Genotype == "Shank3") %>%
	filter(Accession == shank3, Treatment != "QC")
fm2 <-  lm(fx, data = df_shank3)
LT2 <- getContrast(fm2, "KO","WT")
res2 <- lmTestContrast(fm2, LT2) %>% 
	mutate(Contrast=gsub("Condition","",Contrast))
res2 %>% knitr::kable()

# fit a model: Syngap1
df_syngap1 <- tidy_prot %>% 
	filter(Genotype == "Syngap1") %>% 
	filter(Accession == syngap1, Treatment != "QC")
fm3 <- lm(fx, df_syngap1)
LT3 <- getContrast(fm3, "HET","WT")
res3 <- lmTestContrast(fm3, LT3) %>% 
	mutate(Contrast=gsub("Condition","",Contrast))
res3 %>% knitr::kable()

# fit a model: Ube3a
df_ube3a <- tidy_prot %>% 
	filter(Genotype == "Ube3a") %>%
	filter(Accession == ube3a, Treatment != "QC")
fm4 <-  lm(fx, data = df_ube3a)
LT4 <- getContrast(fm4, "KO","WT")
res4 <- lmTestContrast(fm4, LT4) %>% 
	mutate(Contrast=gsub("Condition","",Contrast))
res4 %>% knitr::kable()


## ---- Assess differential abundance with lmTestContrast

## BASICALLY WE CONSIDER EACH EXPERIMENT SEPERATELY
## DO WE NEED IRS THEN?

# for all four tissue specific experiments
# NOTE: we specify a contrast for each call to testProtContrast!
res_list <- list()

# perform statistical analysis for Shank2
res_list[["Shank2"]] <- testProtContrast(tidy_prot %>% filter(Genotype == "Shank2"), fx, LT1)

# perform statistical analysis for Shank3
res_list[["Shank3"]] <- testProtContrast(tidy_prot %>% filter(Genotype == "Shank3"), fx, LT2)

# perform statistical analysis for Syngap1
res_list[["Syngap1"]] <- testProtContrast(tidy_prot %>% filter(Genotype == "Syngap1"), fx, LT3)

# perform statistical analysis for Ube3a
res_list[["Ube3a"]] <- testProtContrast(tidy_prot %>% filter(Genotype == "Ube3a"), fx, LT4)

# NOTE: ^Warnings about unable to moderate s2_prior, therefor s2_prior is set to 0

# collect results as single df, annotate with symbols and entrez ids
#results <- do.call(rbind, res_list) %>% annotProts(gene_map)
results <- bind_rows(res_list,.id="Genotype") %>% annotProts(gene_map)

# check nsig
# basically the same results with/without irs, which makes sense
results %>% 
	group_by(Genotype) %>% 
	summarize(nSig = sum(FDR<0.05), .groups="drop") %>% 
	knitr::kable()


## ---- save output data 
res_list <- results %>% group_by(Genotype) %>% group_split()
names(res_list) <- c("Shank2","Shank3","Syngap1","Ube3a")

myfile <- file.path(root,"tables", paste0("SP_",tissue,"-Results.xlsx"))
write_excel(res_list, myfile)
message("saved:", myfile)

myfile <- file.path(root,"data", paste0(tolower(tissue),"_results.rda"))
save(results, file=myfile, version=2)
message("saved: ", myfile)
