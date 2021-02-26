#!/usr/bin/env Rscript

# title: SynaptopathyProteomics
# author: twab 
# description: analyzing corum protein complexes with linear mixed model


## ---- imports

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(geneLists)
  library(tidyProt)
})

#library(SynaptopathyProteomics)
root <- "~/projects/soderling-lab/SynaptopathyProteomics"
devtools::load_all(root)


## ---- functions

testComplexes <- function(corum_prots, tidy_prot, genotype) {
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore",calc.derivs=FALSE)
  res_list <- list()
  for (path in names(corum_prots)) {
    prots <- corum_prots[[path]]
    fx <- log2(Intensity) ~ 0 + Condition + (1|Accession)
    df <- tidy_prot %>% filter(Genotype == genotype, Accession %in% prots)
    fm <- lmerTest::lmer(fx, data = df, control = lmer_control)
    LT <- getContrast(fm, "KO|HET", "WT")
    res <- lmerTestContrast(fm, LT) %>% 
	    mutate(nProts=length(prots)) %>%
	    mutate(Contrast = gsub("Condition","",Contrast))
    res_list[[path]] <- res
  } #EOL
  results <- res_list %>% 
	  bind_rows(.id="Pathway") %>% 
	  mutate(FDR = p.adjust(Pvalue, method="fdr")) %>%
	  mutate(Padjust = p.adjust(Pvalue, method="bonferroni")) %>%
	  mutate(candidate = Padjust < 0.05) %>%
	  arrange(Pvalue)
  return(results)
} # EOF


## ---- main

tissue = parseArgs()

# load the data
data(gene_map) # gene_map
data(list=tolower(tissue)) # e.g. data(cortex); tidy_prot
data(list=paste0(tolower(tissue),"_results")) # e.g. data(cortex_results); results

# NOTE: we must subset:
# tidy_prot %>% filter(Genotype == "Shank2")

# predicted RAVE complex
data(rogdi)
dmxl2 <- gene_map$uniprot[match("Dmxl2",gene_map$symbol)]
wdr7 <- gene_map$uniprot[match("Wdr7",gene_map$symbol)]
rave <- c(dmxl2, rogdi, wdr7)

# not efficient
#corum_prots <- lapply(corum,function(x) getIDs(x,from="entrez",to="uniprot",species="mouse"))
data(corum_prots) 

corum_prots[["Predicted Complex: RAVE"]] <- rave

## keep complexes with >50% coverage and at least 2 proteins
df <- data.table(pathway=names(corum_prots), 
		 complex_size = sapply(corum_prots,length), 
		 n_identified=sapply(corum_prots,function(x) sum(x %in% tidy_prot$Accession)))
df <- df %>% mutate(coverage = n_identified/complex_size) %>% 
	mutate(keep=coverage>0.5 & n_identified>2)

# we will analyze the following protein complexes (pathways):
paths <- df$pathway[df$keep]

filt_prots <- corum_prots[paths]

## ---- example fit Rave

#fx = log2(rel_Intensity) ~ (1|Condition) + (1|Accession)
#df = tidy_prot %>% filter(Treatment != "QC") %>%
#	filter(Accession %in% rave) %>%
#	group_by(Accession) %>%
#	mutate(rel_Intensity = Intensity/sum(Intensity))
#fm = lmerTest::lmer(fx, df)
#vp = getVariance(fm)
#pve = vp/sum(vp)
#pve


## ---- do the analysis

message("Analyzing ", length(paths), " protein complexes.")

results <- list()
results[["Shank2"]] <- testComplexes(filt_prots,tidy_prot,"Shank2")
results[["Shank3"]] <- testComplexes(filt_prots,tidy_prot,"Shank3")
results[["Syngap1"]] <- testComplexes(filt_prots,tidy_prot,"Syngap1")
results[["Ube3a"]] <- testComplexes(filt_prots,tidy_prot,"Ube3a")

results_df <- results %>% bind_rows(.id="Genotype")


## ---- save the results

myfile <- file.path(root,"tables",paste0("SP_",tissue,"-complex-results.xlsx"))
write_excel(results,myfile)
message("saved: ", myfile)

results <- results_df
myfile <- file.path(root,"data",paste0(tolower(tissue),"_complex_results.rda"))
save(results,file=myfile, version=2)
message("saved: ", myfile)
