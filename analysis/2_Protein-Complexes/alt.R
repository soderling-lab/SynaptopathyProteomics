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
root <- "~/projects/SynaptopathyProteomics"
devtools::load_all(root)


## ---- functions



testComplexes <- function(corum_prots, tidy_prot, LT) {
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore",calc.derivs=FALSE)
  res_list <- list()
  for (path in names(corum_prots)) {
    prots <- corum_prots[[path]]
    #fx <- log2(Intensity) ~ 0 + Condition + (1|Accession)
    fx <- log2(Intensity) ~ Condition + (1|Accession)
    fm <- lmerTest::lmer(fx, tidy_prot %>% subset(Accession %in% prots), control = lmer_control)
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

# predicted RAVE complex
data(list=c("rogdi", "wdr7", "dmxl2"))
rave <- c(dmxl2, rogdi, wdr7)

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


## data munge
idx = grepl("QC",tidy_prot$Condition)
tidy_prot$Condition[idx] <- "SPQC"
tidy_prot$Condition <- factor(tidy_prot$Condition, levels=unique(tidy_prot$Condition))
stopifnot(levels(tidy_prot$Condition)[1]=="SPQC")

# fit a model and get contrasts
fm <- lmerTest::lmer(fx, tidy_prot %>% filter(Accession %in% rave))
LT1 <- getContrast(fm, "Shank2.KO", "Shank2.WT")
LT2 <- getContrast(fm, "Shank3.KO", "Shank3.WT")
LT3 <- getContrast(fm, "Syngap1.HET", "Syngap1.WT")
LT4 <- getContrast(fm, "Ube3a.WT", "Ube3a.KO")


## ---- do the analysis

message("Analyzing ", length(paths), " protein complexes.")

results <- list()
results[["Shank2"]] <- testComplexes(filt_prots,tidy_prot, LT1)
results[["Shank3"]] <- testComplexes(filt_prots,tidy_prot,LT2)
results[["Syngap1"]] <- testComplexes(filt_prots,tidy_prot,LT3)
results[["Ube3a"]] <- testComplexes(filt_prots,tidy_prot,LT3)

# combine
results_df <- results %>% bind_rows(.id="Genotype")

#results[["Shank2"]] %>% filter(FDR<0.05) # 1 glurdelta (6prots)
#results[["Shank3"]] %>% filter(FDR<0.05) # none for cortex
#results[["Syngap1"]] %>% filter(FDR<0.05) # none for cortex
#results[["Ube3a"]] %>% filter(FDR<0.05) # none for cortex


## ---- save the results

myfile <- file.path(root,"tables",paste0("SP_",tissue,"-complex-results.xlsx"))
write_excel(results,myfile)
message("saved: ", myfile)

results <- results_df
myfile <- file.path(root,"data",paste0(tolower(tissue),"_complex_results.rda"))
save(results,file=myfile, version=2)
message("saved: ", myfile)
