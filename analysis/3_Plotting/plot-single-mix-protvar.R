#!/usr/bin/env Rscript

# title: SynaptopathyProteomics
# author: twab
# description: examine variance attributable to major experimental covariates

## ---- prepare the env

#library(SynaptopathyProteomics)
root <- "~/projects/SynaptopathyProteomics"
devtools::load_all(root)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(doParallel)
})


## ---- inputs

data(cortex) # tidy_prot

prot_col <- "Accession" # protein identifier column

# loop through all proteins, fit the model with all
# experimental covariates modeled as random effects
fx <- log2(Intensity) ~ (1|Condition) + (1|Batch)

# FIXME: see what the results look like for a single mixture, try dropping Genotype

# * Genotype = which genetic background
# * Batch = which Synaptosome purification batch
# * Condition = Shank2.WT and Syngap1.HET for example

# it don't matter how you you model things, residuals are very high for some
# reason

tidy_prot  <- tidy_prot %>% 
	# drop QC
	filter(!grepl("QC",Condition)) %>% 
	# input munge
	mutate(Condition=interaction(Genotype,Condition)) %>%
	# subset
	filter(Genotype == "Ube3a")



## ---- loop to do work

# register parallel backend
doParallel::registerDoParallel(parallel::detectCores() - 1)

# all prots
proteins <- unique(tidy_prot[[prot_col]])

# NOTE: this can take a couple minutes
pve_list <- foreach(prot = proteins) %dopar% {
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
  # FIXME: this line depends upon prot_col!
  fm <- lme4::lmer(fx, tidy_prot %>% subset(Accession == prot), 
		   control = lmer_control)
  vp <- getVariance(fm)
  pve <- vp/sum(vp)
  # FIXME: this line depends upon prot_col!
  pve_df <- as.data.table(t(pve)) %>% mutate(Accession = prot) %>% select(-Fixed)
  return(pve_df)
} #EOL


# collect results
prot_pve <- do.call(rbind, pve_list) %>%
	reshape2::melt(id.var=prot_col,
		       variable.name = "Parameter", 
		       value.name = "Variance")

# examine summary
prot_pve %>% group_by(Parameter) %>%
	summarize(Min=min(Variance), 
		  Median=median(Variance), 
		  Max=max(Variance), .groups="drop") %>%
	knitr::kable()


## ---- generate a plot

df <- prot_pve

# set the order
sort_by <- "Batch" # usually sort by max
xpos <- df %>% filter(Parameter == sort_by) %>%
  # FIXME: this line depends upon prot_col!
	arrange(desc(Variance)) %>% select(Accession) %>% unlist(use.names=FALSE)
  # FIXME: this line depends upon prot_col!
df$xpos <- match(df$Accession,xpos)

# generate the plot
plot <- ggplot(df)
plot <- plot + aes(x=xpos)
plot <- plot + aes(y=Variance)
plot <- plot + aes(group=Parameter)
plot <- plot + aes(stat=Variance)
plot <- plot + aes(fill=Parameter)
plot <- plot + geom_area()
plot <- plot + xlab("Protein")
plot <- plot + ylab("Percentage of Variance")
plot <- plot + theme(axis.line.x=element_line())
plot <- plot + theme(axis.line.y=element_line())
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(axis.text.x = element_text(color="black",size=11))
plot <- plot + theme(axis.text.x = element_text(angle=0,hjust=1,family="Arial"))
plot <- plot + theme(axis.text.y = element_text(color="black",size=11))
plot <- plot + theme(axis.text.y = element_text(angle=0,hjust=1,family="Arial"))
plot <- plot + theme(panel.border = element_rect(colour = "black", fill=NA))
plot <- plot + scale_x_continuous(expand=c(0,0))
plot <- plot + scale_y_continuous(expand=c(0,0))


ggsave("ube3a-batch.pdf",plot)


## --- save the plot

stop()

ggtheme()

setFont("Arial", font_path=file.path(root,"fonts"))

# NOTE: warnings probs bc plot is really huge and gets saved w/o labels
myfile <- file.path(root,"figs","Variance","protein-varpart.pdf")
ggsave(myfile, plot)
