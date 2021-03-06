#!/usr/bin/env Rscript

# title: SynaptopathyProteomics
# author: twab
# description: examine variance attributable to major experimental covariates

# generates plot for a given geno and tissue
geno <- "Shank2"
tissue <- "Cortex"

# define fx specifying experimental covariates
# cortex models:
fx <- log2(Intensity) ~ (1|Condition) + (1|Batch) + (1|Sex) + (1|Age)
#fx <- log2(Intensity) ~ (1|Condition) + (1|Batch) + (1|Sex) # Syngap1, Ube3a (no age)

# striatum models:
#fx <- log2(Intensity) ~ (1|Condition) + (1|Batch) + (1|Sex) + (1|Age) # Shank2
#fx <- log2(Intensity) ~ (1|Condition) + (1|Age) + (1|Sex) # Shank3 (no batch)
#fx <- log2(Intensity) ~ (1|Condition) + (1|Sex) # Syngap1 (no batch or age)
#fx <- log2(Intensity) ~ (1|Condition) + (1|Sex) + (1|Batch) # Ube3a (no age)

# protein identifier column
prot_col <- "Accession" 


## ---- prepare the env

#library(SynaptopathyProteomics)
root <- "~/projects/soderling-lab/SynaptopathyProteomics"
devtools::load_all(root)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyProt)
  library(data.table)
  library(doParallel)
})

# set plotting theme
ggtheme()
setFont("Arial", font_path=file.path(root,"fonts"))


## ---- inputs

data(list=tolower(tissue)) # tidy_prot
data(list=c("shank2","shank3","syngap1","ube3a"))

# loop through all proteins, fit the model with
# experimental covariates modeled as random effects

# * Genotype = which genetic background
# * Batch = which Synaptosome purification batch
# * Condition = Shank2.WT and Syngap1.HET, for example

# it don't matter how you you model things, residuals are very high for some
# reason... meaning there is some source of unexplained variance

# subset the data 
tidy_prot  <- tidy_prot %>% 
	# drop QC - don't include QC samples in modeling
	filter(!grepl("QC",Condition)) %>% 
	# input munge
	mutate(Condition=interaction(Genotype,Condition)) %>%
	# subset
	filter(Genotype == geno)

# examine a goi
fm <- lmerTest::lmer(fx, data=tidy_prot %>% subset(Accession==eval(parse(text=tolower(geno))))) 
vp <- getVariance(fm)
pve <- vp/sum(vp)
data.table(factor=names(pve),pve=pve) %>% knitr::kable()

# examine a random protein
prot = sample(unique(tidy_prot[[prot_col]]),1)
fm <- lmerTest::lmer(fx, data=tidy_prot %>% subset(Accession==prot)) 
vp <- getVariance(fm)
pve <- vp/sum(vp)
data.table(factor=names(pve),pve=pve) %>% knitr::kable()


## ---- loop to do work

# register parallel backend
doParallel::registerDoParallel(parallel::detectCores() - 1)

# all prots
proteins <- unique(tidy_prot[[prot_col]])

# NOTE: this can take a couple minutes
pve_list <- foreach(prot = proteins) %dopar% {
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
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
		  Mean = mean(Variance),
		  Median=median(Variance), 
		  SD=sd(Variance),
		  Max=max(Variance), .groups="drop") %>%
	knitr::kable()


## ---- generate a plot

# set the order
sort_by <- "Residual" # usually sort by max
xpos <- df %>% filter(Parameter == sort_by) %>%
	arrange(desc(Variance)) %>% select(Accession) %>% unlist(use.names=FALSE)
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
plot <- plot + ggtitle(paste(geno,tissue))

## --- save the plot

myfile <- file.path(root,"figs","variance",
		    paste0(tissue,"_",geno,"_varpart.pdf"))
ggsave(myfile, plot, height=5, width=5)
