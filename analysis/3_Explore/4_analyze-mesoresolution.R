#!/usr/bin/env Rscript

source('~/projects/SynaptopathyProteomics/renv.R')

rdatdir <- file.path(root,"rdata")

myfile <- file.path(rdatdir,"TMT_Cortex_CPM_partition.csv")
partitions <- data.table::fread(myfile,drop=1)

modules <- apply(partitions,1,function(x) split(colnames(partitions),x))
nmodules <- sapply(modules,length)

myfile <- file.path(rdatdir,"TMT_Cortex_CPM_quality.csv")
quality_dt <- data.table::fread(myfile,drop=1)
quality_dt$k <- nmodules

library(ggplot2)
ggplot(quality_dt,aes(x=Resolution,y=Modularity)) + geom_point()

ggplot(quality_dt,aes(x=Resolution,y=k)) + geom_point()
