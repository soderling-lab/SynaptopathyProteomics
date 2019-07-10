#!/usr/bin/env Rscript

library(styler)

# Get the users project and code directory.
dir <- unlist(strsplit(getwd(),"/"))
project_dir <- paste(dir[c(1:length(dir)-1)], collapse = "/")
code_dir <- paste(project_dir, "Code", sep = "/")

# Style the code directory. 
style_dir(code_dir) #styles all .R and/or .Rmd files in a directory.
