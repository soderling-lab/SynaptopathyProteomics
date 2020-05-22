#!/usr/bin/env Rscript

# Add header to functions if its missing.

# Directories.
here <- getwd()
root <- dirname(here)
funcdir <- file.path(root,"R")

# Get all functions.
myfiles <- list.files(funcdir,pattern="*.R",full.names=TRUE)
myfun <- lapply(myfiles,readLines)
names(myfun) <- myfiles

# Function to check if a function has appropriate header.
hasHeader <- function(myfun,func_name) {
	func_header <- unlist(strsplit(myfun[1],"\\ "))
	checks <- func_header[1] == "#'" & func_header[2] == func_name
	return(checks)
}

# Fix functions without header.
func_names <- gsub(".R","",basename(names(myfun)))
has_header <- mapply(hasHeader,myfun,func_names)
missing_header <- names(has_header[!has_header])

# Load blank header.
blank_header <- readLines(file.path(here,"function_header.R"))

# Functions missing headers.
path2fun <- names(myfun[missing_header])

sapply(path2fun,function(x) { fixHeader(x,blank_header) })

fixHeader <- function(path2fun,blank_header){
	func_name <- gsub(".R","",basename(path2fun))
	func_header <- blank_header
	func_header[1] <- gsub("myfunction",func_name,func_header[1])
	func_header[18] <- gsub("myfunction", paste0(func_name,"()"),func_header[18])
	func_header[19] <- gsub("myfunction", func_name,func_header[19])
	func_header <- func_header[-20]
	func_body <- readLines(path2fun)
	func_body <- func_body[-grep("\\<-*function\\(",func_body)]
	output <- c(func_header,func_body)
	writeLines(output,path2fun)
}

