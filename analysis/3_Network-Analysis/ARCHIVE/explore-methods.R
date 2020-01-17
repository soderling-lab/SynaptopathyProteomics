#!/usr/bin/env Rscript

## Which method to use for GO and PPI graph clustering?

library(data.table)

# Directories.
here <- getwd()
root <- dirname(dirname(here))
rdatdir <- file.path(root,"rdata")

# Load partitions files.
file_prefix <- c("3_PPI_","3_GO_")
file_suffix <- "VertexPartition_partitions.csv"
partition_methods <- c("CPM","Modularity","RBConfiguration","RBER","Significance")
parts_files <- apply(expand.grid(file_prefix,partition_methods,file_suffix), 1, paste, collapse="")
data <- lapply(parts_files, function(myfile) {
		       # Check if file exists.
		       if (file.exists(file.path(rdatdir,myfile))) { 
			       data <- fread(file.path(rdatdir,myfile),header=TRUE,drop=1) 
			       return(data)
		       } else {
			       message(paste("Warning", myfile,"does not exist!"))
			       return(NA)
		       }
		 })
names(data) <- tools::file_path_sans_ext(parts_files) # Name.
data <- data[order(names(data))] # Sort.

# Loop through data.
# Best partition maximizes the number of clustered proteins.
results <- list()
for (partition in names(data)) {
	parts <- data[[partition]]
	method <- paste(sapply(strsplit(partition,"_"),"[",c(2,3)),collapse="_")
	not_clustered <- apply(parts,1,function(x) sum(x==0)/length(x))
	rbest <- seq(not_clustered)[not_clustered==min(not_clustered)]
	nMods <- sum(names(table(as.numeric(parts[rbest,])))!="0")
	results[[method]] <- c("method" = method,"rbest"=rbest,
	  "nMods"=nMods,"not_clustered"=not_clustered[rbest])
}
	
# Save data for examination.
df <- as.data.frame(do.call(rbind,results))
fwrite(df,"temp.csv")
