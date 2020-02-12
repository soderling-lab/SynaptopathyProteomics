# Is co-expression correlated with GO sem sim score?

# Imports.
library(getPPIs)
library(dplyr)

# Load Frela data.
here = getwd()
root = dirname(dirname(here))
downdir = file.path(root,"downloads")
datadir = file.path(downdir,"frela")
rdatdir = file.path(root,"rdata")

# Load Co-expression matrix.
adjm <- as.matrix(readRDS(file.path(rdatdir,"3_Cortex_Adjm.RData")))

# Load protein id map.
prot_map <- readRDS(file.path(rdatdir,"2_Protein_ID_Map.RData"))

# Remap adjm column names as Entrez Ids.
colNames <- prot_map$entrez[match(colnames(adjm),prot_map$ids)]
colnames(adjm) <- rownames(adjm) <- colNames

# Get Frela input.
#input_files <- list.files(file.path(datadir,"frela_input"),full.names=TRUE)
#input <- data.table::fread(input_files[1],header=FALSE)
#colnames(input) <- c("ProtA","ProtB")

# Get Frela result.
result_files <- list.files(file.path(datadir,"frela_results"),full.names=TRUE)
result <- data.table::fread(result_files[1])

# Map Uniprot IDs to Entrez.
uniprot <- unique(c(result$UniProt1,result$UniProt2))
entrez <- mapIDs(uniprot,from="uniprot",to="entrez",species="mouse")
names(entrez) <- uniprot

# Add Entrez columns to result.
result = tibble::add_column(result,
			    Entrez1 = entrez[result$UniProt1],
			    .after="UniProt1")
result = tibble::add_column(result,
			    Entrez2 = entrez[result$UniProt2],
			    .after="UniProt2")

# Remove unmapped rows.
out <- is.na(result$Entrez1) | is.na(result$Entrez2)
nout <- sum(out)
message(paste("rows removed:",nout))
result_filt = result[!out,]

# Add bicor stats.
bicor <- vector(mode="numeric",nrow(result))
pbar <- txtProgressBar(min=1,max=nrow(result),style=3)
for (i in seq(nrow(result))) {
	     setTxtProgressBar(pbar,i)
	     idx = match(result$Entrez1[i],colnames(adjm))
	     idy = match(result$Entrez2[i],colnames(adjm))
	     bicor[i] = adjm[idx,idy]
	     if (i == nrow(result)) { close(pbar); message("/n") }
}

# Add to results.
result$bicor <- bicor

# Coerce to numeric df.
stats = colnames(result)[c(7:14)]
dm = as.matrix(result %>% select(stats))
df = as.data.frame(apply(dm,2,function(x) as.numeric(x)))

# Calc average z.score.
cols <- grep("z\\.",colnames(df))
df$avg.z = apply(df[,cols],1,function(x) mean(x,na.rm=TRUE))


# Assess correlation...
# Collect data and coerce to numeric.
stats = colnames(df)
for (stat in stats){
	x = bicor
	y=as.numeric(df[[stat]])
	mydat = na.omit(as.matrix(cbind(x,y)))
	rho = cor(mydat[,"x"],mydat[,"y"],method="spearman")
	print(rho)
}


x = df$avg.z
y = df$bicor




