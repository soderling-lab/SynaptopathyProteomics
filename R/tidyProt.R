tidyProt <- function(raw_data,intensity.cols,species=NULL,
		     samples=NULL,summary=FALSE){
	# Function to tidy proteomics data.
	# Imports.
	suppressPackageStartupMessages({
		library(data.table)
		library(TBmiscr)
		library(tibble)
		library(dplyr)
	})
	# Collect relevant columns. 
	#colNames <- colnames(PDdata)
	#colIDs <- c(grep(accessionCol,colNames),
	#	    grep(sequenceCol,colNames),
	#	    grep(modCol,colNames),
	#	    grep(dataCol,colNames))
	# Melt and fix column names.
	#dt <- PDdata %>% dplyr::select(colIDs) %>%  as.data.table()
	dt <- as.data.table(raw_data) %>% 
		melt(id.vars = intensity.cols,
		     variable.name="Sample",
		     value.name="Intensity",
		     variable.factor=FALSE) # Don't coerce to factor.
	# Remove proteins that do not coorespond to species of interest.
	idx <- grepl(paste0("OS=",species),dt$Description)
	if (!all(idx)) {
		n_out <- length(unique(dt$Accession[!idx]))
		msg <- paste(n_out,"proteins are not from", lquote(species),
			     "and will be removed.")
		warning(msg,call.=FALSE)
		dt <- dt %>% filter(grepl(paste0("OS=",species),Description))
	}
	# Insure zeros are NA.
	dt$Intensity[dt$Intensity==0] <- NA
	# Return tidy df.
	return(dt)
}
