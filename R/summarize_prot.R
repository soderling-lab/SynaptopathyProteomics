summarize_prot <- function(tp){
	# Summarize to protein level by summing peptides.

	# Imports.
	suppressPackageStartupMessages({
		library(data.table)
		library(dplyr)
	})

	# Remove NA, and ungroup.
	tp$Intensity[is.na(tp$Intensity)] <- 0
	tp <- ungroup(tp)
	# Sum to protein level.
	proteins <- tp %>% group_by(Genotype,Sample,Channel,
				    Treatment,Accession) %>%
		summarize(Peptides = length(Intensity),
			  Intensity = sum(Intensity,na.rm=TRUE))
	# Convert 0 back to NA.
	proteins$Intensity[proteins$Intensity==0] <- NA

	return(proteins)
}
