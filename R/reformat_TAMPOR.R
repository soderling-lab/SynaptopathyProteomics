reformat_TAMPOR <- function(tp,samples){
	# Coerce tidy data to matrix.
	tp <- filt_protein
	tp_in <- tp <- as.data.table(tp)
	dt <- tp %>% dcast(Accession ~ Sample,value.var="Intensity")
	dm <- as.matrix(dt,rownames="Accession")
	# Row names are Symbol|Accession
	idx <- match(rownames(dm),gene_map$uniprot)
	rownames(dm) <- paste(gene_map$symbol[idx],rownames(dm),sep="|")
	# Column names are batch.channel
	idx <- match(colnames(dm), samples$Sample)
	batch <- paste0("b",as.numeric(as.factor(samples$Genotype)))
	channel <- samples$Channel
	colnames(dm) <- paste(batch,channel,sep=".")[idx]
	# Reorder based on batch.channel.
	dm <- dm[,order(colnames(dm))]
	return(dm)
}
