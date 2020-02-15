enrichR <- function(gene_list,database){
	# Simple wrapper around enrichR.
	# Supress output from cat within enrichr with capture.output().
	capture.output({
		results <- lapply(gene_list, function(x) {
				  enrichr(x,db) 
		 })
	})
	return(results)
}
