#' map_ids
#'
#' map gene identifiers given a dataframe with gene mapping data.
#'
#' @export

map_ids <- function(ids,gene_map,input_format,output_format) {
	ids <- as.character(ids)
	idx <- match(ids,gene_map[[input_format]])
	new_ids <- gene_map[[output_format]][idx]
	return(new_ids)
}
