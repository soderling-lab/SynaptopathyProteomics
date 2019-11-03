#' annotate_Entrez
#'
#' function_description
#'
#' @param
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @import
#'
#' @export
#'
#' @examples
#' function_name(param1, ...)
#' #-------------------------------------------------------------------------------
#' # Annotate topTags table with entrez and gene symbols.
#' # based on rownames being Uniprot IDS.
#' # Function:
annotate_Entrez <- function(y_TT) {
  Uniprot <- rownames(y_TT)
  # Map Uniprot IDs to Entrez IDs:
  Entrez <- mapIds(org.Mm.eg.db, keys = Uniprot, column = "ENTREZID", keytype = "UNIPROT", multiVals = "first")
  Gene <- mapIds(org.Mm.eg.db, keys = Uniprot, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
  # Add columns for Entrez ID and gene symbols.
  y_TT <- add_column(y_TT, Entrez, .before = 1)
  y_TT <- add_column(y_TT, Uniprot, .before = 1)
  y_TT <- add_column(y_TT, Gene, .after = 2)
  return(y_TT)
}
