#' edgeR_GSE
#'
#' A Function for performing GO and KEGG GSE.
#' Requires an internet connection!
#' Can operate on qlf or ET objects.
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
#' edgeR_GSE(qlf, FDR = 0.05, filter = TRUE)
edgeR_GSE <- function(qlf, FDR = 0.05, filter = TRUE) {
  GSE_results <- list()
  # Map Uniprot IDs to Entrez IDs:
  Entrez <- mapIds(org.Mm.eg.db,
    keys = rownames(qlf),
    column = "ENTREZID",
    keytype = "UNIPROT",
    multiVals = "first"
  )
  Entrez[is.na(Entrez)] <- paste("not_mapped", c(1:sum(is.na(Entrez))), sep = "_")
  # Perform GO enrichment testing.
  GO <- goana(qlf, geneid = Entrez, species = "Mm", FDR = FDR)
  # Calculate adjusted p.values.
  GO$P.Adj.Up <- p.adjust(GO$P.Up, method = "hochberg")
  GO$P.Adj.Down <- p.adjust(GO$P.Down, method = "hochberg")
  # Perform KEGG enrichment testing.
  KEGG <- kegga(qlf, geneid = Entrez, species = "Mm", FDR = FDR)
  KEGG$P.Adj.Up <- p.adjust(KEGG$P.Up, method = "hochberg")
  KEGG$P.Adj.Down <- p.adjust(KEGG$P.Down, method = "hochberg")
  # Eliminate insignificant results if filter=TRUE.
  if (filter == TRUE) {
    GO <- GO[GO$P.Adj.Up < 0.05 | GO$P.Adj.Down < 0.05, ]
    KEGG <- KEGG[KEGG$P.Adj.Up < 0.05 | KEGG$P.Adj.Down < 0.05, ]
  }
  GSE_results$GO <- GO
  GSE_results$KEGG <- KEGG
  return(GSE_results)
}
