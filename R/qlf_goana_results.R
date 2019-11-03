#' qlf_goana_results
#'
#' Get GO and KEGG results from qlf object.
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
#' qlf_goana_results(qlf_result)
qlf_goana_results <- function(qlf_result) {
  qlf <- qlf_result$QLF
  Uniprot_IDs <- sapply(strsplit(rownames(qlf), "\\|"), "[", 2)
  Entrez <- mapIds(org.Mm.eg.db, keys = Uniprot_IDs, column = "ENTREZID", keytype = "UNIPROT", multiVals = "first")
  Entrez[is.na(Entrez)] <- paste("not_mapped", c(1:sum(is.na(Entrez))), sep = "_")
  rownames(qlf) <- Entrez
  go <- goana(qlf, species = "Mm")
  go$FDR.Up <- p.adjust(go$P.Up, method = "BH")
  go$FDR.Down <- p.adjust(go$P.Down, method = "BH")
  keg <- kegga(qlf, species = "Mm")
  keg$FDR.Up <- p.adjust(keg$P.Up, method = "BH")
  keg$FDR.Down <- p.adjust(keg$P.Down, method = "BH")
  results <- list(go, keg)
  names(results) <- c("GO", "KEGG")
  return(results)
}
