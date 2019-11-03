#' clean_DEP
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
#' # ## Define function: clean_DEP
#' # This function just compresses the work for cleaning up the DEP limma results.
clean_DEP <- function(data_in, alpha, lfc) {
  # Denote significant proteins: no logFC cut-off, alpha = 0.1
  dep <- add_rejections(data_diff, alpha = alpha, lfc = lfc)

  # Generate a results table
  results_limma <- get_results(dep)

  # Drop un-wanted columns
  keep <- c(1, grep("p.val", colnames(results_limma)), grep("p.adj", colnames(results_limma)))
  results_limma <- results_limma[, keep]

  # Calculate Log2FC from data_in
  data_temp <- assay(data_in)
  log2fc <- matrix(NA, nrow(data_temp), 4)
  group <- c("Shank2", "Syngap1", "Ube3a", "Shank3")

  # Loop to calc log2fc (Avg.KO-Avg.WT)
  for (i in 1:4) {
    sub <- grep(group[i], colnames(data_temp))
    data_sub <- data_temp[, sub]
    log2fc[, i] <- rowMeans(data_sub[, 5:8]) - rowMeans(data_sub[, 1:4])
  }

  # Clean up result
  log2fc <- as.data.frame(log2fc)
  colnames(log2fc) <- c("Shank2_Log2FC", "Syngap1_Log2FC", "Ube3a_Log2FC", "Shank3_Log2FC")
  rownames(log2fc) <- rownames(data_temp)
  order <- match(results_limma$name, rownames(log2fc))
  log2fc <- log2fc[order, ]

  # Combine with limma results
  results_limma <- cbind(results_limma, log2fc)

  # Get Uniprot_IDs
  Uniprot_IDs <- rownames(results_limma)
  head(Uniprot_IDs)

  # Map Uniprot IDs to Gene names
  symbol <- mapIds(org.Mm.eg.db, keys = Uniprot_IDs, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
  results_limma <- add_column(results_limma, symbol, .after = 1)

  # Need to write to results.
  return(results_limma)
}
