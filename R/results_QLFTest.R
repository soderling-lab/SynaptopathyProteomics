#' results_QLFTest
#'
#' A function to perform QLF test with GLM fit.
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
#' results_QLFTest(fit, comparison, alpha = 0.05)
results_QLFTest <- function(fit, comparison, alpha = 0.05) {
  # test comparisons with QLFtest, add FDR, extract results, and sort by pvalue.
  qlf <- glmQLFTest(fit, contrast = comparison)
  res <- qlf
  geno <- strsplit(strsplit(res$comparison, "\\.")[[1]][2], "\\ ")[[1]][1]
  res <- topTags(res, n = Inf, sort.by = "none")$table
  res <- res[order(res$PValue), ]
  # Determine number of significant results.
  summary_table <- summary(decideTests(qlf))
  # Convert logCPM to percent WT.
  res$logCPM <- round(100 * (2^res$logFC), 2)
  colnames(res)[2] <- "%WT"
  colnames(res)[3] <- "F Value"
  # Categorize candidates by FDR.
  res$candidate <- "no"
  res[which(res$FDR <= 0.10 & res$FDR > 0.05), dim(res)[2]] <- "low"
  res[which(res$FDR <= 0.05 & res$FDR > 0.01), dim(res)[2]] <- "med"
  res[which(res$FDR <= 0.01), dim(res)[2]] <- "high"
  res$candidate <- factor(res$candidate, levels = c("high", "med", "low", "no"))
  # Map Uniprot IDs to Gene names
  Uniprot_IDs <- sapply(strsplit(rownames(res), "\\|"), "[", 2)
  if (any(is.na(Uniprot_IDs))) {
    Uniprot_IDs <- rownames(res)
  }
  symbol <- mapIds(org.Mm.eg.db, keys = Uniprot_IDs, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
  res <- add_column(res, symbol, .before = 1)
  colnames(res)[1] <- "Symbol"

  # Extract comparison name from comparison (contrast) matrix
  name <- gsub(" ", "", strsplit(colnames(comparison), "\\.|-")[[1]][2])
  # Subset significant results.
  res_sub <- subset(res, res[, 5] < alpha)
  colnames(res)[c(2:6)] <- paste(geno, colnames(res)[c(2:6)])
  results_list <- list(summary_table, res, qlf, res_sub)
  names(results_list) <- c("summary_table", "QLFTest_results", "QLF", "Sig_results")
  return(results_list)
}
