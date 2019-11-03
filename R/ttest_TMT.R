#' ttest_TMT
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
#' # ## ttest_TMT
#' # This function just does all the work of ttest.
ttest_TMT <- function(data_in, groups) {
  data_temp <- as.data.frame(assay(data_in))
  data_out <- matrix(NA, nrow = nrow(data_temp), ncol = length(groups))
  data_out <- as.data.frame(data_out)
  rownames(data_out) <- rownames(data_temp)
  for (i in 1:length(groups)) {
    data_sub <- data_temp[, grep(groups[i], colnames(data_temp))]
    results_ttest <- apply(data_sub, 1, function(x) t.test(x[5:8], x[9:12]))
    ttest_pvalue <- unlist(lapply(results_ttest, function(x) x$p.value))
    ttest_fdr <- p.adjust(ttest_pvalue, method = "BH")
    data_out[, i] <- ttest_pvalue
    colnames(data_out)[i] <- paste(groups[i], "p.value", sep = "_")
    data_out <- add_column(data_out, ttest_fdr)
    colnames(data_out)[ncol(data_out)] <- paste(groups[i], "fdr", sep = "_")
  }
  return(data_out)
}
