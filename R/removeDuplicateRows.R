#' removeDuplicateRows
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
#' ## Function: removed duplicate rows.
#' # Removes rows that are duplicate. e.g. A|B and B|A
removeDuplicateRows <- function(dataframe, colA = 1, colB = 2) {
  df <- dataframe[, c(colA, colB)]
  df <- as.data.frame(apply(df, 2, function(x) as.character(x)))
  df2 <- df[!duplicated(data.frame(list(do.call(pmin, df), do.call(pmax, df)))), ]
  idx <- as.numeric(rownames(df2))
  result <- dataframe[idx, ]
  return(result)
}
