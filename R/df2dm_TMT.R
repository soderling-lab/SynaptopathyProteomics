#' df2dm_TMT
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
#' # ## df2dm_TMT
#'
#' # This is a simple function to convert data frame to dm for plotting purposes.
df2dm_TMT <- function(df, colID) {
  tmt_cols <- grep(colID, colnames(df))
  dm <- as.matrix(df[, tmt_cols])
  rownames(dm) <- if ("Accession" %in% df) {
    df$Accession
  } else {
    rownames(df)
  } #*** EBD edited
  colnames(dm) <- c(1:ncol(dm))
  return(dm)
}
