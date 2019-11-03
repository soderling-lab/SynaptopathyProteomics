#' normalize_TMM
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
#' ## normalize_TMM(data_in,colID)
#' # This function performs tmm normalization using the edgeR package on the columns
#' # specified by tmt_cols.
normalize_TMM <- function(data_in, groups) {
  for (i in 1:length(groups)) {
    colID <- groups[i]
    # Get data
    tmt_cols <- grep(colID, colnames(data_in))
    dm <- as.matrix(data_in[, tmt_cols])
    # replace NA with 0 (TMM cannot have missing values (NA))
    logic <- is.na(dm)
    if (length(logic[logic == TRUE]) > 0) {
      print("Warning: missing values (NA) are not tolerated. These will be replaced with 0")
    }
    dm[logic] <- 0
    data_work <- dm
    # TMM Normalization.
    factors_tmm <- calcNormFactors(data_work)
    dm_tmm <- sweep(data_work, 2, factors_tmm, FUN = "/")
    dm_tmm[logic] <- NA # Insure 0 values are now NA
    # Write to data frame
    data_in[, tmt_cols] <- dm_tmm
  }
  return(data_in)
}
