#' formatDEP
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
#' # ## formatDEP
#'
#' # Reformats the data for input into DEP package
formatDEP <- function(data_in) {
  #  Add ID and names columns
  data_in$ID <- data_in$Accession
  data_in$name <- data_in$Accession
  ## Create SE object
  TMT_columns <- grep("F", colnames(data_in)) # Intensity column numbers
  long.cols <- colnames(data_in)[TMT_columns]
  long.cols <- gsub("\\ ", "", long.cols)
  long.cols <- gsub("\\..", "", long.cols)
  vars <- colsplit(long.cols, ":", c("A", "B", "C"))
  vars <- vars[, 3]
  vars <- colsplit(vars, ",", c("channel", "sample", "condition", "proteomicsID", "genotype")) # Cortex format
  vars$sample <- NULL
  vars$replicate <- rep(c(1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4), 4)
  label <- paste(vars$genotype, vars$condition, ".", vars$replicate, sep = "")
  # Create SE object
  colnames(data_in)[TMT_columns] <- label
  data_se <- make_se_parse(data_in, TMT_columns)
  return(data_se)
}
