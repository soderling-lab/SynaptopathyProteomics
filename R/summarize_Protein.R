#' summarize_Protein
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
#' ## summarize_Protein
#' # Summarize proteins by summing all peptides for a unique Accession identifier.
summarize_Protein <- function(peptide_data) {
  # Add column for peptides, summarize using dplyr::summarize_all(sum)
  Peptides <- rep(1, nrow(peptide_data))
  temp_data <- add_column(peptide_data, Peptides, .after = 5)
  tmt_cols <- getCols(temp_data, "Abundance")
  temp_data <- temp_data[, c(2, 6, tmt_cols)]
  prot_data <- temp_data %>%
    group_by(Accession) %>%
    summarise_all(funs(sum), na.rm = TRUE)
  # Replace 0 with NA
  prot_data[prot_data == 0] <- NA
  return(prot_data)
}
