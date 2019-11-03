#' plotpdf
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
#' # ## Define function: plotpdf(figure_name)
#' # This is a simple function for plotting figure as PDF format using CairoPDF.
#' # Create plot object with pryr %<a-%.
plotpdf <- function(figure_name, plot_object) {
  CairoPDF(file = paste(savefigs_path, "/", figure_name, sep = ""), width = 12, height = 8, paper = "special")
  plot_object
  dev.off()
}
