#' plot2pdf
#'
#' This is a simple function for plotting figure as PDF format using CairoPDF.
#' You can link a standar plot object with pryr %<a-%.
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
#' plot2pdf(figure_name, plot_object)
plot2pdf <- function(figure_name, plot_object) {
  CairoPDF(
    file = file.path(savefigs_path, figure_name),
    width = 12, height = 8, paper = "special"
  )
  plot_object
  dev.off()
}
