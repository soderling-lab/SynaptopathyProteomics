#' ggplotMeanSdPlot
#'
#' plot mean SD plot with ggplot
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
#' ggplotMeanSdPlot(data_in, colID, title, log)
ggplotMeanSdPlot <- function(data_in, colID, title, log) {
  data_in <- data_in[, grepl(colID, colnames(data_in))]
  if (log == TRUE) {
    data_in <- as.matrix(log2(data_in))
  } else {
    data_in <- as.matrix(data_in)
  }
  ff <- tempfile()
  png(filename = ff)
  plot <- vsn::meanSdPlot(data_in)
  plot <- plot$gg
  dev.off()
  unlink(ff)
  plot <- plot + ggtitle(title) + theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold")
  )
  return(plot)
}
