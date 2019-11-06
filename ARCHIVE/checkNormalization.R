#' checkNormalization
#'
#' This function calls several of the custom ggplot functions for checking
#' the data:
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
#' checkNormalization(data_in, traits, colors, title)
checkNormalization <- function(data_in, traits, colors, title) {
  # Generate plots.
  p1 <- ggplotBoxPlotv2(log2(data_in),
    colID = "b", traits = traits,
    colors = colors,
    title = title
  )
  p2 <- ggplotDensityv2(log2(data_in),
    colID = "b", traits = traits,
    colors = colors,
    title = title
  )
  p3 <- ggplotMDS(log2(data_in),
    colID = "b", traits,
    title = title
  ) + theme(legend.position = "none")
  p4 <- ggplotMeanSdPlot(log2(data_in), title = title)
  # Store in list.
  plot_list <- list(p1, p2, p3, p4)
  # Return plots.
  return(plot_list)
}
