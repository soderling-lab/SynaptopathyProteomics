#' ggplotBCV
#'
#' plot biological coefficients of variation
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
#' ggplotBCV(y_DGE)
ggplotBCV <- function(y_DGE) {
  y <- y_DGE
  tag_x <- y$AveLogCPM
  tag_y <- sqrt(y$tagwise.dispersion)
  common <- sqrt(y$common.dispersion)
  trend_y <- sqrt(y$trended.dispersion)
  trend_x <- y$AveLogCPM
  o <- order(trend_x)
  fill <- rep(1, nrow = length(trend_x), ncol = 1)
  df <- as.data.frame(cbind(tag_x, tag_y, trend_x, trend_y), fill)

  plot <- ggplot(df, aes(tag_x, tag_y, colour = "Tagwise")) + geom_point(size = 1) +
    geom_line(mapping = aes(trend_x[o], trend_y[o], colour = "Trend"), size = 0.75) +
    geom_hline(aes(yintercept = common, colour = "Common"), linetype = "solid", size = 0.75) +
    scale_colour_manual(name = "Dispersion", values = c("blue", "black", "red")) +
    xlab("Average log CPM") + ylab("Biological coefficient of variation") +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}
