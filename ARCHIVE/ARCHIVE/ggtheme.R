#' ggtheme
#'
#' Set a theme for ggplots.
#'
#' @param none
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' ggtheme()
ggtheme <- function() {
  suppressPackageStartupMessages({
    require(ggplot2)
  })
  ggtheme <- theme_gray() +
    theme(
      plot.title = element_text(color = "black", size = 11, face = "bold", hjust = 0.5),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold"),
      axis.text.x = element_text(color = "black", size = 11, angle = 0, hjust = 1.0)
    )
  theme_set(ggtheme)
}
