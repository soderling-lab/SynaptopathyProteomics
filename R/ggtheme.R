#' ggtheme
#'
#' Default theme for gg plots.
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
#' ggtheme()
ggtheme <- function() {
  require(ggplot2)
  ggtheme <- theme(
    plot.title = element_text(color = "black", size = 11, face = "bold", hjust = 0.5),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold")
  )
}
