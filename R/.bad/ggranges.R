#' ggranges
#'
#' get x and y limits of plot and top right and top left corner positions
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
#' ggranges(plot)
ggranges <- function(plot) {
  # Calculate ranges.
  df <- data.frame(
    ggplot_build(plot)$layout$panel_params[[1]][1],
    ggplot_build(plot)$layout$panel_params[[1]][8]
  )
  df <- rbind(df, df[2, ] - df[1, ])
  rownames(df) <- c("min", "max", "delta")
  # TopRight
  TopRight <- data.frame(
    xmin = df$x.range[1] + 0.75 * df$x.range[3],
    xmax = df$x.range[2],
    ymin = df$y.range[1] + 0.55 * df$y.range[3],
    ymax = df$y.range[2]
  )
  # BottomRight
  BottomRight <- data.frame(
    xmin = df$x.range[1] + 0.75 * df$x.range[3],
    xmax = df$x.range[2],
    ymin = df$y.range[1] - 0.55 * df$y.range[3],
    ymax = df$y.range[2]
  )
  # Combine in a list.
  out <- list(df, TopRight, BottomRight)
  names(out) <- c("limits", "TopRight", "BottomRight")
  return(out)
}
