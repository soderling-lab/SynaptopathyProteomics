#' col2hex
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
#' function_name(param1, ... )
#-------------------------------------------------------------------------------
#' ## A Function for mixing colors.
#' Returns hex code for new color and prints plot showing the color in the console.

col2hex <- function(color, maxValue = 255) {
  z <- col2rgb(color)
  hex <- rgb(z[1], z[2], z[3], maxColorValue = maxValue)
  return(hex)
}

mixcolors <- function(color1, color2, ratio1 = 1, ratio2 = 1, plot = FALSE) {
  # Convert colors to RGB and mix (average).
  x <- col2rgb(color1)
  y <- col2rgb(color2)
  z <- (ratio1 * x + ratio2 * y) / 2
  # Convert to hex format.
  hex <- rgb(z[1], z[2], z[3], maxColorValue = 255)
  # Plot for visualizing new color.
  df <- as.data.frame(cbind(x = 1, y = 1, hex))
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_text(label = hex, colour = hex, size = 25) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  if (plot == TRUE) print(p)
  return(hex)
}

