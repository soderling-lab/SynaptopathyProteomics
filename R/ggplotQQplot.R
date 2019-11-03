#' ggplotQQplot
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
#' Define function: ggplotQQplot()
ggplotQQplot <- function(data_in, title) {
  df <- na.omit(melt(as.data.frame(data_in)))
  plot <- ggplot(df, aes(sample = value)) + stat_qq() + stat_qq_line() +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}

