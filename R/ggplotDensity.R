#' ggplotDensity
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
#' ## Define function: ggplotDensity(data_in,title)
ggplotDensity <- function(data_in, colID, title) {
  # The function ggplotDensity plots a density plot using ggplot.
  # The input data should be a data matrix. This data is reshaped with reshape2::melt
  # The title provided as input is used as the ggplot title.
  dm <- df2dm_TMT(data_in, colID)
  data_temp <- melt(log2(dm))
  colnames(data_temp) <- c("Accession", "Run", "Intensity")
  data_temp$Run <- as.factor(data_temp$Run)
  data_temp <- na.omit(data_temp)
  plot <- ggplot(data_temp, aes(x = Intensity, color = Run)) + geom_density() +
    ggtitle(title) +
    xlab("Log2 Intensity") +
    ylab("Density") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}
