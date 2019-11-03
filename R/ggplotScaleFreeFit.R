#' ggplotScaleFreeFit
#'
#' Function for plotting WGCNA powers.
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
#' ggplotScaleFreeFit(sft)
ggplotScaleFreeFit <- function(sft) {
  require(ggplot2)
  # Gather the data, calculate scale free fit.
  data <- sft$fitIndices
  data$fit <- -sign(data$slope) * data$SFT.R.sq
  # Generate Scale free topology plot.
  plot1 <- ggplot(data, aes(x = Power, y = fit)) +
    geom_text(aes(label = Power), color = "red") +
    ggtitle("Scale independence") +
    xlab(expression(Soft ~ Threshold ~ Power ~ (beta))) +
    ylab(expression(Scale ~ Free ~ Topology ~ (R^2))) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", size = 0.6) +
    # geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray", size = 0.6) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  # Generate mean connectivity plot.
  plot2 <- ggplot(data, aes(x = Power, y = mean.k.)) +
    geom_text(aes(label = Power), color = "red") +
    ggtitle("Mean Connectivity") +
    xlab(expression(Soft ~ Threshold ~ Power ~ (beta))) +
    ylab(expression(Mean ~ Connectivity ~ (k))) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  data_return <- list(plot1, plot2)
  names(data_return) <- c("ScaleFreeFit", "MeanConnectivity")
  return(data_return)
}
