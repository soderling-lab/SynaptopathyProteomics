#' ggplotBoxPlotv2
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
#' ## boxplot function.
ggplotBoxPlotv2 <- function(data_in, colID, traits, colors, title) {
  sampleIndex <- as.data.frame(do.call(rbind, strsplit(colnames(data_in), "\\.")))
  colnames(sampleIndex) <- c("batch", "channel")
  batchIndex <- unique(sampleIndex$batch)
  model_colors <- as.data.frame(cbind(batchIndex,
    Model = unique(traits$Model),
    Colors = colors
  ))
  colors <- rep(model_colors$Colors, as.vector(rowSums(table(sampleIndex))))

  dm <- df2dm_TMT(data_in, colID)
  data_temp <- melt(dm)
  colnames(data_temp) <- c("Accession", "Run", "Intensity")
  data_temp$Run <- as.factor(data_temp$Run)
  data_temp <- na.omit(data_temp)

  plot <- ggplot(data_temp, aes(x = Run, y = Intensity, fill = Run)) +
    geom_boxplot(outlier.colour = "black", outlier.shape = 20, outlier.size = 1) +
    scale_fill_manual(
      values = colors,
      name = "Genotype",
      breaks = c(cumsum(as.vector(rowSums(table(sampleIndex))))),
      labels = paste(model_colors$Model, " (", model_colors$batchIndex, ")", sep = "")
    ) +
    ggtitle(title) +
    xlab("TMT Run") +
    ylab("Log2 Intensity") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}
