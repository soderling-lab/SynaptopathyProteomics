#' ggplotDensityv2
#'
#' plot sample density plots
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
#' ggplotDensity(data_in, title)
ggplotDensityv2 <- function(data_in, colID, colors, traits, title) {
  sampleIndex <- as.data.frame(do.call(rbind, strsplit(colnames(cleanDat), "\\.")))
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
  plot <- ggplot(data_temp, aes(x = Intensity, color = Run)) + geom_density() +
    ggtitle(title) +
    xlab("Log2 Intensity") +
    ylab("Density") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    ) +
    scale_color_manual(values = colors)
  return(plot)
}
