#' ggplotMDSv3
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
#' ## Define function: ggplotMDS(data_in, colID, traits, title)
ggplotMDSv3 <- function(data_in, colID, traits, title) {
  idx <- match(colnames(data_in), rownames(traits))
  # colnames(data_in) <- traits$SampleType[idx]
  colnames(data_in) <- traits$Sample.Model[idx]
  ff <- tempfile()
  png(filename = ff)
  data_MDS <- plotMDS(data_in)
  x <- data_MDS$x
  y <- data_MDS$y
  dev.off()
  unlink(ff)
  dm_MDS <- cbind(x, y)
  Condition <- rownames(dm_MDS)
  df_MDS <- as.data.frame(cbind(x, y))
  plot <- ggplot(df_MDS, aes(x, y, color = Condition)) + geom_text(aes(label = Condition)) +
    ggtitle(title) + xlab("Leading LogFC dim 1") + ylab("Leading LogFC dim 2") +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}
