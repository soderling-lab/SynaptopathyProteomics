#' ggplotPCAv2
#'
#' Generates a PCA plot.
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
#' ggplotPCA(data_in, traits)
ggplotPCAv2 <- function(data_in, traits, title = "2D PCA Plot") {
  # Perform PCA
  PC <- prcomp(t(na.omit(data_in)))$x[, 1:2]
  # Add annotations to PC data frame.
  PC <- as.data.frame(PC)
  PC$Label <- traits$Sample.Model[match(rownames(PC), rownames(traits))]
  colors <- traits$Color[match(rownames(PC), rownames(traits))]
  # Generate plot.
  plot <- ggplot(PC, aes(x = PC1, y = PC2)) +
    geom_text(aes(label = Label), color = colors) +
    ggtitle(title) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )

  # Generate dendrogram
  distMat <- dist(PC)
  sampleTree <- flashClust(t(distMat), method = "complete")
  dendro <- ggdendrogram(sampleTree, rotate = TRUE, theme_dendro = FALSE) + xlab("Sample") +
    ylab("Height") + ggtitle("Sample Clustering based on PCA") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )

  # Generate output
  data_out <- list(PC, plot, dendro)
  names(data_out) <- c("data_PCA", "plot", "dendro")

  return(data_out)
}
