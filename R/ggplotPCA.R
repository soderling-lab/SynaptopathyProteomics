#' ggplotPCA
#'
#' plot a PCA plot with ggplot
#'
#' @param data_in - the expression data
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @export
#'
#' @examples
#' ggplotPCA(data_in, traits, colors)
ggplotPCA <- function(data_in, traits, colors, 
		      colID = "Abundance", title = "2D PCA Plot") {

  # Get the numeric data.
  data_in <- as.data.frame(data_in)
  idy <- grepl(colID, colnames(data_in))
  data_in <- data_in[, idy]

  # Perform PCA. Remove any NA.
  pca <- prcomp(t(na.omit(data_in)))
  pca_summary <- as.data.frame(t(summary(pca)$importance))
  idx <- order(pca_summary[["Proportion of Variance"]],decreasing=TRUE)
  top2 <- head(rownames(pca_summary[idx,]),2)
  pve <- pca_summary[top2,"Proportion of Variance"]
  names(pve) <- top2
  PC <- pca$x[,top2]
  colnames(PC) <- c("x","y")

  # Add annotations to PC data frame.
  PC <- as.data.frame(PC)
  idx <- match(rownames(PC), traits$ColumnName)
  PC$Label <- paste(traits$Model, traits$SampleType, sep = "_")[idx]

  # X and y labs for plot:
  x_label <- paste0(names(pve)[1]," (",round(100*pve[1],2),"%)")
  y_label <- paste0(names(pve)[2]," (",round(100*pve[2],2),"%)")

  # Generate plot.
  plot <- ggplot(PC, aes(x,y)) +
    geom_text(aes(label = Label), color = colors) +
    ggtitle(title) + xlab(x_label) + ylab(y_label) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", 
				size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}
