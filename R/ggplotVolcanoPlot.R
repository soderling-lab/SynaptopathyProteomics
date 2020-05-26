# Function for producing volcano plots.
ggplotVolcanoPlot <- function(df) {
  df$x <- df[, grep("FC", colnames(df))]
  df$y <- -log10(df[, grep("PValue", colnames(df))])
  logic <- df$FDR < 0.05
  df$Color[!logic] <- "gray"
  df$Color <- as.factor(df$Color)
  y_int <- -1 * log10(max(df$PValue[df$FDR < 0.05]))
  # Generate plot.
  plot <- ggplot(data = df, aes(x = x, y = y, color = Color)) +
    geom_point(size = 3, alpha = 0.5) +
    scale_color_manual(values = levels(df$Color)) +
    geom_hline(
      yintercept = y_int, linetype = "dashed",
      color = "black", size = 0.6
    ) +
    geom_vline(
      xintercept = 0, linetype = "dashed",
      color = "black", size = 0.6
    ) +
    xlab(expression(bold(Log[2](Fold ~ Change)))) +
    ylab(expression(bold(-Log[10](P - value)))) +
    theme(
      plot.title = element_text(
        hjust = 0.5, color = "black",
        size = 12, face = "bold"
      ),
      axis.title.y = element_text(
        color = "black", face = "bold",
        size = 11, angle = 90, vjust = 0.5
      ),
      axis.title.x = element_text(
        color = "black", face = "bold",
        size = 11, angle = 0,
        hjust = 0.5, vjust = 0.5
      ),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = "none"
    )
  # Add annotation.
  ypos <- unlist(ggplot_build(plot)$layout$panel_params[[1]][8])
  xpos <- range(plot$data$x, na.rm = TRUE)
  plot <- plot + annotate("text",
    x = xpos[1] + 0.3 * (xpos[2] - xpos[1]),
    y = y_int + 0.04 * (ypos[2] - ypos[1]),
    label = "FDR < 0.05", size = 4
  )
  return(plot)
}
