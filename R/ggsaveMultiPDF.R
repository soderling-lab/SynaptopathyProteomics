#' ggsaveMultiPDF
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
#' ## Define function: ggsave_plots(plot_list,file_ext)
ggsaveMultiPDF <- function(plot_list, savefigs_path, file_ext, verbose = TRUE) {
  # Saves plots in plot_list to file as given by savefigs_path.
  # Saves figures in format given by file_ext.
  # Will name plots as plot1, plot2, ect if not named in list.
  # If name exists in directory, then plot name will be plotn+1
  # Coerce to list.
  if (!inherits(plot_list, "list")) {
    plot_list <- list(plot_list)
  }
  # Loop to save plots.
  for (i in 1:length(plot_list)) {
    # Check if name exists.
    if (is.null(names(plot_list))) {
      save_name <- paste("plot", i, file_ext, sep = "")
      # If name exists, rename as plot(N+1)
      if (save_name %in% list.files(savefigs_path)) {
        idx <- grepl("plot", list.files(savefigs_path))
        num <- max(as.numeric(sapply(strsplit(list.files(savefigs_path)[idx], "plot|\\."), "[", 2))) + 1
        save_name <- paste("plot", num, file_ext, sep = "")
      }
    } else {
      save_name <- paste(names(plot_list)[i], file_ext, sep = "")
    }
    # Save the plots.
    suppressMessages(ggsave(save_name,
      plot = plot_list[[i]],
      scale = 1,
      width = NA,
      height = NA,
      path = savefigs_path
    ))
  }
  # Progress message:
  if (verbose == TRUE) {
    print(paste(length(plot_list), "plot(s) saved as", file_ext, "in", savefigs_path))
  }
}
