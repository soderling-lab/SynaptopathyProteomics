#' store_ggplot
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
#' ## Define function: store_plot
#' #  Saves ggplot(s) to list.
store_ggplot <- function(plot_list, plot, name) {
  k <- length(plot_list)
  plot_list[[k + 1]] <- plot
  names(plot_list)[k + 1] <- name
  return(plot_list)
}
