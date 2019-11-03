#' store_ggplot
#'
#' Saves ggplot(s) to list.
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
#' store_ggplot(plot_list, plot, name)
store_ggplot <- function(plot_list, plot, name) {
  k <- length(plot_list)
  plot_list[[k + 1]] <- plot
  names(plot_list)[k + 1] <- name
  return(plot_list)
}
