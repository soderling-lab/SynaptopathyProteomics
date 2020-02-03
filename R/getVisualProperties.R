#' getVisualProperties
#'
#' description
#'
#' @param
#'
#' @return
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords
#'
#' @export
#'
#' @examples
#' getVisualProperties()()
getVisualProperties <- function() {
  visual.properties <- getVisualPropertyNames()
  graph.properties <- lapply(visual.properties, getVisualPropertyDefault)
  names(graph.properties) <- visual.properties
  if (!is.null(property)) {
    idx <- grep(toupper(property), names(graph.properties))
    return(graph.properties[idx])
  } else {
    return(graph.properties)
  }
}
