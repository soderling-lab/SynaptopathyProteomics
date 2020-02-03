#' emptyVisualProperty
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
#' emptyVisualProperty()()
emptyVisualProperty <- function() {
# Empty visual property.
  visual.property <- list(
    visual.prop = NULL,
    table.column = NULL,
    mapping.type = NULL,
    table.column.values = NULL,
    visual.prop.values = NULL
  )
  return(visual.property)
}
