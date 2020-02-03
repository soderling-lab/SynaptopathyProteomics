#' moduleJS
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
#' moduleJS()()
moduleJS <- function() {
    x <- contrasts[idx,]
    idm1 <- x[["M1"]]
    idm2 <- x[["M2"]]
    m1 <- all_modules[[idm1]]
    m2 <- all_modules[[idm2]]
    if (idm1 == idm2) {
      modulejs <- 1
    } else {
      modulejs <- js(m1,m2)
    }
    return(modulejs)
}
