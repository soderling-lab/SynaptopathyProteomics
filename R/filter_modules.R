#' filter_modules
#'
#' remove small modules from partition
#'
#' @param partition - named vector 
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
#'filter_modules(partition,cutoff=5)
filter_modules <- function(partition, cutoff = 1) {
modules <- split(partition, partition)
out <- names(modules)[sapply(modules, function(x) table(x) < cutoff)]
partition[partition %in% out] <- 0
return(partition)
}
