#' rip
#'
#' convienence wrapper around several R installation methods
#'
#' @param package - the package to be downloaded
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
#' rip("WGCNA")
rip <- function(package, method = "utils", ...) {
  # Install a R package. Supports the following methods:
  #     utils::install.packages()
  #     BiocManager::install()
  #     devtools::install_github()
  #     source - installs the package from Cran provided its source url,
  #              this method depends upon the Linux bash utility, rip..
  # If method is source, parse the package name from its url.
  if (method == "source") {
    url <- package
    package <- strsplit(strsplit(url, "/")[[1]][6], "_")[[1]][1]
  }
  # Insure that the package is not already installed.
  if (requireNamespace(package, quietly = TRUE)) {
    message(paste(package, "is already installed!"))
  } else if (method == "BiocManager") {
    BiocManager::install(package, ...)
  } else if (method == "utils") {
    utils::install.packages(package, ...)
  } else if (method == "devtools") {
    devtools::install_github(package, ...)
  } else if (method == "source") {
    cmd <- paste("rip", url, ...)
    system(cmd)
  } else {
    stop("problem installing package")
  }
}
