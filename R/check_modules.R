#' check_modules
#'
#' Check modules for evidence of preservation.
#'
#' @param selfPreservation result returned by NetRep.
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords none
#'
#' @examples
#' check_modules(selfPreservation)
check_modules <- function(selfPreservation,
                          strength = "strong",
                          stats = c(1:7),
                          alpha = 0.05) {

  # collect observed values, mean(nulls)
  spres <- selfPreservation
  obs <- spres$observed[, stats]

  # dim(spres$nulls) [1] 155 7 10000 # [nmodules, nstats, nperm]
  nulls <- apply(spres$nulls, 2, function(x) apply(x, 1, mean))[, stats]
  sd_null <- apply(spres$nulls, 2, function(x) apply(x, 1, sd))[, stats]

  # adjust pval for n module comparisons
  q <- apply(spres$p.values, 2, function(x) p.adjust(x, "bonferroni"))[, stats]
  q[is.na(q)] <- 1

  # eval(fx) ~ all() | weak()
  fx <- c("strong" = "all", "weak" = "any")[strength]

  if (length(stats) > 1) {
    # If testing more than one stat, evaluate 'strong' or 'weak' preservation:
    sig <- apply(q < alpha, 1, eval(fx))
    greater <- apply(obs > nulls, 1, eval(fx))
    less <- apply(obs < nulls, 1, eval(fx))
  } else {
    # If testing a single statistic don't use apply:
    sig <- q < alpha
    greater <- obs > nulls
    less <- obs < nulls
  }

  # Define preserved (greater), divergent (less), and ns modules.
  nModules <- length(spres$nVarsPresent)
  v <- rep("ns", nModules)

  v[greater & sig] <- "preserved"

  v[less & sig] <- "divergent"
  names(v) <- names(spres$nVarsPresent)

  # return partition with preservation enforced
  return(v)
}
