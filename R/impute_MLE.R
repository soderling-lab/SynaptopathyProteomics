#' impute_MLE
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
#' function_name(param1, ... )
#-------------------------------------------------------------------------------
#' ## impute_MLE
#' Updated 021419. Previosly data was returned with no imputed values!!!
# This function replaces missing values using MLE algorithm for MAR data.
# qc_threshold - number of tolerated missing qc values within an experiment
# bio_threshold - number of tolerated missing biological replicate valutes within an exp
# group - column identifier for experiments

impute_MLE <- function(data_in, groups, qc_threshold = 0, bio_threshold = 2) {
  for (i in 1:length(groups)) {
    data_work <- data_in
    rownames(data_work) <- paste(data_in$Accession, c(1:nrow(data_in)), sep = "_")
    group <- groups[i]
    cols <- grep(group, colnames(data_work))
    data_sub <- data_work[, cols]
    data_sub$Sort <- c(1:nrow(data_sub))
    # Ignore rows with too many missing values.
    data_sub$qc_NA <- apply(data_sub[, 1:3], 1, function(x) sum(is.na(x)))
    unique(data_sub$qc_NA)
    no_impute_qc <- data_sub$qc_NA > qc_threshold
    data_sub$bio1_NA <- apply(data_sub[, 4:7], 1, function(x) sum(is.na(x)))
    unique(data_sub$bio1_NA)
    data_sub$bio2_NA <- apply(data_sub[, 8:11], 1, function(x) sum(is.na(x)))
    unique(data_sub$bio2_NA)
    no_impute_bio <- data_sub$bio1_NA > bio_threshold | data_sub$bio2 > bio_threshold
    rows_out <- no_impute_bio == TRUE | no_impute_qc == TRUE
    length(rows_out[rows_out == TRUE])
    # Replace values in rows to ignore with NA.
    data_out <- data_sub[rows_out, ]
    data_out[, c(1:11)] <- NA
    # Data for imputing
    data_temp <- data_sub[!rows_out, ]
    data_temp[, c(1:11)] <- log2(data_temp[, c(1:11)])
    # Number of missing values
    num_NA <- is.na(data_temp[, c(1:11)])
    num_NA <- length(num_NA[num_NA == TRUE])
    cat(paste(num_NA, "values from", groups[i], "are missing and will be replaced by imputing.\n", sep = " "))
    # Try MLE impute
    conditions <- as.factor(c(rep(1, 3), rep(2, 4), rep(3, 4)))
    # MLE Impute
    data_imp <- data_temp
    data_imp[, c(1:11)] <- impute.mle(data_temp[, c(1:11)], conditions)
    # Put back together with data_out
    data <- rbind(data_out, 2^data_imp) # Un-log.
    # Sort to original order.
    order <- as.numeric(do.call(rbind, strsplit(rownames(data), "_"))[, 2])
    data$Order <- order
    data <- data[order(data$Order), ]
    # Output
    data_return <- data_sub[, c(1:11)]
    data_return[, c(1:11)] <- data[, c(1:11)]
    data_in[, cols] <- data_return[, c(1:11)]
  }
  return(data_in)
}

