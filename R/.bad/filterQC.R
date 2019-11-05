#' filterQC
#'
#' Removes peptides that were not quantified in all three replicates.
#' Remove peptides from n x sd away from mean sd of given intensity bin.
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
#' filterQC(data_in, groups, nbins, threshold)
filterQC <- function(data_in, groups, nbins, threshold) {
  for (i in 1:length(groups)) {
    group <- groups[i]
    cols <- getCols(data_in, group)
    data_sub <- data_in[, cols]
    tmt_cols <- getCols(data_sub, "Abundance")
    qc_cols <- getCols(data_sub, "QC")

    # Remove peptide if QC was not quantified in all QC replicates.
    out <- is.na(data_sub[, qc_cols])
    logic <- matrix(NA, nrow = nrow(data_sub), ncol = 1)
    for (j in 1:nrow(data_sub)) {
      logic[j] <- any(out[j, ])
    }
    x <- length(logic[logic == TRUE])

    # Write to data_sub
    data_sub[logic, ] <- NA

    # Calcuate QC CV
    data_sub$QC_Average <- apply(data_sub[, qc_cols], 1, FUN = mean, na.rm = TRUE)
    data_sub$QC_SD <- apply(data_sub[, qc_cols], 1, FUN = sd, na.rm = TRUE)
    data_sub$QC_CV <- data_sub$QC_SD / data_sub$QC_Average * 100

    # Mean and SD
    mean_cv <- round(mean(data_sub$QC_CV, na.rm = TRUE), 2)
    sd_cv <- round(sd(data_sub$QC_CV, na.rm = TRUE), 2)

    # Print mean CV
    print(paste("The CV of", groups[i], "QC samples is", mean_cv, "\u00B1", sd_cv, "%", sep = " "))

    #  Filter based on 3 SDs away from mean SD of bins (5x)
    #  Get intensity bins
    rows_ignore <- is.nan(data_sub$QC_Average)
    data_sub$Intensity_bin[!rows_ignore] <- BurStMisc::ntile(data_sub$QC_Average[!rows_ignore], nbins,
      na.rm = TRUE,
      checkBleed = FALSE, result = "numeric"
    )

    # Calculate summary statistics of Intensity bins
    sd_info <- subset(data_sub) %>%
      group_by(Intensity_bin) %>%
      summarize(
        min = min(QC_Average, na.rm = TRUE),
        max = max(QC_Average, na.rm = TRUE),
        min_sd = min(QC_SD, na.rm = TRUE),
        max_sd = max(QC_SD, na.rm = TRUE),
        av_sd = mean(QC_SD, na.rm = TRUE),
        sd_sd = sd(QC_SD, na.rm = TRUE)
      )

    # Loop through peptide_data and determine if SD is outside threshold for intensity bin...
    logic <- matrix(NA, nrow = nrow(data_sub), ncol = 1)
    for (k in 1:nrow(data_sub)) {
      kbin <- data_sub$Intensity_bin[k]
      kmean <- sd_info$av_sd[kbin]
      ksd <- sd_info$sd_sd[kbin]
      logic[k] <- data_sub$QC_SD[k] > threshold * ksd
    }

    # Insure NA in logical vector are TRUE
    logic[is.na(logic)] <- TRUE

    # Print number of peptides that will be removed:
    y <- length(logic[logic == TRUE])
    print(paste(y - x, "peptides will be removed from", groups[i], "because of QC variability"))

    # write to data_filt
    data_filt <- data_sub[, grep("Abundance", colnames(data_sub))]
    data_filt[logic, ] <- NA
    cols <- grep(groups[i], colnames(data_in))
    data_in[, cols] <- data_filt
  }
  return(data_in)
}
