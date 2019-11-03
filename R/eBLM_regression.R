#' eBLM_regression
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
#' # A function for performing moderated empirical bayes regression for removal of covariate effects.
eBLM_regression <- function(data_in, traits, group, colID, cov, OLS = FALSE) {
  # Prepare the expression data.
  cols <- grep(group, colnames(data_in))
  data_sub <- data_in[, c(1, 2, cols)]
  out <- grepl("QC", colnames(data_sub))
  data_sub <- data_sub[, !out]
  traits_sub <- subset(traits, traits$ColumnName %in% colnames(data_sub))
  cols <- grep("Abundance", colnames(data_sub))
  data <- log2(as.matrix(data_sub[, cols]))
  rownames(data) <- data_sub$Accession
  data <- t(na.omit(data))
  data[1:5, 1:5]
  dim(data)
  dim(traits_sub)

  # Prepare the design df.
  sex <- as.factor(traits_sub$Sex)
  age <- as.numeric(traits_sub$Age)
  batch <- as.factor(traits_sub$PrepDate)
  status <- traits_sub$SampleType
  design <- as.data.frame(cbind(status, batch, sex, age))

  # Define covariates to be removed. Which model to choose?
  if (cov == 1) {
    covariates <- cbind(design$batch) # Batch
  } else if (cov == 2) {
    covariates <- cbind(design$batch, design$sex) # Batch + Sex
  } else if (cov == 3) {
    covariates <- cbind(design$sex) # Sex
  }

  covs <- c("Batch", "Batch+Sex", "Sex")[cov]

  # Correct for the batch effect using empiricalBayesLM from WGCNA package.
  fit.eblm <- empiricalBayesLM(data,
    removedCovariates = covariates,
    fitToSamples = design$status == "WT"
  )
  if (OLS == TRUE) {
    data.eblm <- fit.eblm$adjustedData.OLS
    print(paste("OLS used for regression of ", covs, ".", sep = ""))
  } else {
    data.eblm <- fit.eblm$adjustedData
    print(paste("eBLM used for regression of ", covs, ".", sep = ""))
  }

  # Check PCA after regression.
  plot1 <- ggplotPCA(t(data), traits_sub,
    colors = traits_sub$CustomColor,
    title = "2D PCA Plot (Pre-Regression)"
  )
  plot2 <- ggplotPCA(t(data.eblm), traits_sub,
    colors = traits_sub$CustomColor,
    title = "2D PCA Plot (Post-Regression)"
  )

  # Return regressed data, un-logged.
  data_out <- as.data.frame(2^t(data.eblm))

  # Add Peptides column.
  idx <- match(rownames(data_out), data_in$Accession)
  Peptides <- data_in$Peptides[idx]
  Accession <- rownames(data_out)

  data_out <- add_column(data_out, Accession, .before = 1)
  data_out <- add_column(data_out, Peptides, .before = 2)
  rownames(data_out) <- NULL

  results_list <- list(data_out, plot1, plot2, fit.eblm)
  names(results_list) <- c("adjustedData", "plot1", "plot2", "fit.eblm")
  return(results_list)
}
