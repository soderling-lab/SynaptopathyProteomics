#' eBLM_glmQLFit_QLFTest
#'
#' A Function for performing eBLM regression and EdgeR GLM and QLFTest.
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
#' eBLM_glmQLFit_QLFTest(data_in, geno, traits, cov, col, contrast)
eBLM_glmQLFit_QLFTest <- function(data_in, geno, traits, cov, col, contrast) {

  # Subset data
  cols <- grep(geno, colnames(data_in))
  data_sub <- data_in[, c(1, 2, cols)]

  # Filter proteins.
  filt_protein <- filter_proteins(data_sub, "Abundance")

  # TMM Normalization
  TMM_protein <- normalize_TMM(filt_protein, "Abundance")

  # Remove QC from expression data and traits.
  out <- grepl("QC", colnames(TMM_protein))
  data_sub <- TMM_protein[, !out]
  traits <- sample_info[grepl(geno, sample_info$Model), ]
  traits <- subset(traits, SampleType != "QC")
  traits <- traits[order(traits$Order), ]
  dim(traits)

  # Prepare the expression data.
  cols <- grep("Abundance", colnames(data_sub))
  data <- log2(as.matrix(data_sub[, cols]))
  rownames(data) <- data_sub$Accession
  data <- t(data)
  data[1:5, 1:5]

  # Prepare the design df.
  sex <- as.factor(traits$Sex)
  age <- as.numeric(traits$Age)
  batch <- as.factor(traits$PrepDate)
  status <- traits$SampleType
  design <- as.data.frame(cbind(status, batch, sex, age))
  design

  # Define covariates to be removed. Which model to choose?
  if (cov == 1) {
    covariates <- cbind(design$batch) # Batch
  } else if (cov == 2) {
    covariates <- cbind(design$batch, design$sex) # Batch + Sex
  } else if (cov == 3) {
    covariates <- cbind(design$sex) # Sex
  }

  # Correct for the batch effect using empiricalBayesLM from WGCNA package.
  fit.eblm <- empiricalBayesLM(data,
    removedCovariates = covariates,
    fitToSamples = design$status == "WT"
  )
  data.eblm <- fit.eblm$adjustedData

  # Check PCA after regression.
  colors <- c(rep("gray", 4), rep("black", 4))
  plot1 <- ggplotPCA(t(data), traits, colors, title = "2D PCA Plot (Pre-Regression)")
  plot2 <- ggplotPCA(t(data.eblm), traits, colors, title = "2D PCA Plot (Post-Regression)")
  # plot_grid(plot1,plot2)

  # The counts data
  cleanDat <- 2^t(data.eblm)
  # cleanDat[1:5,1:5]

  # Prepare traits matrix
  rows <- grepl(geno, traits$Model)
  traits_sub <- traits[rows, ]

  if (!all(traits_sub$ColumnName == colnames(cleanDat))) {
    print("Warning:design does not match data!")
  }

  # Create DGEList object.
  y_DGE <- DGEList(counts = cleanDat)

  # Create sample mapping.
  group <- traits_sub$SampleType
  sex <- as.factor(traits_sub$Sex)
  age <- as.numeric(traits_sub$Age)
  batch <- as.factor(traits_sub$PrepDate)

  # Add sample mapping to DGE object.
  y_DGE$samples$group <- as.factor(group)

  # Create Design matrix (model.matrix(~batch+group)).
  design <- model.matrix(~ 0 + group, data = y_DGE$samples)
  colnames(design)[c(1, 2)] <- levels(y_DGE$samples$group)
  # design

  # Estimate dispersion:
  y_DGE <- estimateDisp(y_DGE)

  # Plot BCV.
  plotBCV(y_DGE)

  # Fit a general linear model.
  fit <- glmQLFit(y_DGE, design, robust = TRUE)

  # Create a list of contrasts
  if (contrast == 1) {
    contrasts <- list(WTvKO <- makeContrasts(KO - WT, levels = design))
  } else {
    contrasts <- list(WTvHET <- makeContrasts(HET - WT, levels = design))
  }

  # Call function results_QLFTest to perform QLFtest for desired contrasts.
  res1 <- results_QLFTest(fit, comparison = contrasts[[1]], alpha = 0.1)

  # Pvalue histogram
  plot3 <- ggplotPvalHist(res1$QLFTest_results, color = col, geno)

  results_list <- list(res1, plot1, plot2, plot3)
  names(results_list) <- c("QLFTest_Results", "plot1", "plot2", "plot3")
  return(results_list)
}
