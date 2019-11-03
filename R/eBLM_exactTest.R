#' eBLM_exactTest
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
#' # Function for doing eBLM regression and exactTest.
eBLM_exactTest <- function(data_in, geno, sample_info, cov, OLS, col) {
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

  # Prepare the expression data.
  cols <- grep("Abundance", colnames(data_sub))
  data <- log2(as.matrix(data_sub[, cols]))
  rownames(data) <- data_sub$Accession
  data <- t(data)
  # data[1:5,1:5]
  # dim(data)

  # Prepare the design df.
  sex <- as.factor(traits$Sex)
  age <- as.numeric(traits$Age)
  batch <- as.factor(traits$PrepDate)
  status <- traits$SampleType
  design <- as.data.frame(cbind(status, batch, sex, age))
  covariates <- cbind(design$age, design$sex)
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
  if (OLS == TRUE) {
    data.eblm <- fit.eblm$adjustedData.OLS
    print("OLS used for regression.")
  } else {
    data.eblm <- fit.eblm$adjustedData
    print("eBLM used for regression.")
  }

  # Check PCA after regression.
  colors <- c(rep(col[1], 4), rep(col[2], 4))
  plot1 <- ggplotPCA(t(data), traits, colors, title = "2D PCA Plot (Pre-Regression)")
  plot2 <- ggplotPCA(t(data.eblm), traits, colors, title = "2D PCA Plot (Post-Regression)")

  # Check MDS after regression.
  colors <- c("gold1", "black")
  colID <- "Abundance"
  plot3 <- ggplotMDS(t(data), colID, colors, "MDS Pre-Regression", sample_info, labels = TRUE) + theme(legend.position = "none")
  plot4 <- ggplotMDS(t(data.eblm), colID, colors, "MDS Post-Regression", sample_info, labels = TRUE) + theme(legend.position = "none")

  ## EdgeR Exact Test
  data_EdgeR <- 2^t(data.eblm)

  # Create DGEList object with mapping to genotype (group).
  group <- factor(c(rep("WT", 4), rep("KO", 4)))
  y_DGE <- DGEList(counts = data_EdgeR, group = group)

  # Estimate dispersion.
  y_DGE <- estimateDisp(y_DGE)
  y_DGE <- estimateCommonDisp(y_DGE)
  y_DGE <- estimateTagwiseDisp(y_DGE)

  # Figure
  plot <- ggplotBCV(y_DGE)

  # Perform exactTest.
  y_ET <- exactTest(y_DGE, pair = c("WT", "KO"))

  # Map Uniprot IDs to Gene symbols and Entrez IDs.
  Uniprot_IDs <- rownames(y_ET)
  Entrez <- mapIds(org.Mm.eg.db,
    keys = Uniprot_IDs,
    column = "ENTREZID", keytype = "UNIPROT", multiVals = "first"
  )
  # Calculate GO And KEGG enrichment results.
  go <- goana(y_ET, geneid = Entrez, species = "Mm", FDR = 0.1)
  go <- go[go$P.Up < 0.05 | go$P.Down < 0.05, ]
  keg <- kegga(y_ET, geneid = Entrez, species = "Mm", FDR = 0.1)
  keg <- keg[keg$P.Up < 0.05 | keg$P.Down < 0.05, ]

  # Call topTags to add FDR. Keep the data the same order by sort.by="none".
  y_TT <- topTags(y_ET, n = Inf, sort.by = "none")

  # Extract the results from the result.
  y_TT <- y_TT$table

  # Categorize candidates by FDR.
  y_TT$candidate <- "no"
  y_TT[which(y_TT$FDR <= 0.10 & y_TT$FDR > 0.05), dim(y_TT)[2]] <- "low"
  y_TT[which(y_TT$FDR <= 0.05 & y_TT$FDR > 0.01), dim(y_TT)[2]] <- "med"
  y_TT[which(y_TT$FDR <= 0.01), dim(y_TT)[2]] <- "high"
  y_TT$candidate <- factor(y_TT$candidate, levels = c("high", "med", "low", "no"))

  # Map Uniprot IDs to Gene symbols and Entrez IDs.
  y_TT <- add_column(y_TT, rownames(y_TT), .before = 1)
  Uniprot_IDs <- rownames(y_TT)
  Gene <- mapIds(org.Mm.eg.db, keys = Uniprot_IDs, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
  Entrez <- mapIds(org.Mm.eg.db, keys = Uniprot_IDs, column = "ENTREZID", keytype = "UNIPROT", multiVals = "first")
  y_TT <- add_column(y_TT, Gene, .after = 1)
  y_TT <- add_column(y_TT, Entrez, .after = 2)
  colnames(y_TT)[1] <- "Uniprot"
  rownames(y_TT) <- NULL

  # nsig results
  nsig <- y_TT[, 7] < 0.1
  nsig <- length(nsig[nsig == TRUE])

  # Replace CPM column with percent WT (2^log2FC).
  y_TT$logCPM <- 100 * (2.^y_TT$logFC)
  colnames(y_TT)[5] <- "%WT"
  colnames(y_TT)[c(4:7)] <- paste(geno, "_", colnames(y_TT)[c(4:7)], sep = "")

  results_list <- list(fit.eblm, y_DGE, y_TT, go, keg, nsig, filt_protein, TMM_protein, plot1, plot2, plot3, plot4)
  names(results_list) <- c(
    "fit.eblm", "y_DGE", "y_TT", "GO", "KEGG", "nsig",
    "filt_protein", "TMM_protein", "plot1", "plot2", "plot3", "plot4"
  )
  return(results_list)
}
