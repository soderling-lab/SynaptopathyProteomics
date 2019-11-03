#' TMT_exactTest
#'
#' A Function to evaluate DAPs using the exactTest from edgeR.
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
#' function_name(data_in, geno)
TMT_exactTest <- function(data_in, geno) {
  keep <- grep(geno, colnames(data_in))
  data_sub <- data_in[, keep]
  rownames(data_sub) <- rownames(data_in)

  # Insure that there are no rows with na, required for edgeR function
  data_sub <- na.omit(data_sub)

  # set up the sample mapping
  group <- c(rep("QC", 3), rep("WT", 4), rep("KO", 4))

  # make group into factors and set the order
  group <- factor(group, levels = c("QC", "WT", "KO"))

  # create a DGEList object with our data
  y_sl <- DGEList(counts = data_sub, group = group)
  y_sl <- estimateDisp(y_sl)

  # the exact test object has columns like fold-change, CPM, and p-values
  et_sl <- exactTest(y_sl, pair = c("WT", "KO"))

  # the topTags function adds the BH FDR values to an exactTest data frame. Make sure we do not change the row order!
  tt_sl <- topTags(et_sl, n = Inf, sort.by = "none")
  tt_sl <- tt_sl$table # tt_sl is a list. We just need the data frame table

  # add the default value as a new column
  tt_sl$candidate <- "no"
  tt_sl[which(tt_sl$FDR <= 0.10 & tt_sl$FDR > 0.05), dim(tt_sl)[2]] <- "low"
  tt_sl[which(tt_sl$FDR <= 0.05 & tt_sl$FDR > 0.01), dim(tt_sl)[2]] <- "med"
  tt_sl[which(tt_sl$FDR <= 0.01), dim(tt_sl)[2]] <- "high"
  tt_sl$candidate <- factor(tt_sl$candidate, levels = c("high", "med", "low", "no"))

  # Get Uniprot_IDs
  tt_sl <- add_column(tt_sl, rownames(tt_sl), .before = 1)
  Uniprot_IDs <- rownames(tt_sl)

  # Map Uniprot IDs to Gene names
  symbol <- mapIds(org.Mm.eg.db, keys = Uniprot_IDs, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
  tt_sl <- add_column(tt_sl, symbol, .after = 1)
  colnames(tt_sl)[1] <- "ID"
  rownames(tt_sl) <- NULL
  tt_sl$logCPM <- NULL
  tt_sl$Percent_WT <- 100 * (2.^tt_sl$logFC)
  colnames(tt_sl)[c(3:5, 7)] <- paste(geno, "_", colnames(tt_sl)[c(3:5, 7)], sep = "")
  return(tt_sl)
}
