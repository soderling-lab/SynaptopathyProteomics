#' description: functions for TMT analysis
#' authors: Tyler W Bradshaw
#' ---

################################################################################
# Note:
# Please note,  these functions are not robust. They do not typically check the
# inputs they are provided. They were made to operate on inputs specific to this
# analysis. Their general use should be applied with caution.
################################################################################

#' fill down
#'
#' Fill in a dataframe with missing values.
#'   Missing values are replaced with the value above them.
#'   From StackOverflow user [nacnudus](https://stackoverflow.com/users/937932/nacnudus).
#'
#' @param x column vector with blank values.
#' @param blank logic vector specifying blank values.
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{https://stackoverflow.com/questions/10554741/fill-in-data-frame-with-values-from-rows-above}
#' @keywords fill down blank missing values
#'
#' @examples
#' fill_down()
#' @export
# @importFrom grDevices rgb2hsv
#'
fill_down <- function(x, blank = is.na) {
  # Find the values
  if (is.function(blank)) {
    isnotblank <- !blank(x)
  } else {
    isnotblank <- x != blank
  }
  # Fill down
  x[which(isnotblank)][cumsum(isnotblank)]
}

#-------------------------------------------------------------------------------
#' cleanPD
#'
#' Clean up raw TMT data exported as xlsx from PD.
#'   Utilizes the fill_down function.
#'
#' @param data_in data frame imported from read_excel.
#' @param sample_info data frame with sample information.
#'
#' @return reformated data frame with raw peptide data from PD.
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords PD ProteomeDiscover raw intensity peptides
#'
#' @examples
#' cleanPD()
#' @export
# @importFrom
#'
cleanPD <- function(data_in, sample_info) {
  logic <- is.na(data_in$`Protein FDR Confidence: Mascot`)
  data_in[logic, c(1, 3, 4)] <- NA
  data_in[, c(1, 3, 4)] <- apply(data_in[, c(1, 3, 4)], 2, function(x) fill_down(x))
  data_in <- subset(data_in, Master == "High")
  colnames(data_in) <- colnames(data_in)[-6]
  colnames(data_in)[ncol(data_in)] <- "Quan Usage"
  data_in <- subset(data_in, data_in$`Quan Usage` == "Used")
  cols_out <- dim(data_in)[2] - length(grep("Abundance", colnames(data_in))) - 6
  data_in <- data_in[, c(1:(ncol(data_in) - cols_out))]
  tmt_cols <- grep("Abundance", colnames(data_in))
  data_in[, tmt_cols] <- sapply(data_in[, tmt_cols], as.numeric)
  data_sort <- data_in[, tmt_cols]
  target <- sample_info$ColumnName[order(sample_info$Order)]
  data_sort <- data_sort[, match(target, colnames(data_sort))]
  data_out <- cbind(data_in[, 1:6], data_sort)
  data_out <- data_out[, -1]
  colnames(data_out)[c(1, 4, 5)] <- c("Confidence", "Sequence", "Modifications")
  return(data_out)
}

#-------------------------------------------------------------------------------
#' getCols
#'
# Gets a vector of column numbers for a given ID. An alternative to grep().
#'
#' @param data_in a data frame.
#' @param ID character specifying column of interest.
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords grep getCols index data frame column name
#'
#' @examples
#' cleanPD()
#' @export
# @importFrom
#'
getCols <- function(data_in, ID) {
  cols <- grep(ID, colnames(data_in))
  return(cols)
}

#-------------------------------------------------------------------------------
# ## peptide_overlap
#-------------------------------------------------------------------------------
# Calculates peptide identification overlap for given contrasts matrix.

peptide_overlap_TMT <- function(data_in, contrasts, info_cols) {
  num_iter <- dim(contrasts)[2]
  contrasts <- rbind(contrasts, matrix(NA, nrow = 3, ncol = num_iter))
  for (i in 1:num_iter) {
    IDa <- contrasts[1, i]
    tmt_cols <- grep(IDa, colnames(data_in))
    data_sub <- na.omit(data_in[, c(info_cols, tmt_cols)])
    listA <- unique(data_sub$Sequence)

    IDb <- contrasts[2, i]
    tmt_cols <- grep(IDb, colnames(data_in))
    data_sub <- na.omit(data_in[, c(info_cols, tmt_cols)])
    listB <- unique(data_sub$Sequence)

    overlap <- length(unique(intersect(listA, listB)))
    total <- length(unique(union(listA, listB)))
    percent <- 100 * (overlap / total)

    contrasts[3, i] <- overlap
    contrasts[4, i] <- total
    contrasts[5, i] <- round(percent, 4)
  }
  return(contrasts)
}

#-------------------------------------------------------------------------------
#' ## corQC
#-------------------------------------------------------------------------------
# GGplots are saved as objects in a list.
# Scatter plot showing correlation of QC samples. All comparisons.
# nbin number of histograms showind distribution of QC ratios for all intensity bins.

ggplotcorQC <- function(data_in, groups, colID, nbins, annotate = TRUE) {
  plot_list <- list()
  for (i in 1:length(groups)) {
    cols <- grep(groups[i], colnames(data_in))
    data_sub <- data_in[, cols]
    QCcols <- grep(colID, colnames(data_sub))
    data_work <- na.omit(data_sub[, QCcols])
    contrasts <- combn(c(1:ncol(data_work)), 2) # QC comparisons.
    num_iter <- dim(contrasts)[2]
    data_list <- list()
    for (k in 1:num_iter) {
      x <- contrasts[1, k]
      y <- contrasts[2, k]
      Log2QC1 <- log2(data_work[, x])
      Log2QC2 <- log2(data_work[, y])
      Ratio <- Log2QC1 - Log2QC2
      data_list[[k]] <- cbind(Log2QC1, Log2QC2, Ratio)
    }
    # merge data frames in list.
    data <- as.data.frame(do.call(rbind, data_list))
    # Bin by mean intensity.
    mu <- rowMeans(data_work)
    data$bins <- rep(BurStMisc::ntile(mu, nbins, na.rm = TRUE, checkBleed = FALSE, result = "numeric"), num_iter)
    # Determine best fit line.
    fit <- lm(data$Log2QC1 ~ data$Log2QC2)
    # Calculate Pearson P-Value.
    corTest <- cor.test(~ data$Log2QC1 + data$Log2QC2,
      data = cbind(data$Log2QC1, data$Log2QC2),
      method = "pearson", conf.level = 0.95
    )
    Slope <- paste("Slope =", round(coef(fit)[2], 4))
    R2 <- paste("R2 =", round(corTest$estimate, 4))
    # Generate scatter plot.
    plot <- ggplot(data, aes(x = Log2QC1, y = Log2QC2, color = bins)) + geom_point() +
      scale_color_continuous(name = "Intensity Bin") +
      # geom_abline(intercept=coef(fit)[1],slope=coef(fit)[2], color = "black", linetype = "dashed") +
      ggtitle(groups[i]) +
      xlab(expression(Log[2] ~ QC1)) +
      ylab(expression(Log[2] ~ QC2)) +
      theme(
        plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
        axis.title.x = element_text(color = "black", size = 11, face = "bold"),
        axis.title.y = element_text(color = "black", size = 11, face = "bold")
      )
    # Add annotation layer.
    mytable <- rbind(R2, Slope)
    xrange <- unlist(ggplot_build(plot)$layout$panel_params[[1]][1])
    yrange <- unlist(ggplot_build(plot)$layout$panel_params[[1]][8])
    xmin <- min(xrange)
    xmax <- max(xrange)
    xdelta <- xmax - xmin
    ymin <- min(yrange)
    ymax <- max(yrange)
    ydelta <- ymax - ymin
    tt <- ttheme_default(base_size = 11, core = list(bg_params = list(fill = "white")))
    tab <- tableGrob(mytable, rows = NULL, theme = tt)
    g <- gtable_add_grob(tab,
      grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
      t = 1, b = nrow(tab), l = 1, r = ncol(tab)
    )
    if (annotate == TRUE) {
      plot <- plot + annotation_custom(g,
        xmin = xmin - 0.65 * xdelta, xmax,
        ymin = ymin + 0.8 * ydelta, ymax
      )
    }
    plot_list[[i]] <- plot
    names(plot_list)[[i]] <- groups[[i]]
  }
  return(plot_list)
}

#-------------------------------------------------------------------------------
# ## df2dm_TMT
#-------------------------------------------------------------------------------
# This is a simple function to convert data frame to dm for plotting purposes.

df2dm_TMT <- function(df, colID) {
  tmt_cols <- grep(colID, colnames(df))
  dm <- as.matrix(df[, tmt_cols])
  rownames(dm) <- if ("Accession" %in% df) {
    df$Accession
  } else {
    rownames(df)
  } #*** EBD edited
  colnames(dm) <- c(1:ncol(dm))
  return(dm)
}

#-------------------------------------------------------------------------------
# ## calc_CV_TMT
#-------------------------------------------------------------------------------
# A function to calculate the Coefficient of variation (STDEV/Mean) of log2
# transformed input.
# Groups defined by groups
# Columns with numerical data defined by ColID.
# row identifier should be provided by "Accession" column

calc_CV_TMT <- function(data_in, groups, colID) {
  data_work <- df2dm_TMT(data_in, colID)
  colnames(data_work) <- colnames(data_in)[grep(colID, colnames(data_in))]
  names <- as.vector(data_in$Accession)
  nums <- as.numeric(c(1:nrow(data_in)))
  rownames(data_work) <- paste(names, nums, sep = "_")
  data_cv <- matrix(NA, nrow = dim(data_work)[1], ncol = length(groups))
  for (i in 1:length(groups)) {
    id <- grep(groups[i], colnames(data_work))
    mu <- apply(data_work[, id], 1, mean, na.rm = TRUE)
    stdev <- apply(data_work[, id], 1, sd, na.rm = TRUE)
    data_cv[, i] <- 100 * (stdev / mu)
  }
  colnames(data_cv) <- paste(groups, "CV", sep = "_")
  rownames(data_cv) <- rownames(data_work)
  return(data_cv)
}

#-------------------------------------------------------------------------------
#' ## impute_MLE
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
#' ## normalize_SL
#-------------------------------------------------------------------------------
# Performs sample loading normalization based on column sums of TMT-runs.
# This function also removes rows in which qc was not quantified in all three replicates.

normalize_SL <- function(data_in, colID = "", group = "") {

  # Seperate data columns and info columns
  tmt_cols <- grep(colID, colnames(data_in))
  info_cols <- data_in[, c(1:ncol(data_in))[-tmt_cols]]
  data_cols <- data_in[, tmt_cols]
  # Insure 0 is NA
  data_cols[data_cols == 0] <- NA
  # Loop through groups, calculate SL norm
  for (i in 1:length(group)) {
    group_cols <- grep(group[i], colnames(data_cols))
    sub_data <- data_cols[, group_cols]
    target <- mean(colSums(sub_data, na.rm = TRUE), na.rm = TRUE)
    norm_facs <- target / colSums(sub_data, na.rm = TRUE)
    sl_data <- sweep(sub_data, 2, norm_facs, FUN = "*")
    data_cols[, group_cols] <- sl_data
  }
  # Bind info columns and data columns
  data_out <- cbind(info_cols, data_cols)
  return(data_out)
}

#-------------------------------------------------------------------------------
# ## filterQC
#-------------------------------------------------------------------------------
# Removes peptides that were not quantified in all three replicates.
# Remove peptides from n x sd away from mean sd of given intensity bin.

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

#-------------------------------------------------------------------------------
#' ## filterQCv2
#-------------------------------------------------------------------------------
# This function bins peptides by mean intensity of QC samples, calculates the
# log ratio of all QC sample comparisons and then removes a peptides if its mean QC ratio
# is outside NxSD the mean its intensity bin.

filterQCv2 <- function(data_in, groups, nbins, threshold) {
  for (i in 1:length(groups)) {
    group <- groups[i]
    # Get the data for the specified group.
    cols <- grep(group, colnames(data_in))
    data_sub <- data_in[, cols]
    rownames(data_sub) <- paste(data_in$Accession, data_in$Sequence, c(1:nrow(data_in)), sep = "_")
    tmt_cols <- grep("Abundance", colnames(data_sub))
    qc_cols <- grep("QC", colnames(data_sub))

    # Remove row if QC was not quantified in all QC replicates.
    out <- is.na(data_sub[, qc_cols])
    logic_out <- matrix(NA, nrow = nrow(data_sub), ncol = 1)
    for (j in 1:nrow(data_sub)) {
      logic_out[j] <- any(out[j, ])
    }

    # Write to data_sub
    data_sub[logic_out, ] <- NA

    # Log2 transform and define tmt and qc columns
    data_work <- log2(data_sub)
    cols_qc <- grep("QC", colnames(data_work))
    cols_tmt <- grep("Abundance", colnames(data_work))

    # Calculate the ratio of QC samples (QC1-QC2, QC1-QC3, QC2-QC3)
    contrasts <- combn(cols_qc, 2)
    ratio_dm <- NULL
    for (m in 1:ncol(contrasts)) {
      idx <- contrasts[1, m]
      idy <- contrasts[2, m]
      ratio_dm <- cbind(ratio_dm, data_work[, idx] - data_work[, idy])
    }

    colnames(ratio_dm) <- paste("ratio_", c(1:ncol(ratio_dm)), sep = "")
    data_work <- cbind(data_work, ratio_dm)
    cols <- grep("ratio_", colnames(data_work))
    if (length(cols) == 1) {
      data_work$ratio_avg <- data_work[, cols]
    } else {
      data_work$ratio_avg <- rowMeans(data_work[, cols], 1, na.rm = TRUE)
    }
    qc_cols <- grep("QC", colnames(data_work))
    data_work$avgQC <- rowMeans(data_work[, qc_cols], 1, na.rm = TRUE)
    data_work$sdQC <- apply(data_work[, qc_cols], 1, FUN = sd, na.rm = TRUE)

    # Calculate bins based on mean intensity of QC replicates.
    rows_ignore <- is.nan(data_work$avgQC)
    data_work$bins[!rows_ignore] <- BurStMisc::ntile(data_work$avgQC[!rows_ignore], nbins,
      na.rm = TRUE,
      checkBleed = FALSE, result = "numeric"
    )

    sdQC <- aggregate(data_work$ratio_avg, by = list(bin = data_work$bins), FUN = sd)
    avgQC <- aggregate(data_work$ratio_avg, by = list(bin = data_work$bins), FUN = mean)

    # Loop through bins and determine if mean ratio is outside N*SD away from mean precision.
    logic <- matrix(NA, nrow = nrow(data_work), ncol = 1)
    for (k in 1:nrow(data_work)) {
      kbin <- data_work$bins[k]
      kmean <- avgQC$x[kbin]
      ksd <- sdQC$x[kbin]
      logic[k] <- data_work$ratio_avg[k] > kmean + threshold * ksd || data_work$ratio_avg[k] < kmean - 1 * threshold * ksd
    }

    # Print number of peptides that will be removed:
    num_filt <- length(logic[logic == TRUE]) - length(logic[is.na(logic)])
    print(paste(num_filt, "peptides will be removed from ", group, " because of QC imprecision"))

    # Insure NA in logical vector are TRUE
    logic[is.na(logic)] <- TRUE

    # Add logic column.
    data_work$out <- logic

    # write to data_filt (unlog)
    data_filt <- 2.^data_work[, tmt_cols]
    data_filt[logic, ] <- NA

    # write to data_in
    cols <- grep(group, colnames(data_in))
    data_in[, cols] <- data_filt
  }
  return(data_in)
}

#-------------------------------------------------------------------------------
#' ## summarize_Protein
#-------------------------------------------------------------------------------
# Summarize proteins by summing all peptides for a unique Accession identifier.

summarize_Protein <- function(peptide_data) {
  # Add column for peptides, summarize using dplyr::summarize_all(sum)
  Peptides <- rep(1, nrow(peptide_data))
  temp_data <- add_column(peptide_data, Peptides, .after = 5)
  tmt_cols <- getCols(temp_data, "Abundance")
  temp_data <- temp_data[, c(2, 6, tmt_cols)]
  prot_data <- temp_data %>%
    group_by(Accession) %>%
    summarise_all(funs(sum), na.rm = TRUE)
  # Replace 0 with NA
  prot_data[prot_data == 0] <- NA
  return(prot_data)
}

#-------------------------------------------------------------------------------
# ## formatDEP
#-------------------------------------------------------------------------------

# Reformats the data for input into DEP package

formatDEP <- function(data_in) {
  #  Add ID and names columns
  data_in$ID <- data_in$Accession
  data_in$name <- data_in$Accession
  ## Create SE object
  TMT_columns <- grep("F", colnames(data_in)) # Intensity column numbers
  long.cols <- colnames(data_in)[TMT_columns]
  long.cols <- gsub("\\ ", "", long.cols)
  long.cols <- gsub("\\..", "", long.cols)
  vars <- colsplit(long.cols, ":", c("A", "B", "C"))
  vars <- vars[, 3]
  vars <- colsplit(vars, ",", c("channel", "sample", "condition", "proteomicsID", "genotype")) # Cortex format
  vars$sample <- NULL
  vars$replicate <- rep(c(1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4), 4)
  label <- paste(vars$genotype, vars$condition, ".", vars$replicate, sep = "")
  # Create SE object
  colnames(data_in)[TMT_columns] <- label
  data_se <- make_se_parse(data_in, TMT_columns)
  return(data_se)
}

#-------------------------------------------------------------------------------
# ## TMT_exactTest
#-------------------------------------------------------------------------------
# A Function to evaluate DEP using the exactTest from edgeR.

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

#-------------------------------------------------------------------------------
#' ## normalize_IRS
#-------------------------------------------------------------------------------
# Function for IRS normalization.
# Supports geometric or robust mean.

normalize_IRS <- function(data_in, IRS_ID, groups, robust) {

  # Subset the data.
  cols_qc <- grep(IRS_ID, colnames(data_in))
  data_qc <- data_in[, cols_qc]

  # Empty list for output of loop.
  data_list <- list()

  # Loop to calculate rowMeans of QC channels within an experiment.
  for (i in 1:length(groups)) {
    exp_cols <- grep(groups[i], colnames(data_qc))
    df <- as.data.frame(rowMeans(data_qc[, exp_cols], na.rm = TRUE, dims = 1))
    colnames(df) <- paste(groups[1], "QC", "rowMeans", sep = "_")
    data_list[[i]] <- df
  }

  # Bind the data frames in list.
  df_IRS <- do.call(cbind, data_list)

  # Calculate mean of QC means (supports Robust, geometric mean).
  if (robust == TRUE) {
    # Calculate geometric rowMeans, ignore missing values.
    df_IRS$Avg <- apply(df_IRS, 1, function(x) exp(mean(log(x), na.rm = TRUE)))
    print("Used robust (geometric) mean.")
  } else {
    df_IRS$Avg <- apply(df_IRS, 1, function(x) mean(x, na.rm = TRUE))
    print("Used arithmetic mean.")
  }

  # Compute scaling factors for each experiment.
  factors_list <- list()
  for (j in 1:length(groups)) {
    factors_list[[j]] <- df_IRS$Avg / df_IRS[, j]
  }

  # Loop through factors list and generate matrix of factors, store in a new list.
  new_list <- list()
  for (k in 1:length(groups)) {
    num_channels <- length(grep(groups[k], colnames(data_in)))
    new_list[[k]] <- matrix(factors_list[[k]], nrow = nrow(data_in), ncol = num_channels)
  }

  # Bind factors in list.
  dm_factors <- do.call(cbind, new_list)

  # Perform IRS (factors * data_in)
  tmt_cols <- grep("Abundance", colnames(data_in))
  data_IRS <- data_in
  data_IRS[, tmt_cols] <- dm_factors * data_in[, tmt_cols]

  # Return IRS data.
  return(data_IRS)
}

#-------------------------------------------------------------------------------
# ## ttest_TMT
#-------------------------------------------------------------------------------
# This function just does all the work of ttest.
ttest_TMT <- function(data_in, groups) {
  data_temp <- as.data.frame(assay(data_in))
  data_out <- matrix(NA, nrow = nrow(data_temp), ncol = length(groups))
  data_out <- as.data.frame(data_out)
  rownames(data_out) <- rownames(data_temp)
  for (i in 1:length(groups)) {
    data_sub <- data_temp[, grep(groups[i], colnames(data_temp))]
    results_ttest <- apply(data_sub, 1, function(x) t.test(x[5:8], x[9:12]))
    ttest_pvalue <- unlist(lapply(results_ttest, function(x) x$p.value))
    ttest_fdr <- p.adjust(ttest_pvalue, method = "BH")
    data_out[, i] <- ttest_pvalue
    colnames(data_out)[i] <- paste(groups[i], "p.value", sep = "_")
    data_out <- add_column(data_out, ttest_fdr)
    colnames(data_out)[ncol(data_out)] <- paste(groups[i], "fdr", sep = "_")
  }
  return(data_out)
}

#-------------------------------------------------------------------------------
# ## Define function: clean_DEP
# This function just compresses the work for cleaning up the DEP limma results.

clean_DEP <- function(data_in, alpha, lfc) {
  # Denote significant proteins: no logFC cut-off, alpha = 0.1
  dep <- add_rejections(data_diff, alpha = alpha, lfc = lfc)

  # Generate a results table
  results_limma <- get_results(dep)

  # Drop un-wanted columns
  keep <- c(1, grep("p.val", colnames(results_limma)), grep("p.adj", colnames(results_limma)))
  results_limma <- results_limma[, keep]

  # Calculate Log2FC from data_in
  data_temp <- assay(data_in)
  log2fc <- matrix(NA, nrow(data_temp), 4)
  group <- c("Shank2", "Syngap1", "Ube3a", "Shank3")

  # Loop to calc log2fc (Avg.KO-Avg.WT)
  for (i in 1:4) {
    sub <- grep(group[i], colnames(data_temp))
    data_sub <- data_temp[, sub]
    log2fc[, i] <- rowMeans(data_sub[, 5:8]) - rowMeans(data_sub[, 1:4])
  }

  # Clean up result
  log2fc <- as.data.frame(log2fc)
  colnames(log2fc) <- c("Shank2_Log2FC", "Syngap1_Log2FC", "Ube3a_Log2FC", "Shank3_Log2FC")
  rownames(log2fc) <- rownames(data_temp)
  order <- match(results_limma$name, rownames(log2fc))
  log2fc <- log2fc[order, ]

  # Combine with limma results
  results_limma <- cbind(results_limma, log2fc)

  # Get Uniprot_IDs
  Uniprot_IDs <- rownames(results_limma)
  head(Uniprot_IDs)

  # Map Uniprot IDs to Gene names
  symbol <- mapIds(org.Mm.eg.db, keys = Uniprot_IDs, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
  results_limma <- add_column(results_limma, symbol, .after = 1)

  # Need to write to results.
  return(results_limma)
}

#-------------------------------------------------------------------------------
#' ## normalize_TMM(data_in,colID)
# This function performs tmm normalization using the edgeR package on the columns
# specified by tmt_cols.

normalize_TMM <- function(data_in, groups) {
  for (i in 1:length(groups)) {
    colID <- groups[i]
    # Get data
    tmt_cols <- grep(colID, colnames(data_in))
    dm <- as.matrix(data_in[, tmt_cols])
    # replace NA with 0 (TMM cannot have missing values (NA))
    logic <- is.na(dm)
    if (length(logic[logic == TRUE]) > 0) {
      print("Warning: missing values (NA) are not tolerated. These will be replaced with 0")
    }
    dm[logic] <- 0
    data_work <- dm
    # TMM Normalization.
    factors_tmm <- calcNormFactors(data_work)
    dm_tmm <- sweep(data_work, 2, factors_tmm, FUN = "/")
    dm_tmm[logic] <- NA # Insure 0 values are now NA
    # Write to data frame
    data_in[, tmt_cols] <- dm_tmm
  }
  return(data_in)
}
#-------------------------------------------------------------------------------
# ## Define function: plotpdf(figure_name)
# This is a simple function for plotting figure as PDF format using CairoPDF.
# Create plot object with pryr %<a-%.
plotpdf <- function(figure_name, plot_object) {
  CairoPDF(file = paste(savefigs_path, "/", figure_name, sep = ""), width = 12, height = 8, paper = "special")
  plot_object
  dev.off()
}

#-------------------------------------------------------------------------------
## ggplotMeanSdPlot.

ggplotMeanSdPlot <- function(data_in, colID, title, log) {
  data_in <- data_in[, grepl(colID, colnames(data_in))]
  if (log == TRUE) {
    data_in <- as.matrix(log2(data_in))
  } else {
    data_in <- as.matrix(data_in)
  }
  ff <- tempfile()
  png(filename = ff)
  plot <- vsn::meanSdPlot(data_in)
  plot <- plot$gg
  dev.off()
  unlink(ff)
  plot <- plot + ggtitle(title) + theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold")
  )
  return(plot)
}

#-------------------------------------------------------------------------------
#' ## checkNormalization
# This function calls several of the custom ggplot functions for checking the data:
checkNormalization <- function(data_in, traits, colors, title) {
  # Generate plots.
  p1 <- ggplotBoxPlotv2(log2(data_in),
    colID = "b", traits = traits,
    colors = colors,
    title = title
  )
  p2 <- ggplotDensityv2(log2(data_in),
    colID = "b", traits = traits,
    colors = colors,
    title = title
  )
  p3 <- ggplotMDS(log2(data_in),
    colID = "b", traits,
    title = title
  ) + theme(legend.position = "none")
  p4 <- ggplotMeanSdPlot(log2(data_in), title = title)
  # Store in list.
  plot_list <- list(p1, p2, p3, p4)
  # Return plots.
  return(plot_list)
}

#-------------------------------------------------------------------------------
# Function for plotting sample connectivity.
ggplotSampleConnectivityv2 <- function(data_in, log = TRUE, colID, threshold = -2.5) {
  cols <- grep(colID, colnames(data_in))
  dm <- as.matrix(data_in[, cols])
  if (log == TRUE) {
    dm <- log2(dm)
  } else {
    dm <- dm
  }
  # Calcualte adjacency matrix.
  adjm <- (0.5 + 0.5 * bicor(dm, use = "pairwise.complete.obs")^2)

  # Generate a dendrogram.
  sampleTree <- hclust(as.dist(1 - adjm), method = "average")
  dendro <- ggdendrogram(sampleTree, theme_dendro = FALSE) + xlab("Sample") +
    ylab("Height") + ggtitle("Sample Clustering to detect outliers") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )

  # Calculate Ki (standardized connectivity)
  # The standardized connectivity (Z.K; Oldham et al.,) is a quantity that describes
  # the overall strength of connections between a given node (sample) and all of the other
  # nodes (samples) in a network.

  # The total connectivity of a Node (sample) is the sum of all of its connections (colSum).
  ki <- fundamentalNetworkConcepts(adjm)$Connectivity
  kmax <- max(ki)
  Ki <- ki / kmax # Normalized ki by maximum.
  Kmean <- mean(Ki)
  Kvar <- var(Ki)
  Z.Ki <- (Ki - Kmean) / sqrt(Kvar)

  # Gather data in df for plotting.
  data <- as.data.frame(Z.Ki)
  data$Sample <- sampleTree$order
  rownames(data) <- colnames(adjm)

  # Sort by standardized connectivity.
  data <- data[order(data$Z.Ki), ]
  data$label <- rownames(data)
  data$label[!data$Z.Ki < threshold] <- ""

  # Plot Z.K
  plot <- ggplot(data, aes(Sample, Z.Ki)) + geom_point(size = 2) + ggtitle("Sample Connectivity") +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red", size = 0.25) +
    geom_label_repel(aes(label = label), nudge_y = 0.3, colour = "red", alpha = 0.85) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(hjust = 0.5, color = "black", size = 11),
      axis.title.y = element_text(hjust = 0.5, color = "black", size = 11)
    )
  data$label <- NULL
  result <- list(data, dendro, plot)
  names(result) <- c("table", "dendrogram", "connectivityplot")
  return(result)
}

#-------------------------------------------------------------------------------
# Function for plotting sample connectivity.

outliers_Oldham <- function(data_in, colID, threshold = 2.5) {

  # get the data
  cols <- grep(colID, colnames(data_in))
  dm <- as.data.frame(data_in[, cols])
  rownames(dm) <- data_in$Accession
  dm <- log2(dm)

  # Calcualte adjacency matrix.
  adjm <- (0.5 + 0.5 * bicor(dm, use = "pairwise.complete.obs")^2)

  # Calculate Ki (standardized connectivity)
  # The standardized connectivity (Z.K; Oldham et al.,) is a quantity that describes
  # the overall strength of connections between a given node (sample) and all of the other
  # nodes (samples) in a network.

  # The total connectivity of a Node (sample) is the sum of all of its connections (colSum).
  ki <- fundamentalNetworkConcepts(adjm)$Connectivity
  kmax <- max(ki)
  Ki <- ki / kmax # Normalized ki by maximum.
  Kmean <- mean(Ki)
  Kvar <- var(Ki)
  Z.Ki <- (Ki - Kmean) / sqrt(Kvar)

  # Gather data in df for plotting.
  data <- as.data.frame(Z.Ki)
  data$Sample <- c(1:ncol(adjm))
  rownames(data) <- colnames(adjm)

  # Sort by standardized connectivity.
  data <- data[order(data$Z.Ki), ]

  idx <- data$Z.Ki < -1 * threshold
  out <- colnames(dm)[idx]
  print(out)
  return(data)
}

#-------------------------------------------------------------------------------
## Define function: store_plot
#  Saves ggplot(s) to list.
store_ggplot <- function(plot_list, plot, name) {
  k <- length(plot_list)
  plot_list[[k + 1]] <- plot
  names(plot_list)[k + 1] <- name
  return(plot_list)
}

#-------------------------------------------------------------------------------
#' ggplotMDS
#'
#' This function utilizes limma::plotMDS to generate a MDS plot which is then plotted with ggplot2.
#' The column names of the input data (an expression data frame) are used as geom_point() labels.
#' To supress plot output which results from calling limma::plotMDS, a temporary file
#' is created (see references).
#'
#' @param data_in the expression data frame.
#' @param colors colors for geom_point().
#' @param title a title for the plot.
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{https://stackoverflow.com/questions/20363266/how-can-i-suppress-the-creation-of-a-plot-while-calling-a-function-in-r}
#' @keywords fill down blank missing values
#'
#' @examples
#' ggplotMDS(data_in, colors, title)
#' @export
#' @importFrom limma
#'

ggplotMDS <- function(data_in,
                      colID,
                      colors,
                      title,
                      sample_info,
                      labels = FALSE) {
	require(limma, quietly = TRUE)
  # get the data
  cols <- grep(colID, colnames(data_in))
  dm <- as.matrix(data_in[, cols])
  idx <- match(colnames(dm), sample_info$ColumnName)
  simple_cols <- paste(sample_info$Model, sample_info$SampleType, sep = "_")
  # simple_cols <- sample_info$Genotype
  simple_cols <- gsub(" ", "", simple_cols)
  colnames(dm) <- simple_cols[idx]
  ff <- tempfile()
  png(filename = ff)
  data_MDS <- plotMDS(log2(dm))
  x <- data_MDS$x
  y <- data_MDS$y
  dev.off()
  unlink(ff)
  dm_MDS <- cbind(x, y)
  Condition <- rownames(dm_MDS)
  df_MDS <- as.data.frame(cbind(x, y))

  # Plot with no labels.
  plot <- ggplot(df_MDS, aes(x, y, color = Condition)) + geom_point(size = 3) +
    scale_color_manual(values = colors) +
    ggtitle(title) + xlab("Leading LogFC dim 1") + ylab("Leading LogFC dim 2") +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )

  # Plot with text labels.
  if (labels == TRUE) {
    plot <- NULL
    labs <- paste(sample_info$Model, sample_info$SampleType, sep = "_")[idx]
    plot <- ggplot(df_MDS, aes(x, y, color = labs)) + geom_text(aes(label = labs)) +
      scale_color_manual(values = colors) +
      ggtitle(title) + xlab("Leading LogFC dim 1") + ylab("Leading LogFC dim 2") +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
        axis.title.x = element_text(color = "black", size = 11, face = "bold"),
        axis.title.y = element_text(color = "black", size = 11, face = "bold")
      )
  }
  return(plot)
}

#-------------------------------------------------------------------------------
## Define function: ggplotMDS(data_in, colID, traits, title)

ggplotMDSv3 <- function(data_in, colID, traits, title) {
  idx <- match(colnames(data_in), rownames(traits))
  # colnames(data_in) <- traits$SampleType[idx]
  colnames(data_in) <- traits$Sample.Model[idx]
  ff <- tempfile()
  png(filename = ff)
  data_MDS <- plotMDS(data_in)
  x <- data_MDS$x
  y <- data_MDS$y
  dev.off()
  unlink(ff)
  dm_MDS <- cbind(x, y)
  Condition <- rownames(dm_MDS)
  df_MDS <- as.data.frame(cbind(x, y))
  plot <- ggplot(df_MDS, aes(x, y, color = Condition)) + geom_text(aes(label = Condition)) +
    ggtitle(title) + xlab("Leading LogFC dim 1") + ylab("Leading LogFC dim 2") +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}

#-------------------------------------------------------------------------------
# Annotate topTags table with entrez and gene symbols.
# based on rownames being Uniprot IDS.
# Function:
annotate_Entrez <- function(y_TT) {
  Uniprot <- rownames(y_TT)
  # Map Uniprot IDs to Entrez IDs:
  Entrez <- mapIds(org.Mm.eg.db, keys = Uniprot, column = "ENTREZID", keytype = "UNIPROT", multiVals = "first")
  Gene <- mapIds(org.Mm.eg.db, keys = Uniprot, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
  # Add columns for Entrez ID and gene symbols.
  y_TT <- add_column(y_TT, Entrez, .before = 1)
  y_TT <- add_column(y_TT, Uniprot, .before = 1)
  y_TT <- add_column(y_TT, Gene, .after = 2)
  return(y_TT)
}

#-------------------------------------------------------------------------------
## Define function: ggplotMDS(data_in, colID, traits, title)
# Does not log transform the data.

ggplotMDSv2 <- function(data_in, colID, traits, title) {
  idx <- match(colnames(data_in), rownames(traits))
  # colnames(data_in) <- traits$SampleType[idx]
  colnames(data_in) <- traits$Sample.Model[idx]
  # colnames(data_in) <- traits$PrepDate[idx]
  ff <- tempfile()
  png(filename = ff)
  data_MDS <- plotMDS(data_in)
  x <- data_MDS$x
  y <- data_MDS$y
  dev.off()
  unlink(ff)
  dm_MDS <- cbind(x, y)
  Condition <- rownames(dm_MDS)
  df_MDS <- as.data.frame(cbind(x, y))
  plot <- ggplot(df_MDS, aes(x, y, color = Condition)) + geom_text(aes(label = Condition), size = 6) +
    ggtitle(title) + xlab("Leading LogFC dim 1") + ylab("Leading LogFC dim 2") +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )

  # Dendrogram
  distMat <- dist(data_MDS$cmdscale.out)
  sampleTree <- flashClust(distMat, method = "complete")
  dendro <- ggdendrogram(sampleTree, rotate = TRUE, theme_dendro = FALSE) + xlab("Sample") +
    ylab("Height") + ggtitle("Sample Clustering based on MDS") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )

  data_out <- list(data_MDS, plot, dendro)
  names(data_out) <- c("data_MDS", "plot", "dendro")
  return(data_out)
}

#-------------------------------------------------------------------------------
## Define function: ggplotDensity(data_in,title)

ggplotDensity <- function(data_in, colID, title) {
  # The function ggplotDensity plots a density plot using ggplot.
  # The input data should be a data matrix. This data is reshaped with reshape2::melt
  # The title provided as input is used as the ggplot title.
  dm <- df2dm_TMT(data_in, colID)
  data_temp <- melt(log2(dm))
  colnames(data_temp) <- c("Accession", "Run", "Intensity")
  data_temp$Run <- as.factor(data_temp$Run)
  data_temp <- na.omit(data_temp)
  plot <- ggplot(data_temp, aes(x = Intensity, color = Run)) + geom_density() +
    ggtitle(title) +
    xlab("Log2 Intensity") +
    ylab("Density") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}

#-------------------------------------------------------------------------------
# Function for plotting distribution (density) of missing values).
ggplotDetect <- function(data_in, group, log = TRUE) {
  # Subset the data.
  cols <- grepl(group, colnames(data_in))
  data_work <- as.data.frame(data_in[, cols])
  # Remove rows will all missing values.
  data_work$NA_count <- apply(data_work, 1, function(x) sum(is.na(x)))
  data_work <- subset(data_work, !NA_count == (ncol(data_work) - 1))
  # Remove rows with missing QC data.
  # cols <- grep("QC",colnames(data_work))
  # data_work$NA_count <- apply(data_work[,cols],1,function(x) sum(is.na(x)))
  # data_work <- subset(data_work,!NA_count==length(cols))
  # Number of remaining missing values.
  data_work$NA_count <- apply(data_work[, c(1:ncol(data_work) - 1)], 1, function(x) sum(is.na(x)))
  # sum(data_work$NA_count)
  # Add Row avg.
  data_work$RowAverage <- apply(data_work[, c(1:ncol(data_work) - 1)], 1, function(x) mean(x, na.rm = TRUE))
  # Get the data for plotting.
  if (log == TRUE) {
    df <- as.data.frame(cbind(Avg = log2(data_work$RowAverage), Group = data_work$NA_count))
  } else {
    df <- as.data.frame(cbind(Avg = data_work$RowAverage, Group = data_work$NA_count))
  }
  # Groups are 1 (TRUE; No missing values), and 0 (FALSE; contains missing values)
  df$Group <- !df$Group == 0
  # df$Group[df$Group==1] <- "Complete"
  # df$Group[df$Group==0] <- "Missing"
  df$Group <- as.factor(df$Group)
  plot <- ggplot(df, aes(x = Avg, fill = Group, colour = Group)) +
    geom_density(alpha = 0.1, size = 1) + ggtitle(paste(group, "missing values")) +
    scale_fill_discrete(name = "Missing Values") +
    scale_color_discrete(name = "Missing Values") +
    xlab(expression(Log[2] ~ Intensity)) +
    ylab(expression(Density)) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}

#-------------------------------------------------------------------------------
## Define function: ggplotDensity(data_in,title)

ggplotDensityv2 <- function(data_in, colID, colors, traits, title) {
  sampleIndex <- as.data.frame(do.call(rbind, strsplit(colnames(cleanDat), "\\.")))
  colnames(sampleIndex) <- c("batch", "channel")
  batchIndex <- unique(sampleIndex$batch)
  model_colors <- as.data.frame(cbind(batchIndex,
    Model = unique(traits$Model),
    Colors = colors
  ))
  colors <- rep(model_colors$Colors, as.vector(rowSums(table(sampleIndex))))

  dm <- df2dm_TMT(data_in, colID)
  data_temp <- melt(dm)
  colnames(data_temp) <- c("Accession", "Run", "Intensity")
  data_temp$Run <- as.factor(data_temp$Run)
  data_temp <- na.omit(data_temp)
  plot <- ggplot(data_temp, aes(x = Intensity, color = Run)) + geom_density() +
    ggtitle(title) +
    xlab("Log2 Intensity") +
    ylab("Density") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    ) +
    scale_color_manual(values = colors)
  return(plot)
}

#-------------------------------------------------------------------------------
## Define function: ggplotBoxPlot(data_in,colors,title)

ggplotBoxPlot <- function(data_in, colID, colors, title) {
  # The function ggplotBoxPlot plots Run level boxplots with ggplot.
  dm <- df2dm_TMT(data_in, colID)
  data_temp <- melt(log2(dm))
  colnames(data_temp) <- c("Accession", "Run", "Intensity")
  data_temp$Run <- as.factor(data_temp$Run)
  data_temp <- na.omit(data_temp)

  # Discrete x-axis labels.
  #v <- seq(1,44,1)
  #v <- v[rep(c(FALSE,TRUE), 22)] <- ""

  plot <- ggplot(data_temp, aes(x = Run, y = Intensity, fill = Run)) +
    geom_boxplot(outlier.colour = "black", outlier.shape = 20, outlier.size = 1) +
    scale_fill_manual(
      values = colors,
      name = "Genotype",
      breaks = c(1, 12, 24, 35),
      labels = c("Syngap1", "Ube3a", "Shank2", "Shank3")
    ) +
    ggtitle(title) +
    xlab("TMT Run") +
    ylab("Log2 Intensity") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.text.x = element_text(color = "black", size = 8, angle = 45),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}



#-------------------------------------------------------------------------------
## boxplot function.

ggplotBoxPlotv2 <- function(data_in, colID, traits, colors, title) {
  sampleIndex <- as.data.frame(do.call(rbind, strsplit(colnames(data_in), "\\.")))
  colnames(sampleIndex) <- c("batch", "channel")
  batchIndex <- unique(sampleIndex$batch)
  model_colors <- as.data.frame(cbind(batchIndex,
    Model = unique(traits$Model),
    Colors = colors
  ))
  colors <- rep(model_colors$Colors, as.vector(rowSums(table(sampleIndex))))

  dm <- df2dm_TMT(data_in, colID)
  data_temp <- melt(dm)
  colnames(data_temp) <- c("Accession", "Run", "Intensity")
  data_temp$Run <- as.factor(data_temp$Run)
  data_temp <- na.omit(data_temp)

  plot <- ggplot(data_temp, aes(x = Run, y = Intensity, fill = Run)) +
    geom_boxplot(outlier.colour = "black", outlier.shape = 20, outlier.size = 1) +
    scale_fill_manual(
      values = colors,
      name = "Genotype",
      breaks = c(cumsum(as.vector(rowSums(table(sampleIndex))))),
      labels = paste(model_colors$Model, " (", model_colors$batchIndex, ")", sep = "")
    ) +
    ggtitle(title) +
    xlab("TMT Run") +
    ylab("Log2 Intensity") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}


#-------------------------------------------------------------------------------
## Define function: ggsave_plots(plot_list,file_ext)

ggsaveMultiPDF <- function(plot_list, savefigs_path, file_ext, verbose = TRUE) {
  # Saves plots in plot_list to file as given by savefigs_path.
  # Saves figures in format given by file_ext.
  # Will name plots as plot1, plot2, ect if not named in list.
  # If name exists in directory, then plot name will be plotn+1
  # Coerce to list.
  if (!inherits(plot_list, "list")) {
    plot_list <- list(plot_list)
  }
  # Loop to save plots.
  for (i in 1:length(plot_list)) {
    # Check if name exists.
    if (is.null(names(plot_list))) {
      save_name <- paste("plot", i, file_ext, sep = "")
      # If name exists, rename as plot(N+1)
      if (save_name %in% list.files(savefigs_path)) {
        idx <- grepl("plot", list.files(savefigs_path))
        num <- max(as.numeric(sapply(strsplit(list.files(savefigs_path)[idx], "plot|\\."), "[", 2))) + 1
        save_name <- paste("plot", num, file_ext, sep = "")
      }
    } else {
      save_name <- paste(names(plot_list)[i], file_ext, sep = "")
    }
    # Save the plots.
    suppressMessages(ggsave(save_name,
      plot = plot_list[[i]],
      scale = 1,
      width = NA,
      height = NA,
      path = savefigs_path
    ))
  }
  # Progress message:
  if (verbose == TRUE) {
    print(paste(length(plot_list), "plot(s) saved as", file_ext, "in", savefigs_path))
  }
}

#-------------------------------------------------------------------------------
# Function for saving ggplots to single pdf.
ggsavePDF <- function(plots, file) {
  pdf(file, onefile = TRUE)
  # If not a list, coerce to list.
  if (!inherits(plots, "list")) {
    plot_list <- list(plots)
  } else {
    plot_list <- plots
  }
  # Loop through list, save plots to pdf.
  for (i in 1:length(plot_list)) {
    # If ggplot, then print.
    plot <- plot_list[[i]]
    if (inherits(plot, "ggplot")) {
      print(plot)
      # If else, coerce to ggplot and print.
    } else {
      print(as.ggplot(plot))
    }
  }
  quiet(dev.off())
}

#-------------------------------------------------------------------------------
#' ## ggplotQCHist(data_in,groups,nbins,threshold)

ggplotQCHist <- function(data_in, group, nbins, threshold) {
  hist_list <- list()
  # Get the data for the specified group.
  cols <- grep(group, colnames(data_in))
  data_sub <- data_in[, cols]
  tmt_cols <- grep("Abundance", colnames(data_sub))
  qc_cols <- grep("QC", colnames(data_sub))

  # Remove row if QC was not quantified in all QC replicates.
  out <- is.na(data_sub[, qc_cols])
  logic_out <- matrix(NA, nrow = nrow(data_sub), ncol = 1)
  for (j in 1:nrow(data_sub)) {
    logic_out[j] <- any(out[j, ])
  }

  # Write to data_sub
  data_sub[logic_out, ] <- NA

  # Log2 transform and define tmt and qc columns
  data_work <- log2(data_sub)
  cols_qc <- grep("QC", colnames(data_work))
  cols_tmt <- grep("Abundance", colnames(data_work))

  # Calculate the ratio of QC samples (QC1-QC2, QC1-QC3, QC2-QC3)
  data_work$ratio_1 <- data_work[, 1] - data_work[, 2]
  data_work$ratio_2 <- data_work[, 1] - data_work[, 3]
  data_work$ratio_3 <- data_work[, 2] - data_work[, 3]

  cols <- grep("ratio_", colnames(data_work))
  cols_ratio <- (last(cols_tmt) + 1):(last(cols_tmt) + 3)
  data_work$ratio_avg <- rowMeans(data_work[, cols_ratio], 1, na.rm = TRUE)
  data_work$avgQC <- rowMeans(data_work[, cols_qc], 1, na.rm = TRUE)
  data_work$sdQC <- apply(data_work[, cols_ratio], 1, FUN = sd, na.rm = TRUE)

  # Calculate bins based on mean intensity of QC replicates.
  rows_ignore <- is.nan(data_work$avgQC)
  data_work$bins[!rows_ignore] <- BurStMisc::ntile(data_work$avgQC[!rows_ignore], nbins,
    na.rm = TRUE,
    checkBleed = FALSE, result = "numeric"
  )

  # Calculate summary statistics of Intensity bins
  data_temp <- data_work[!rows_ignore, ]
  data_QC <- subset(data_temp) %>%
    group_by(bins) %>%
    dplyr::summarize(
      min = min(ratio_avg, na.rm = TRUE),
      max = max(ratio_avg, na.rm = TRUE),
      avg = mean(ratio_avg, na.rm = TRUE),
      stdev = sd(ratio_avg, na.rm = TRUE)
    )
  # Add bin avg and std to data_temp.
  data_temp$bin_avg <- NA
  data_temp$bin_std <- NA
  for (i in 1:nbins) {
    data_temp$bin_avg[data_temp$bins == i] <- data_QC$avg[data_QC$bins == i]
    data_temp$bin_std[data_temp$bins == i] <- data_QC$stdev[data_QC$bins == i]
  }
  # Add upper and lower bounds.
  data_temp$upper <- 4 * data_temp$bin_std + data_temp$bin_avg
  data_temp$lower <- -4 * data_temp$bin_std + data_temp$bin_avg
  # Out if ratio is outside lower/upper bounds.
  data_temp$out1 <- data_temp$ratio_avg > data_temp$upper
  data_temp$out2 <- data_temp$ratio_avg < data_temp$lower
  # Check if value is outside.
  data_temp$out <- as.numeric(apply(data_temp, 1, function(x) any(x[c(23, 24)])))
  # Count number of peptides to be removed.
  data_QC2 <- subset(data_temp) %>%
    group_by(bins) %>%
    dplyr::summarize(
      min = min(ratio_avg, na.rm = TRUE),
      max = max(ratio_avg, na.rm = TRUE),
      avg = mean(ratio_avg, na.rm = TRUE),
      stdev = sd(ratio_avg, na.rm = TRUE),
      nout = sum(out, na.rm = TRUE)
    )
  # Loop through bins and generate plots.
  palette <- c("#132B43", "#22496C", "#336A98", "#448DC6", "#56B1F7")
  for (j in 1:nbins) {
    data_sub <- subset(data_work, bins == j)
    dm <- df2dm_TMT(data_sub, "ratio_avg")
    data_temp <- na.omit(melt(dm))
    mu <- mean(data_temp$value)
    stdev <- sd(data_temp$value)
    plot <- ggplot(data = data_temp, aes(value)) + geom_histogram(bins = 100, fill = palette[j]) +
      ggtitle(paste("Intensity bin =", j)) +
      geom_vline(xintercept = mu + 4 * stdev, linetype = "dashed", color = "red", size = 0.75) +
      geom_vline(xintercept = mu - 4 * stdev, linetype = "dashed", color = "red", size = 0.75) +
      theme(
        plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
        axis.title.x = element_text(color = "black", size = 11, face = "bold"),
        axis.title.y = element_text(color = "black", size = 11, face = "bold")
      )
    plot <- plot + xlim(c(-1.5, 1.5))
    # Add annotation layer.
    xrange <- unlist(ggplot_build(plot)$layout$panel_params[[1]][1])
    yrange <- unlist(ggplot_build(plot)$layout$panel_params[[1]][8])
    xmin <- min(xrange)
    xmax <- max(xrange)
    xdelta <- xmax - xmin
    ymin <- min(yrange)
    ymax <- max(yrange)
    ydelta <- ymax - ymin
    tt <- ttheme_default(base_size = 11, core = list(bg_params = list(fill = "white")))
    mytable <- t(subset(data_QC2, bins == j))
    mytable <- round(mytable, 2)
    rownames(mytable) <- c("Bin", "Min", "Max", "Mean", "SD", "Nout")
    mytable <- paste(rownames(mytable), mytable, sep = " = ")
    tab <- tableGrob(mytable, rows = NULL, theme = tt)
    g <- gtable_add_grob(tab,
      grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
      t = 1, b = nrow(tab), l = 1, r = ncol(tab)
    )
    plot <- plot + annotation_custom(g,
      xmin = xmin + 0.80 * xdelta, xmax,
      ymin = ymin + 0.6 * ydelta, ymax
    )
    hist_list[[j]] <- plot
  }

  return(hist_list)
}

#-------------------------------------------------------------------------------
## Define function: ggplotHist
# Plot a histogram using ggplot.

ggplotHist <- function(data_in, colID, title) {
  dm <- df2dm_TMT(data_in, colID)
  data_temp <- na.omit(melt(log2(dm)))
  colnames(data_temp) <- c("Accession", "Run", "Intensity")
  plot <- ggplot(data = data_temp, aes(Intensity)) + geom_histogram(bins = 100, fill = "gray") +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
}

#-------------------------------------------------------------------------------
#' Define function: ggplotPvalHist()
ggplotPvalHist <- function(data_in, color, title) {
  col <- grep("PValue|pvalue|p value| P value| P Value", colnames(data_in))
  plot <- ggplot(data = data_in, aes(data_in[, col])) +
    geom_histogram(bins = 100, fill = color, col = I("black")) +
    ggtitle(title) + xlab("PValue") + ylab("Count") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}

#-------------------------------------------------------------------------------
#' Define function: ggplotQQplot()
ggplotQQplot <- function(data_in, title) {
  df <- na.omit(melt(as.data.frame(data_in)))
  plot <- ggplot(df, aes(sample = value)) + stat_qq() + stat_qq_line() +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}

#-------------------------------------------------------------------------------
## ggplotPCA

ggplotPCA <- function(data_in, traits, colors, title = "2D PCA Plot") {
  # Perform PCA
  PC <- prcomp(t(na.omit(data_in)))$x[, 1:2]
  # Add annotations to PC data frame.
  PC <- as.data.frame(PC)
  PC$Label <- paste(traits$Model, traits$SampleType, sep = "_")[match(rownames(PC), traits$ColumnName)]

  # Generate plot.
  plot <- ggplot(PC, aes(x = PC1, y = PC2)) +
    geom_text(aes(label = Label), color = colors) +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}

#-------------------------------------------------------------------------------
# Define function: ggplotPCA
# Generates a PCA plot.

ggplotPCAv2 <- function(data_in, traits, title = "2D PCA Plot") {
  # Perform PCA
  PC <- prcomp(t(na.omit(data_in)))$x[, 1:2]
  # Add annotations to PC data frame.
  PC <- as.data.frame(PC)
  PC$Label <- traits$Sample.Model[match(rownames(PC), rownames(traits))]
  colors <- traits$Color[match(rownames(PC), rownames(traits))]
  # Generate plot.
  plot <- ggplot(PC, aes(x = PC1, y = PC2)) +
    geom_text(aes(label = Label), color = colors) +
    ggtitle(title) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )

  # Generate dendrogram
  distMat <- dist(PC)
  sampleTree <- flashClust(t(distMat), method = "complete")
  dendro <- ggdendrogram(sampleTree, rotate = TRUE, theme_dendro = FALSE) + xlab("Sample") +
    ylab("Height") + ggtitle("Sample Clustering based on PCA") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )

  # Generate output
  data_out <- list(PC, plot, dendro)
  names(data_out) <- c("data_PCA", "plot", "dendro")

  return(data_out)
}

#-------------------------------------------------------------------------------
#' Define function: ggplotVolcanoPlot()

ggplotVolcanoPlot <- function(data_in, title, cutoff = log2(1.25)) {
  df <- data_in
  df$x <- df[, grep("FC", colnames(df))]
  df$y <- -log10(df[, grep("PValue", colnames(df))])
  logic <- (df$y > 1.30103) & (df$x > cutoff || df$x < cutoff)
  df$color <- "ns"
  df$color[logic] <- "sig"
  plot <- ggplot(data = df, aes(x = x, y = y, color = "blue")) +
    geom_point(aes(color = df$color)) +
    geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black", size = 0.6) +
    geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", size = 0.6) +
    geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 0.6) +
    ggtitle(title) + xlab("Log2FoldChange") + ylab("-Log10PValue") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold"),
      legend.position = "none"
    )
  return(plot)
}

#-------------------------------------------------------------------------------
# Function to make boxplots for proteins of interst.
ggplotProteinBoxes <- function(data_in, interesting.proteins, dataType, traits, order = NULL, scatter = FALSE) {
  proteinBoxes <- list()
  for (protein in interesting.proteins) {
    # Subset the data.
    idx <- match(protein, rownames(data_in))
    data_sub <- as.data.frame(data_in[idx, ])
    colnames(data_sub) <- "Intensity"
    data_sub$Group <- traits$Sample.Model[match(rownames(data_sub), rownames(traits))]
    if (is.numeric(order)) {
      levels <- unique(data_sub$Group)[order]
      data_sub$Group <- factor(data_sub$Group, levels = levels)
    } else {
      data_sub$Group <- factor(data_sub$Group)
    }
    plot <- ggplot(data_sub, aes(x = Group, y = Intensity, fill = Group)) +
      geom_boxplot(outlier.colour = "black", outlier.shape = 20, outlier.size = 1) +
      ggtitle(protein) + ylab(paste0("Log2", "(", dataType, ")")) + xlab("") +
      theme(
        plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
        axis.title.x = element_text(color = "black", size = 10),
        axis.title.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    # Add scatter
    if (scatter == TRUE) {
      plot <- plot + geom_point(aes(x = Group, y = Intensity), color = "white", size = 2, pch = 21, fill = "black")
    }
    # Store in list.
    proteinBoxes[[protein]] <- plot
  }
  return(proteinBoxes)
}

#-------------------------------------------------------------------------------
## A function to Write to excel.
write.excel <- function(data, file, rowNames = FALSE, colNames = TRUE) {
  if (class(data) == "list") {
    list <- data
  } else {
    # Coerce to list.
    list <- list(data)
  }
  # Insure there are names.
  if (is.null(names(list))) {
    names(list) <- paste("Sheet", c(1:length(list)))
  }
  wb <- createWorkbook()
  # Loop to add a worksheets:
  for (i in 1:length(list)) {
    df <- as.data.frame(list[[i]])
    addWorksheet(wb, sheetName = names(list[i]))
    writeData(wb, sheet = i, df, rowNames = rowNames, colNames = colNames)
  }
  # Save workbook.
  saveWorkbook(wb, file, overwrite = TRUE)
}

#-------------------------------------------------------------------------------
# Function for doing eBLM regression and exactTest.

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


#-------------------------------------------------------------------------------
# Filter proteins.
filter_proteins <- function(data_in, colID) {
  # Removing one hit wonders...
  out <- data_in$Peptides == 1
  print(paste(
    length(out[out == TRUE]),
    "proteins are identified by only one peptide and will be removed."
  ))
  filt_protein <- data_in[!out, ]

  # Removing proteins that are not identified in at least 50% of samples.
  tmt_cols <- grep(colID, colnames(filt_protein))
  threshold <- length(colnames(filt_protein)[tmt_cols]) / 2
  out <- apply(filt_protein[, tmt_cols], 1, function(x) sum(is.na(x))) > threshold
  # Number of proteins identified in less than 50% of samples
  print(paste(
    length(out[out == TRUE]),
    "proteins are identified in less than 50% of samples and are removed."
  ))
  filt_protein <- filt_protein[!out, ]
  return(filt_protein)
}

#-------------------------------------------------------------------------------
#' ## impute_KNN(data_in)
#' A function for imputing TMT protein expression matrix with the KNN algorithm
#' from package impute. Note data is log transformed before imputing and then
#' return un-logged.

# KNN impute.
impute_KNN <- function(data_in, colID) {
  cols <- grep(colID, colnames(data_in))
  data_work <- log2(as.matrix(data_in[, cols]))
  rownames(data_work) <- data_in$Accession
  data_imp <- 2.^impute.knn(data_work)$data
  data_out <- data_in
  data_out[, cols] <- data_imp
  return(data_out)
}

#-------------------------------------------------------------------------------
# A function for performing moderated empirical bayes regression for removal of covariate effects.
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

#-------------------------------------------------------------------------------
# A Function for performing GO and KEGG GSE.
# Requires an internet connection!
# Can operate on qlf or ET objects.

edgeR_GSE <- function(qlf, FDR = 0.05, filter = TRUE) {
  GSE_results <- list()
  # Map Uniprot IDs to Entrez IDs:
  Entrez <- mapIds(org.Mm.eg.db,
    keys = rownames(qlf),
    column = "ENTREZID",
    keytype = "UNIPROT",
    multiVals = "first"
  )
  Entrez[is.na(Entrez)] <- paste("not_mapped", c(1:sum(is.na(Entrez))), sep = "_")
  # Perform GO enrichment testing.
  GO <- goana(qlf, geneid = Entrez, species = "Mm", FDR = FDR)
  # Calculate adjusted p.values.
  GO$P.Adj.Up <- p.adjust(GO$P.Up, method = "hochberg")
  GO$P.Adj.Down <- p.adjust(GO$P.Down, method = "hochberg")
  # Perform KEGG enrichment testing.
  KEGG <- kegga(qlf, geneid = Entrez, species = "Mm", FDR = FDR)
  KEGG$P.Adj.Up <- p.adjust(KEGG$P.Up, method = "hochberg")
  KEGG$P.Adj.Down <- p.adjust(KEGG$P.Down, method = "hochberg")
  # Eliminate insignificant results if filter=TRUE.
  if (filter == TRUE) {
    GO <- GO[GO$P.Adj.Up < 0.05 | GO$P.Adj.Down < 0.05, ]
    KEGG <- KEGG[KEGG$P.Adj.Up < 0.05 | KEGG$P.Adj.Down < 0.05, ]
  }
  GSE_results$GO <- GO
  GSE_results$KEGG <- KEGG
  return(GSE_results)
}

#-------------------------------------------------------------------------------
# Function to annotate DE candidates:
annotateTopTags <- function(y_TT) {
  y_TT$logCPM <- 100 * (2^y_TT$logFC)
  colnames(y_TT)[2] <- "%WT"
  y_TT$candidate <- "no"
  y_TT[which(y_TT$FDR <= 0.10 & y_TT$FDR > 0.05), dim(y_TT)[2]] <- "low"
  y_TT[which(y_TT$FDR <= 0.05 & y_TT$FDR > 0.01), dim(y_TT)[2]] <- "med"
  y_TT[which(y_TT$FDR <= 0.01), dim(y_TT)[2]] <- "high"
  y_TT$candidate <- factor(y_TT$candidate, levels = c("high", "med", "low", "no"))
  return(y_TT)
}

#-------------------------------------------------------------------------------
# A function to perform QLF test with GLM fit.
results_QLFTest <- function(fit, comparison, alpha = 0.05) {
  # test comparisons with QLFtest, add FDR, extract results, and sort by pvalue.
  qlf <- glmQLFTest(fit, contrast = comparison)
  res <- qlf
  geno <- strsplit(strsplit(res$comparison, "\\.")[[1]][2], "\\ ")[[1]][1]
  res <- topTags(res, n = Inf, sort.by = "none")$table
  res <- res[order(res$PValue), ]
  # Determine number of significant results.
  summary_table <- summary(decideTests(qlf))
  # Convert logCPM to percent WT.
  res$logCPM <- round(100 * (2^res$logFC), 2)
  colnames(res)[2] <- "%WT"
  colnames(res)[3] <- "F Value"
  # Categorize candidates by FDR.
  res$candidate <- "no"
  res[which(res$FDR <= 0.10 & res$FDR > 0.05), dim(res)[2]] <- "low"
  res[which(res$FDR <= 0.05 & res$FDR > 0.01), dim(res)[2]] <- "med"
  res[which(res$FDR <= 0.01), dim(res)[2]] <- "high"
  res$candidate <- factor(res$candidate, levels = c("high", "med", "low", "no"))
  # Map Uniprot IDs to Gene names
  Uniprot_IDs <- sapply(strsplit(rownames(res), "\\|"), "[", 2)
  if (any(is.na(Uniprot_IDs))) {
    Uniprot_IDs <- rownames(res)
  }
  symbol <- mapIds(org.Mm.eg.db, keys = Uniprot_IDs, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
  res <- add_column(res, symbol, .before = 1)
  colnames(res)[1] <- "Symbol"

  # Extract comparison name from comparison (contrast) matrix
  name <- gsub(" ", "", strsplit(colnames(comparison), "\\.|-")[[1]][2])
  # Subset significant results.
  res_sub <- subset(res, res[, 5] < alpha)
  colnames(res)[c(2:6)] <- paste(geno, colnames(res)[c(2:6)])
  results_list <- list(summary_table, res, qlf, res_sub)
  names(results_list) <- c("summary_table", "QLFTest_results", "QLF", "Sig_results")
  return(results_list)
}


#-------------------------------------------------------------------------------
# A Function for performing eBLM regression and EdgeR GLM and QLFTest.
# Function for doing all the work.
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

#-------------------------------------------------------------------------------
#' ## ggplotBCV

ggplotBCV <- function(y_DGE) {
  y <- y_DGE
  tag_x <- y$AveLogCPM
  tag_y <- sqrt(y$tagwise.dispersion)
  common <- sqrt(y$common.dispersion)
  trend_y <- sqrt(y$trended.dispersion)
  trend_x <- y$AveLogCPM
  o <- order(trend_x)
  fill <- rep(1, nrow = length(trend_x), ncol = 1)
  df <- as.data.frame(cbind(tag_x, tag_y, trend_x, trend_y), fill)

  plot <- ggplot(df, aes(tag_x, tag_y, colour = "Tagwise")) + geom_point(size = 1) +
    geom_line(mapping = aes(trend_x[o], trend_y[o], colour = "Trend"), size = 0.75) +
    geom_hline(aes(yintercept = common, colour = "Common"), linetype = "solid", size = 0.75) +
    scale_colour_manual(name = "Dispersion", values = c("blue", "black", "red")) +
    xlab("Average log CPM") + ylab("Biological coefficient of variation") +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}


#-------------------------------------------------------------------------------
# Get GO and KEGG results from qlf object.
qlf_goana_results <- function(qlf_result) {
  qlf <- qlf_result$QLF
  Uniprot_IDs <- sapply(strsplit(rownames(qlf), "\\|"), "[", 2)
  Entrez <- mapIds(org.Mm.eg.db, keys = Uniprot_IDs, column = "ENTREZID", keytype = "UNIPROT", multiVals = "first")
  Entrez[is.na(Entrez)] <- paste("not_mapped", c(1:sum(is.na(Entrez))), sep = "_")
  rownames(qlf) <- Entrez
  go <- goana(qlf, species = "Mm")
  go$FDR.Up <- p.adjust(go$P.Up, method = "BH")
  go$FDR.Down <- p.adjust(go$P.Down, method = "BH")
  keg <- kegga(qlf, species = "Mm")
  keg$FDR.Up <- p.adjust(keg$P.Up, method = "BH")
  keg$FDR.Down <- p.adjust(keg$P.Down, method = "BH")
  results <- list(go, keg)
  names(results) <- c("GO", "KEGG")
  return(results)
}

#-------------------------------------------------------------------------------
# Function for slicing data into groups.
# Requires plyr
slice <- function(input, by = 2) {
  starts <- seq(1, length(input), by)
  tt <- lapply(starts, function(y) input[y:(y + (by - 1))])
  llply(tt, function(x) x[!is.na(x)])
}

#-------------------------------------------------------------------------------
# Function for supressing printed messages from a function.
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#-------------------------------------------------------------------------------
# Function to impute peptide level missing values.
# Supports MLE for MAR and KNN for MNAR.

impute_Peptides <- function(data_in, groups, method, qc_threshold = 0, bio_threshold = 2) {
  n_out <- list()
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
    n_out[[i]] <- num_NA
    # KNN Impute
    if (method == "knn") {
      data_imp <- as.matrix(data_temp)
      data_imp[, c(1:11)] <- quiet(impute.knn(data_imp[, c(1:11)])$data)
      data_imp <- as.data.frame(data_imp)
    } else if (method == "mle") {
      # MLE Impute
      conditions <- as.factor(c(rep(1, 3), rep(2, 4), rep(3, 4)))
      data_imp <- data_temp
      data_imp[, c(1:11)] <- impute.mle(data_temp[, c(1:11)], conditions)
    }
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
  names(n_out) <- groups
  results <- list(n_out, data_in)
  names(results) <- c("n_out", "data_imputed")
  return(results)
}

#-------------------------------------------------------------------------------
# Function ggplotVerboseBoxplot
ggplotVerboseBoxplot <- function(x, g, levels, contrasts, color, stats = FALSE,
                                 method = "dunnett", correction_factor = 8) {

  # Bind data as data frame for ggplot.
  df <- as.data.frame(cbind(x, g))
  df$g <- factor(df$g, levels = levels)

  # Calculate Kruskal Wallis pvalue.
  KWtest <- kruskal.test(as.numeric(x), as.factor(g))
  pvalue <- round(KWtest$p.value, 3)

  # If p-value is significant, print the title in red.
  if (as.numeric(KWtest$p.value) < 0.05) {
    sigcolor <- "red"
  } else {
    sigcolor <- "black"
  }

  # Post-hoc comparisons module. Dunn or Dunnett.
  if (method == "dunn") {
    # print("Dunn's test used.")
    # Dunn's post-hoc test.
    Dtest <- dunnTest(as.numeric(x) ~ as.factor(g), kw = FALSE, method = "none")$res
    # Duplicate comparisons in reverse order. Keep rows that have comparison of interest.
    dupDtest <- Dtest
    dupDtest$Comparison <- do.call(
      rbind,
      lapply(
        strsplit(Dtest$Comparison, " - "),
        function(x) paste(c(x[2], x[1]), collapse = " - ")
      )
    )
    Dtest <- rbind(Dtest, dupDtest)
    # Keep only contrasts of interest.
    Dtest <- Dtest[Dtest$Comparison %in% contrasts, ]
  } else if (method == "dunnett") {
    # print("Dunnett's test used.")
    # Dunnett's post-hoc test for comparison to controls.
    Dtest <- DunnettTest(as.numeric(x), as.factor(g), control = c("WT.Cortex", "WT.Striatum"))
    Dtest_list <- list()
    # Extract results from PostHocTest object.
    for (i in 1:length(Dtest)) {
      Dtest_list[[names(Dtest)[i]]] <- Dtest[[i]]
    }
    Dtest <- as.data.frame(do.call(rbind, Dtest_list))
    # Reverse order of contrasts.
    Dtest <- add_column(Dtest,
      Comparison = paste(sapply(strsplit(rownames(Dtest), "-"), "[", 2),
        sapply(strsplit(rownames(Dtest), "-"), "[", 1),
        sep = " - "
      ), .before = 1
    )
    # Discard row names.
    rownames(Dtest) <- NULL

    # keep only contrasts of interest.
    Dtest <- Dtest[Dtest$Comparison %in% contrasts, ]
    # Add p.adj
    Dtest$P.adj <- p.adjust(Dtest$pval, method = "BH")
  } else {
    print("Please specify a post-hoc test (posthoc = dunn or dunnett).")
  }

  # Prepare annotation df.
  stats.df <- Dtest
  stats.df$P.adj <- as.numeric(stats.df$P.unadj) * correction_factor # BH adjustment for 8 comparisons.
  stats.df$symbol <- ""
  stats.df$symbol[stats.df$P.adj < 0.05] <- "*"
  stats.df$symbol[stats.df$P.adj < 0.01] <- "**"
  stats.df$symbol[stats.df$P.adj < 0.001] <- "***"
  stats.df$group1 <- gsub(" ", "", sapply(strsplit(stats.df$Comparison, "-"), "[", 1))
  stats.df$group2 <- gsub(" ", "", sapply(strsplit(stats.df$Comparison, "-"), "[", 2))
  stats.df$ypos <- 1.05 * max(as.numeric(df$x))

  # Generate boxplot.
  plot <- ggplot(df, aes(x = g, y = as.numeric(x), fill = g)) +
    geom_boxplot() +
    scale_fill_manual(values = rep(color, length(unique(g)))) + geom_point(color = "white", size = 1, pch = 21, fill = "black") +
    ggtitle(paste(color, " (P = ", pvalue, ")", sep = "")) + xlab(NULL) +
    ylab("Summary Expression") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  # Extract ggplot build information to adjust y axis.
  build <- ggplot_build(plot)
  y_lims <- build$layout$panel_scales_y[[1]]$range$range
  y_lims[2] <- 1.15 * y_lims[2]
  plot <- plot + ylim(y_lims)

  # Add asterisks indicating significance.
  plot <- plot + annotate("text",
    x = stats.df$group2, y = stats.df$ypos,
    label = stats.df$symbol, size = 6
  )

  results <- list(plot, KWtest, stats.df)
  names(results) <- c("plot", "Kruskal-Wallis", method)
  if (stats == TRUE) {
    return(results)
  } else {
    return(plot)
  }
}

#-------------------------------------------------------------------------------
# Function for plotting overlap frequency for peptides and proteins:
ggplotFreqOverlap <- function(data_in, colID, groups) {
  # Subset the data.
  cols <- grepl(colID, colnames(data_in))
  data_work <- data_in[, cols]
  rownames(data_work) <- paste(data_in$Accession,
    c(1:nrow(data_in)),
    sep = "_"
  )
  # Logical matrix. 1 if expressed. 0 if NA, missing.
  data_logic <- !is.na(data_work)

  # Determine if protein was expressed in an experiment (all=TRUE)
  cols_list <- lapply(as.list(groups), function(x) grepl(x, colnames(data_logic)))
  logic <- lapply(cols_list, function(x) apply(data_logic[, x], 1, function(y) all(y)))
  all_expressed <- as.data.frame(do.call(cbind, logic))
  colnames(all_expressed) <- groups

  # Sum rows. The number of experiments a protein/peptide was expressed in.
  all_expressed$Frequency <- rowSums(all_expressed)

  # Summarize with table, and convert to df.
  df <- as.data.frame(table(all_expressed$Frequency))
  # Ignore row with 0, not expressed in any experiment.
  df <- df[!df$Var1 == 0, ]
  # Calculate percentage of total.
  df$Percent <- paste(round(100 * (df$Freq / sum(df$Freq)), 1), "%", sep = "")
  df$ypos <- 0.1 * max(df$Freq)

  plot <- ggplot(df, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_col() +
    scale_fill_grey(start = 0.8, end = 0.2) +
    labs(
      title = "Identification Overlap",
      x = "Experiment Overlap",
      y = "Frequency"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold"),
      legend.position = "none"
    )

  # Add percent of total annotation.
  plot <- plot + annotate("text",
    x = df$Var1, y = df$ypos,
    label = df$Percent, size = 3, color = "black"
  )
  return(plot)
}

#-------------------------------------------------------------------------------
ggplotSigOverlap <- function(data, protID, colID) {
  cols <- grepl(colID, colnames(data))
  data_work <- data[, cols]
  idcol <- grepl(protID, colnames(data))
  rownames(data_work) <- data[, idcol]
  logic <- as.data.frame(data_work < 0.05)
  rownames(logic) <- data[, idcol]
  logic$Freq <- rowSums(logic)
  df <- as.data.frame(table(logic$Freq))
  df <- df[!df$Var1 == 0, ]
  df$ypos <- 0.1 * max(df$Freq)
  df$label <- paste("N=", df$Freq)
  head(df)
  plot <- ggplot(df, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_col() +
    scale_fill_grey(start = 0.8, end = 0.2) +
    labs(
      title = "Differentially Expressed Protein Overlap",
      x = "Experiment Overlap",
      y = "Protein Frequency"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold"),
      legend.position = "none"
    )
  # Add percent of total annotation.
  plot <- plot + annotate("text",
    x = df$Var1, y = df$ypos,
    label = df$label, size = 4, color = "black"
  )
  return(plot)
}

#-------------------------------------------------------------------------------

# Define function: ggplotVerboseScatterPlot
ggplotVerboseScatterPlot <- function(MMdata, GSdata, moduleGenes, module, trait, stats = TRUE) {

  # Subset the data.
  col <- grep("color|Color", colnames(moduleGenes))
  moduleGenes <- subset(moduleGenes, moduleGenes[, col] == module)

  # Add trait annotation.
  moduleGenes$trait <- trait

  # Y data = GeneSignificance.
  idx <- match(moduleGenes$geneNames, rownames(GSdata))
  ydat <- GSdata[idx, ]
  col <- grepl(trait, colnames(ydat))
  moduleGenes$GS <- ydat[, col]

  # X data = ModuleMembership.
  col <- match(module, colnames(MMdata))
  MM_sub <- as.matrix(MMdata[, col])
  colnames(MM_sub) <- module
  rownames(MM_sub) <- rownames(MMdata)
  idx <- match(moduleGenes$geneNames, rownames(MM_sub))
  moduleGenes$MM <- MM_sub[idx, ]

  # Determine best fit line.
  fit <- lm(moduleGenes$GS ~ moduleGenes$MM)
  # Determine intercepts and slope.
  # coef(fit)
  # Calculate Pearson P-Value.
  corTest <- cor.test(~ moduleGenes$GS + moduleGenes$MM,
    data = cbind(moduleGenes$MM, moduleGenes$GS),
    method = "pearson", conf.level = 0.95
  )
  # Gather line stats for an annotation layer.
  lmstats <- cbind(corTest$p.value, cor(moduleGenes$GS, moduleGenes$MM), coef(fit)[1], coef(fit)[2])
  rownames(lmstats) <- paste(module, trait, sep = "|")
  colnames(lmstats) <- c("p.value", "R2", "intercept", "slope")
  pvalue <- paste("P-value =", formatC(corTest$p.value, format = "e", digits = 2))
  slope <- paste("Slope =", round(as.numeric(coef(fit)[2]), 3))
  R2 <- paste("R2 =", round(as.numeric(cor(moduleGenes$GS, moduleGenes$MM)), 3))
  mytable <- rbind(R2, pvalue, slope)
  if (corTest$p.value < 0.05) {
    sigcolor <- "red"
  } else {
    sigcolor <- "black"
  }
  # Generate plot.
  plot <- ggplot(data = moduleGenes, aes(x = MM, y = GS)) +
    geom_point(color = "black", pch = 21, fill = module, size = 2) +
    geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color = "black", linetype = "dashed") +
    ggtitle(paste(module, "module", "|", trait)) +
    xlab("Module Membership") +
    ylab("Gene Significance") +
    theme(
      plot.title = element_text(hjust = 0.5, color = sigcolor, size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  # Calculate ranges for positioning the annotation layer at the top right corner.
  xrange <- unlist(ggplot_build(plot)$layout$panel_params[[1]][1])
  yrange <- unlist(ggplot_build(plot)$layout$panel_params[[1]][8])
  xmin <- min(xrange)
  xmax <- max(xrange)
  xdelta <- xmax - xmin
  ymin <- min(yrange)
  ymax <- max(yrange)
  ydelta <- ymax - ymin
  tt <- ttheme_default(base_size = 11, core = list(bg_params = list(fill = "white")))
  tab <- tableGrob(mytable, rows = NULL, theme = tt)
  g <- gtable_add_grob(tab,
    grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
    t = 1, b = nrow(tab), l = 1, r = ncol(tab)
  )
  plot <- plot + annotation_custom(g,
    xmin = xmin + 0.75 * xdelta, xmax,
    ymin = ymin + 0.75 * ydelta, ymax
  )
  results <- list(plot, lmstats)
  names(results) <- c("plot", "stats")
  if (stats == TRUE) {
    return(results)
  } else {
    return(plot)
  }
}

#-------------------------------------------------------------------------------
# Function for plotting WGCNA powers.

ggplotScaleFreeFit <- function(sft) {
  # Gather the data, calculate scale free fit.
  data <- sft$fitIndices
  data$fit <- -sign(data$slope) * data$SFT.R.sq
  # Generate Scale free topology plot.
  plot1 <- ggplot(data, aes(x = Power, y = fit)) +
    geom_text(aes(label = Power), color = "red") +
    ggtitle("Scale independence") +
    xlab(expression(Soft ~ Threshold ~ Power ~ (beta))) +
    ylab(expression(Scale ~ Free ~ Topology ~ (R^2))) +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "red", size = 0.6) +
    # geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray", size = 0.6) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  # Generate mean connectivity plot.
  plot2 <- ggplot(data, aes(x = Power, y = mean.k.)) +
    geom_text(aes(label = Power), color = "red") +
    ggtitle("Mean Connectivity") +
    xlab(expression(Soft ~ Threshold ~ Power ~ (beta))) +
    ylab(expression(Mean ~ Connectivity ~ (k))) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  plot3 <- plot_grid(plot1, plot2, labels = "auto")
  data_return <- list(plot1, plot2, plot3)
  names(data_return) <- c("ScaleFreeFit", "MeanConnectivity", "Grid")
  return(data_return)
}

#-------------------------------------------------------------------------------
#' ## A Function for mixing colors.
#' Returns hex code for new color and prints plot showing the color in the console.

col2hex <- function(color, maxValue = 255) {
  z <- col2rgb(color)
  hex <- rgb(z[1], z[2], z[3], maxColorValue = maxValue)
  return(hex)
}

mixcolors <- function(color1, color2, ratio1 = 1, ratio2 = 1, plot = FALSE) {
  # Convert colors to RGB and mix (average).
  x <- col2rgb(color1)
  y <- col2rgb(color2)
  z <- (ratio1 * x + ratio2 * y) / 2
  # Convert to hex format.
  hex <- rgb(z[1], z[2], z[3], maxColorValue = 255)
  # Plot for visualizing new color.
  df <- as.data.frame(cbind(x = 1, y = 1, hex))
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_text(label = hex, colour = hex, size = 25) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  if (plot == TRUE) print(p)
  return(hex)
}

#-------------------------------------------------------------------------------
# Define function: ggplotProteinScatterPlot
ggplotProteinScatterPlot <- function(exprDat, prot1, prot2, annotate_stats = FALSE) {

  # Get data for proteins of interest.
  x <- as.numeric(exprDat[rownames(exprDat) == prot1, ])
  y <- as.numeric(exprDat[rownames(exprDat) == prot2, ])

  # Bind as df for ggplot.
  df <- as.data.frame(cbind(prot1 = x, prot2 = y))

  # Determine best fit line.
  fit <- lm(df$prot2 ~ df$prot1)

  # Calculate bicor statistics.
  stats <- bicorAndPvalue(df$prot1, df$prot2)

  # Build annotation table.
  pvalue <- paste("P-value =", formatC(stats$p, format = "e", digits = 2))
  slope <- paste("Slope =", round(as.numeric(coef(fit)[2]), 3))
  R2 <- paste("R2 =", round(stats$bicor, 3))
  mytable <- rbind(R2, pvalue, slope)
  text <- paste0("R² = ", round(stats$bicor, 2), ", P = ", formatC(stats$p, format = "e", digits = 2))
  # Generate plot with best fit line.
  plot <- ggplot(df, aes(x = prot1, y = prot2)) +
    geom_point(color = "white", pch = 21, fill = "black", size = 2) +
    geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color = "black", linetype = "dashed") +
    ggtitle(text) +
    xlab(paste("Log₂(Expression ", strsplit(prot1, "\\|")[[1]][1], ")", sep = "")) +
    ylab(paste("Log₂(Expression ", strsplit(prot2, "\\|")[[1]][1], ")", sep = "")) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )

  # Calculate ranges for positioning the annotation layer at the top left corner.
  xrange <- unlist(ggplot_build(plot)$layout$panel_params[[1]][1])
  yrange <- unlist(ggplot_build(plot)$layout$panel_params[[1]][8])
  xmin <- min(xrange)
  xmax <- max(xrange)
  xdelta <- xmax - xmin
  ymin <- min(yrange)
  ymax <- max(yrange)
  ydelta <- ymax - ymin
  tt <- ttheme_default(base_size = 11, core = list(bg_params = list(fill = "white")))
  tab <- tableGrob(mytable, rows = NULL, theme = tt)
  g <- gtable_add_grob(tab,
    grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
    t = 1, b = nrow(tab), l = 1, r = ncol(tab)
  )
  if (annotate_stats) {
    plot <- plot + annotation_custom(g,
      xmin = xmin - 0.75 * xdelta, xmax,
      ymin = ymin + 0.75 * ydelta, ymax
    )
  }

  # Add simple annotation in top left corner.
  # text <- paste0("R² = ", round(stats$bicor,2),", P = ", formatC(stats$p,format="e",digits=2))
  # plot <- plot +
  #  annotate("text", x = xmin+0.15*xdelta, y = ymax-0.05*ydelta, label = text, size = 4)

  return(plot)
}

#-------------------------------------------------------------------------------
# Function to add significance stars given a protein boxplot,
# stats with FDR column and the column to be labeled.
annotate_sig <- function(plot, stats, group, annotate = TRUE) {
  # Add symbols.
  idx <- match(plot$labels$title, rownames(stats))
  label.df <- as.data.frame(cbind(
    rownames(stats),
    stats$FDR, group,
    stats$logFC
  ))[idx, ]
  colnames(label.df) <- c("Protein", "FDR", "Group", "logFC")
  label.df$logFC <- as.numeric(label.df$logFC)
  label.df$FDR <- as.numeric(label.df$FDR)
  label.df$percentWT <- 100 * (2^label.df$logFC)
  label.df$symbol <- ""
  label.df$symbol[label.df$FDR < 0.1] <- "*"
  label.df$symbol[label.df$FDR < 0.05] <- "**"
  label.df$symbol[label.df$FDR < 0.001] <- "***"
  # Add ypos.
  groupMax <- subset(plot$data) %>%
    group_by(Group) %>%
    dplyr::summarize(max = max(Intensity))
  label.df$ypos <- 1.01 * groupMax$max[match(label.df$Group, groupMax$Group)]
  # Create annotation table.
  p.adj <- paste("P.adj =", formatC(label.df$FDR, format = "e", digits = 2))
  percentWT <- paste("Percent WT =", round(label.df$percentWT, 2))
  mytable <- rbind(p.adj, percentWT)
  # Add asterisks indicating significance.
  plot <- plot + annotate("text",
    x = label.df$Group, y = label.df$ypos,
    label = label.df$symbol, size = 8
  )

  # Calculate ranges for positioning the annotation layer at the top right corner.
  xrange <- unlist(ggplot_build(plot)$layout$panel_params[[1]][1])
  yrange <- unlist(ggplot_build(plot)$layout$panel_params[[1]][8])
  xmin <- min(xrange)
  xmax <- max(xrange)
  xdelta <- xmax - xmin
  ymin <- min(yrange)
  ymax <- max(yrange)
  ydelta <- ymax - ymin

  # Create annotation table.
  tt <- ttheme_default(base_size = 11, core = list(bg_params = list(fill = "white")))
  tab <- tableGrob(mytable, rows = NULL, theme = tt)
  g <- gtable_add_grob(tab,
    grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
    t = 1, b = nrow(tab), l = 1, r = ncol(tab)
  )
  # Add to plot.
  if (annotate == TRUE) {
    plot <- plot + annotation_custom(g,
      xmin = xmin + 0.75 * xdelta, xmax,
      ymin = ymin + 0.75 * ydelta, ymax
    )
  }
  return(plot)
}

#-------------------------------------------------------------------------------
# Function ggplotModuleSignificanceBoxplot
ggplotModuleSignificanceBoxplot <- function(x, g, trait, stats = TRUE) {

  # Bind data as data frame for ggplot.
  df <- as.data.frame(cbind(x, g))
  # df$g <- factor(df$g,levels=unique(g))

  # Order columns based on median.
  modRanks <- subset(df) %>%
    group_by(g) %>%
    dplyr::summarise(median = median(as.numeric(x)))
  modRanks <- modRanks[order(modRanks$median, decreasing = TRUE), ]
  modRanks$Order <- c(1:nrow(modRanks))
  # Sort the data
  idx <- match(df$g, modRanks$g)
  df$Rank <- modRanks$Order[idx]
  df <- df[order(df$Rank), ]
  df$g <- factor(df$g, levels = unique(df$g))

  # Calculate Kruskal Wallis pvalue.
  KWtest <- kruskal.test(as.numeric(x), as.factor(g))
  pvalue <- formatC(KWtest$p.value, format = "e", digits = 2)

  # Dunn's post-hoc test.
  Dtest <- dunnTest(as.numeric(x) ~ as.factor(g), kw = FALSE, method = "bh")$res

  # Generat boxplot.
  plot <- ggplot(df, aes(x = g, y = as.numeric(x), fill = g)) +
    geom_boxplot() +
    scale_fill_manual(values = levels(df$g)) +
    geom_point(color = "white", size = 2, pch = 21, fill = "black") +
    ggtitle(paste("Module Significance", " (", trait, ")", sep = "")) + xlab(NULL) +
    ylab("Gene Significance") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  results <- list(plot, KWtest, Dtest)
  names(results) <- c("plot", "Kruskal-Wallis", "Dunn")
  if (stats == TRUE) {
    return(results)
  } else {
    return(plot)
  }
}

#-------------------------------------------------------------------------------
# Function to add significance stars given a protein boxplot,
# stats with FDR column and the column to be labeled.
annotate_stars <- function(plot, stats) {
  data <- plot$data
  # Add symbols.
  idx <- match(plot$labels$title, rownames(stats))
  label.df <- as.data.frame(t(stats[idx, ]))
  colnames(label.df) <- "FDR"
  label.df$FDR <- as.numeric(label.df$FDR)
  label.df$symbol <- ""
  label.df$symbol[label.df$FDR < 0.1] <- "*"
  label.df$symbol[label.df$FDR < 0.05] <- "**"
  label.df$symbol[label.df$FDR < 0.001] <- "***"
  label.df$Group <- rownames(label.df)
  # Add ypos.
  label.df$ypos <- 1.01 * max(plot$data$Intensity)
  # Add asterisks indicating significance.
  plot <- plot + annotate("text",
    x = label.df$Group, y = label.df$ypos,
    label = label.df$symbol, size = 4
  )
  return(plot)
}

#-------------------------------------------------------------------------------
# Function for making pairwise contrasts.
makePairwiseContrasts <- function(g1, g2, collapse = " - ") {
  # Coerce to list if necessary.
  if (!inherits(g1, "list")) {
    g1 <- as.list(g1)
  }
  if (!inherits(g2, "list")) {
    g2 <- as.list(g2)
  }
  # Loop to generate pairwise contrasts.
  contrasts <- list()
  for (i in 1:length(g1)) {
    contrasts[[i]] <- expand.grid(g1[[i]], g2[[i]])
  }
  # Gather contrasts in a data matrix.
  contrasts <- do.call(rbind, contrasts)
  # If collapse = character, collapse.
  if (is.character(collapse)) {
    contrasts <- apply(contrasts, 1, function(x) paste(x, collapse = collapse))
  }
  return(contrasts)
}

#-------------------------------------------------------------------------------
ggranges <- function(plot) {

  # Calculate ranges.
  df <- data.frame(
    ggplot_build(plot)$layout$panel_params[[1]][1],
    ggplot_build(plot)$layout$panel_params[[1]][8]
  )
  df <- rbind(df, df[2, ] - df[1, ])
  rownames(df) <- c("min", "max", "delta")

  # TopRight
  TopRight <- data.frame(
    xmin = df$x.range[1] + 0.75 * df$x.range[3],
    xmax = df$x.range[2],
    ymin = df$y.range[1] + 0.55 * df$y.range[3],
    ymax = df$y.range[2]
  )

  # BottomRight
  BottomRight <- data.frame(
    xmin = df$x.range[1] + 0.75 * df$x.range[3],
    xmax = df$x.range[2],
    ymin = df$y.range[1] - 0.55 * df$y.range[3],
    ymax = df$y.range[2]
  )

  # Combine in a list.
  out <- list(df, TopRight, BottomRight)
  names(out) <- c("limits", "TopRight", "BottomRight")

  return(out)
}

#-------------------------------------------------------------------------------
## Function: winpath
winpath <- function() {
  x <- readClipboard()
  winpath <- gsub("\\\\", "/", x)
  return(winpath)
}

#-------------------------------------------------------------------------------
## Function: removed duplicate rows.
# Removes rows that are duplicate. e.g. A|B and B|A
removeDuplicateRows <- function(dataframe, colA = 1, colB = 2) {
  df <- dataframe[, c(colA, colB)]
  df <- as.data.frame(apply(df, 2, function(x) as.character(x)))
  df2 <- df[!duplicated(data.frame(list(do.call(pmin, df), do.call(pmax, df)))), ]
  idx <- as.numeric(rownames(df2))
  result <- dataframe[idx, ]
  return(result)
}


#-------------------------------------------------------------------------------
ggplotScaleFreePlot <- function(connectivity, nBreaks = 10, truncated = FALSE,
                                removeFirst = FALSE, main = "", ...) {
  k <- connectivity
  discretized.k <- cut(k, nBreaks)
  dk <- tapply(k, discretized.k, mean)
  p.dk <- as.vector(tapply(k, discretized.k, length) / length(k))
  breaks1 <- seq(from = min(k), to = max(k), length = nBreaks + 1)
  hist1 <- suppressWarnings(hist(k,
    breaks = breaks1, equidist = FALSE,
    plot = FALSE, right = TRUE
  )) # ...
  dk2 <- hist1$mids
  dk <- ifelse(is.na(dk), dk2, dk)
  dk <- ifelse(dk == 0, dk2, dk)
  p.dk <- ifelse(is.na(p.dk), 0, p.dk)
  log.dk <- as.vector(log10(dk))
  if (removeFirst) {
    p.dk <- p.dk[-1]
    log.dk <- log.dk[-1]
  }
  log.p.dk <- as.numeric(log10(p.dk + 1e-09))
  lm1 <- lm(log.p.dk ~ log.dk)
  pvalue <- lmp(lm1)
  print(pvalue)

  title <- paste0(
    main, " Scale Free R2 =", as.character(round(summary(lm1)$adj.r.squared, 2)),
    ", slope =", round(lm1$coefficients[[2]], 2)
  )

  OUTPUT <- data.frame(
    scaleFreeRsquared = round(summary(lm1)$adj.r.squared, 2),
    slope = round(lm1$coefficients[[2]], 2)
  )
  # Generate ggplot.
  df <- as.data.frame(cbind(log.dk, log.p.dk))
  plot <- ggplot(df, aes(x = log.dk, y = log.p.dk)) + geom_point(size = 2) +
    ggtitle(title) +
    geom_abline(intercept = coef(lm1)[1], slope = coef(lm1)[2], color = "black", linetype = "dashed") +
    labs(y = expression(Log[10](p(k)))) +
    labs(x = expression(Log[10](k))) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(hjust = 0.5, color = "black", size = 11),
      axis.title.y = element_text(hjust = 0.5, color = "black", size = 11)
    )
  out <- list(ggplot = plot, stats = OUTPUT)
  return(out)
}


#-------------------------------------------------------------------------------
# Function for visualizing GO terms.
ggplotGOscatter <- function(results_GOenrichment, color, topN = 10) {
  # Collect data in df.
  GOres <- results_GOenrichment[[color]]
  x <- GOres$enrichmentRatio
  y <- -log(GOres$pValue)
  FDR <- as.numeric(GOres$Bonferroni)
  nGenes <- GOres$nCommonGenes
  label <- GOres$shortDataSetName
  df <- data.frame(x, y, FDR, nGenes, label)
  df <- df[order(df$FDR), ]

  # Display only the topN terms.
  df$label[seq(topN + 1, nrow(df))] <- ""
  # df$label[seq(round(topN * nrow(df)), nrow(df))] <- ""

  # Generate plot.
  plot <- ggplot(df, aes(x = x, y = y, colour = FDR, size = nGenes, label = label)) +
    geom_point() + geom_text_repel(colour = "black", alpha = 0.85) +
    scale_colour_gradient(low = color, high = "white") +
    xlab("Fold Enrichment") +
    ylab("-Log(P-value)") +
    ggtitle("Go Enrichment") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 8, face = "bold"),
      axis.title.x = element_text(color = "black", size = 8, face = "bold"),
      axis.title.y = element_text(color = "black", size = 8, face = "bold")
    )
  return(plot)
}


#-------------------------------------------------------------------------------
# A Volcano plot function.
# Function for producing volcano plots.
ggplotVolcanoPlot2 <- function(df) {
  df$x <- df[, grep("FC", colnames(df))]
  df$y <- -log10(df[, grep("PValue", colnames(df))])
  logic <- df$FDR < 0.05
  df$Color[!logic] <- "gray"
  df$Color <- as.factor(df$Color)
  y_int <- -1 * log10(max(df$PValue[df$FDR < 0.05]))
  plot <- ggplot(data = df, aes(x = x, y = y, fill = Color)) +
    geom_point() + scale_color_manual(values = levels(df$Color)) +
    geom_hline(yintercept = y_int, linetype = "dashed", color = "black", size = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.6) +
    # geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 0.6) +
    xlab("Log2(Fold Change ASD vs Control)") + ylab("-Log10(P-Value)") +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold"),
      legend.position = "none"
    )
  return(plot)
}

#-----------------------------------------------------------------------------
# Function to get pvalue from a linear model object.
# Source: https://www.gettinggeneticsdone.com/2011/01/rstats-function-for-extracting-f-test-p.html

lmp <- function(modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = F)
  attributes(p) <- NULL
  return(p)
}

#-----------------------------------------------------------------------------
#' Network connectivity histogram.

ggplotHistK <- function(connectivity) {
  plot <- qplot(connectivity,
    geom = "histogram",
    binwidth = 5,
    main = "Connectivity Histogram",
    xlab = "Connectivity (k)",
    ylab = "Frequency",
    fill = I("black"),
    col = I("black"),
    alpha = 0.2
  ) +
    scale_x_continuous(limits = c(0, 50), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 11, face = "bold")
    )
  return(plot)
}

#-------------------------------------------------------------------------------
#' write.pajek
#'
#' Write network adjacency network to file in Pajek (*.net) format.
#'
#' @param adjm (matrix) symmetric adjacency matrix representing the network graph.
#' @param file (string) name of output file (e.g. 'network.net')
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{https://gephi.org/users/supported-graph-formats/pajek-net-format/}
#' @keywords network graph pajek write
#'
#' @examples
#' write.pajek(adjm, "network.net")
#'
#' @export

write.pajek <- function(adjm, file, ...) {
  # Write network adjacency matrix to .net file in Pajek format.
  # Uses data.table::fwrite for faster performance.
	require(data.table, quietly = TRUE)
	colnames(adjm) <- rownames(adjm) <- c(1:ncol(adjm))
	edge_list <- as.data.table(na.omit(melt(adjm)))
	colnames(edge_list) <- c("protA","protB","weight")
	v <- as.data.table(paste(seq(1,ncol(adjm)), " \"", seq(1,ncol(adjm)), "\"", sep = ""))
	write.table(paste("*Vertices", dim(adjm)[1]), file,
	quote = FALSE, row.names = FALSE, col.names = FALSE)
	fwrite(v, file, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE, append = TRUE)
	write.table("*Edges", file, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
	fwrite(edge_list, file, sep = " ", col.names = FALSE, append = TRUE)
}

#-------------------------------------------------------------------------------
#' silently
#'
#' suppress any unwanted output from a function with sink().
#'
#' @param func (function) symmetric adjacency matrix representing the network graph.
#' @param ... (string) additional arguments passed to func().
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords supress output silent quiet
#'
#' @examples
#' silently(wgcna::bicor, exprDat)
#'
#' @export
## Define a function that can suppress unwanted messages from a function.
silently <- function(func, ...) {
  sink(tempfile())
  out <- func(...)
  sink(NULL)
  return(out)
}

#-------------------------------------------------------------------------------
#' ggplotVerboseBoxPlot
#'
#' suppress any unwanted output from a function with sink().
#'
#' @param func (function) symmetric adjacency matrix representing the network graph.
#' @param ... (string) additional arguments passed to func().
#'
#' @return None
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#' @references \url{}
#' @keywords supress output silent quiet
#'
#' @examples
#' silently(wgcna::bicor, exprDat)
#'
#' @export
## Define a function that can suppress unwanted messages from a function.

#-------------------------------------------------------------------------------
# Function ggplotVerboseBoxplot
ggplotVerboseBoxplot <- function(x, 
				 g, 
				 levels, 
				 contrasts, 
				 color, 
				 stats = FALSE,
                                 method = "dunnett", 
				 correction_factor = 8) {
  # Bind data as data frame for ggplot.
  df <- as.data.frame(cbind(x, g))
  df$g <- factor(df$g, levels = levels)

  # Calculate Kruskal Wallis pvalue.
  KWtest <- kruskal.test(as.numeric(x), as.factor(g))
  pvalue <- round(KWtest$p.value, 3)

  # If p-value is significant, print the title in red.
  if (as.numeric(KWtest$p.value) < 0.05) {
    sigcolor <- "red"
  } else {
    sigcolor <- "black"
  }

  # Post-hoc comparisons module. Dunn or Dunnett.
  if (method == "dunn") {
    # print("Dunn's test used.")
    # Dunn's post-hoc test.
    Dtest <- dunnTest(as.numeric(x) ~ as.factor(g), kw = FALSE, method = "none")$res
    # Duplicate comparisons in reverse order. Keep rows that have comparison of interest.
    dupDtest <- Dtest
    dupDtest$Comparison <- do.call(
      rbind,
      lapply(
        strsplit(Dtest$Comparison, " - "),
        function(x) paste(c(x[2], x[1]), collapse = " - ")
      )
    )
    Dtest <- rbind(Dtest, dupDtest)
    # Keep only contrasts of interest.
    Dtest <- Dtest[Dtest$Comparison %in% contrasts, ]
  } else if (method == "dunnett") {
    # print("Dunnett's test used.")
    # Dunnett's post-hoc test for comparison to controls.
    Dtest <- DunnettTest(as.numeric(x), as.factor(g), control = c("WT.Cortex", "WT.Striatum"))
    Dtest_list <- list()
    # Extract results from PostHocTest object.
    for (i in 1:length(Dtest)) {
      Dtest_list[[names(Dtest)[i]]] <- Dtest[[i]]
    }
    Dtest <- as.data.frame(do.call(rbind, Dtest_list))
    # Reverse order of contrasts.
    Dtest <- add_column(Dtest,
      Comparison = paste(sapply(strsplit(rownames(Dtest), "-"), "[", 2),
        sapply(strsplit(rownames(Dtest), "-"), "[", 1),
        sep = " - "
      ), .before = 1
    )
    # Discard row names.
    rownames(Dtest) <- NULL

    # keep only contrasts of interest.
    Dtest <- Dtest[Dtest$Comparison %in% contrasts, ]
    # Add p.adj
    Dtest$P.adj <- p.adjust(Dtest$pval, method = "BH")
  } else {
    print("Please specify a post-hoc test (posthoc = dunn or dunnett).")
  }

  # Prepare annotation df.
  stats.df <- Dtest
  stats.df$P.adj <- as.numeric(stats.df$P.unadj) * correction_factor # BH adjustment for 8 comparisons.
  stats.df$symbol <- ""
  stats.df$symbol[stats.df$P.adj < 0.05] <- "*"
  stats.df$symbol[stats.df$P.adj < 0.01] <- "**"
  stats.df$symbol[stats.df$P.adj < 0.001] <- "***"
  stats.df$group1 <- gsub(" ", "", sapply(strsplit(stats.df$Comparison, "-"), "[", 1))
  stats.df$group2 <- gsub(" ", "", sapply(strsplit(stats.df$Comparison, "-"), "[", 2))
  stats.df$ypos <- 1.05 * max(as.numeric(df$x))

  # Generate boxplot.
  plot <- ggplot(df, aes(x = g, y = as.numeric(x), fill = g)) +
    geom_boxplot() +
    scale_fill_manual(values = rep(color, length(unique(g)))) + geom_point(color = "white", size = 1, pch = 21, fill = "black") +
    ggtitle(paste(color, " (P = ", pvalue, ")", sep = "")) + xlab(NULL) +
    ylab("Summary Expression") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
      axis.title.x = element_text(color = "black", size = 11, face = "bold"),
      axis.title.y = element_text(color = "black", size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  # Extract ggplot build information to adjust y axis.
  build <- ggplot_build(plot)
  y_lims <- build$layout$panel_scales_y[[1]]$range$range
  y_lims[2] <- 1.15 * y_lims[2]
  plot <- plot + ylim(y_lims)

  # Add asterisks indicating significance.
  plot <- plot + annotate("text",
    x = stats.df$group2, y = stats.df$ypos,
    label = stats.df$symbol, size = 6
  )

  results <- list(plot, KWtest, stats.df)
  names(results) <- c("plot", "Kruskal-Wallis", method)
  if (stats == TRUE) {
    return(results)
  } else {
    return(plot)
  }
}
