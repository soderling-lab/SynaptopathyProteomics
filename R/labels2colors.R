#' labels2colors
#'
#' convert labels to colors
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
#' labels2colors(labels)
labels2colors <- function(labels, zeroIsGrey = TRUE, colorSeq = NULL, naColor = "grey", commonColorCode = TRUE) # small function from WGCNA
{
  if (is.null(colorSeq)) {
    colorSeq <-
      c(
        "turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink", "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue", "lightcyan", "grey60", "lightgreen", "lightyellow", "royalblue",
        "darkred", "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white", "skyblue", "saddlebrown", "steelblue", "paleturquoise", "violet", "darkolivegreen", "darkmagenta", "sienna3", "yellowgreen",
        "skyblue3", "plum1", "orangered4", "mediumpurple3", "lightsteelblue1", "lightcyan1", "ivory", "floralwhite", "darkorange2", "brown4", "bisque4", "darkslateblue", "plum2", "thistle2", "thistle1", "salmon4",
        "palevioletred3", "navajowhite2", "maroon", "lightpink4", "lavenderblush3", "honeydew1", "darkseagreen4", "coral1", "antiquewhite4", "coral2", "mediumorchid", "skyblue2", "yellow4", "skyblue1", "plum", "orangered3",
        "mediumpurple2", "lightsteelblue", "lightcoral", "indianred4", "firebrick4", "darkolivegreen4", "brown2", "blue2", "darkviolet", "plum3", "thistle3", "thistle", "salmon2", "palevioletred2", "navajowhite1", "magenta4",
        "lightpink3", "lavenderblush2", "honeydew", "darkseagreen3", "coral", "antiquewhite2", "coral3", "mediumpurple4", "skyblue4", "yellow3", "sienna4", "pink4", "orangered1", "mediumpurple1", "lightslateblue",
        "lightblue4", "indianred3", "firebrick3", "darkolivegreen2", "blueviolet", "blue4", "deeppink", "plum4", "thistle4", "tan4", "salmon1", "palevioletred1", "navajowhite", "magenta3", "lightpink2", "lavenderblush1",
        "green4", "darkseagreen2", "chocolate4", "antiquewhite1", "coral4", "mistyrose", "slateblue", "yellow2", "sienna2", "pink3", "orangered", "mediumpurple", "lightskyblue4", "lightblue3", "indianred2", "firebrick2",
        "darkolivegreen1", "blue3", "brown1", "deeppink1", "powderblue", "tomato", "tan3", "royalblue3", "palevioletred", "moccasin", "magenta2", "lightpink1", "lavenderblush", "green3", "darkseagreen1", "chocolate3",
        "aliceblue", "cornflowerblue", "navajowhite3", "slateblue1", "whitesmoke", "sienna1", "pink2", "orange4", "mediumorchid4", "lightskyblue3", "lightblue2", "indianred1", "firebrick", "darkgoldenrod4", "blue1",
        "brown3", "deeppink2", "purple2", "tomato2", "tan2", "royalblue2", "paleturquoise4", "mistyrose4", "magenta1", "lightpink", "lavender", "green2", "darkseagreen", "chocolate2", "antiquewhite", "cornsilk",
        "navajowhite4", "slateblue2", "wheat3", "sienna", "pink1", "orange3", "mediumorchid3", "lightskyblue2", "lightblue1", "indianred", "dodgerblue4", "darkgoldenrod3", "blanchedalmond", "burlywood", "deepskyblue", "red1",
        "tomato4", "tan1", "rosybrown4", "paleturquoise3", "mistyrose3", "linen", "lightgoldenrodyellow", "khaki4", "green1", "darksalmon", "chocolate1", "antiquewhite3", "cornsilk2", "oldlace", "slateblue3", "wheat1",
        "seashell4", "peru", "orange2", "mediumorchid2", "lightskyblue1", "lightblue", "hotpink4", "dodgerblue3", "darkgoldenrod1", "bisque3", "burlywood1", "deepskyblue4", "red4", "turquoise2", "steelblue4", "rosybrown3",
        "paleturquoise1", "mistyrose2", "limegreen", "lightgoldenrod4", "khaki3", "goldenrod4", "darkorchid4", "chocolate", "aquamarine", "cyan1", "orange1", "slateblue4", "violetred4", "seashell3", "peachpuff4",
        "olivedrab4", "mediumorchid1", "lightskyblue", "lemonchiffon4", "hotpink3", "dodgerblue1", "darkgoldenrod", "bisque2", "burlywood2", "dodgerblue2", "rosybrown2", "turquoise4", "steelblue3", "rosybrown1",
        "palegreen4", "mistyrose1", "lightyellow4", "lightgoldenrod3", "khaki2", "goldenrod3", "darkorchid3", "chartreuse4", "aquamarine1", "cyan4", "orangered2", "snow", "violetred2", "seashell2", "peachpuff3",
        "olivedrab3", "mediumblue", "lightseagreen", "lemonchiffon3", "hotpink2", "dodgerblue", "darkblue", "bisque1", "burlywood3", "firebrick1", "royalblue1", "violetred1", "steelblue1", "rosybrown", "palegreen3",
        "mintcream", "lightyellow3", "lightgoldenrod2", "khaki1", "goldenrod2", "darkorchid2", "chartreuse3", "aquamarine2", "darkcyan", "orchid", "snow2", "violetred", "seashell1", "peachpuff2", "olivedrab2",
        "mediumaquamarine", "lightsalmon4", "lemonchiffon2", "hotpink1", "deepskyblue3", "cyan3", "bisque", "burlywood4", "forestgreen", "royalblue4", "violetred3", "springgreen3", "red3", "palegreen1", "mediumvioletred",
        "lightyellow2", "lightgoldenrod1", "khaki", "goldenrod1", "darkorchid1", "chartreuse2", "aquamarine3", "darkgoldenrod2", "orchid1", "snow4", "turquoise3", "seashell", "peachpuff1", "olivedrab1", "maroon4",
        "lightsalmon3", "lemonchiffon1", "hotpink", "deepskyblue2", "cyan2", "beige", "cadetblue", "gainsboro", "salmon3", "wheat", "springgreen2", "red2", "palegreen", "mediumturquoise", "lightyellow1", "lightgoldenrod",
        "ivory4", "goldenrod", "darkorchid", "chartreuse1", "aquamarine4", "darkkhaki", "orchid3", "springgreen1", "turquoise1", "seagreen4", "peachpuff", "olivedrab", "maroon3", "lightsalmon2", "lemonchiffon", "honeydew4",
        "deepskyblue1", "cornsilk4", "azure4", "cadetblue1", "ghostwhite", "sandybrown", "wheat2", "springgreen", "purple4", "palegoldenrod", "mediumspringgreen", "lightsteelblue4", "lightcyan4", "ivory3", "gold3",
        "darkorange4", "chartreuse", "azure", "darkolivegreen3", "palegreen2", "springgreen4", "tomato3", "seagreen3", "papayawhip", "navyblue", "maroon2", "lightsalmon1", "lawngreen", "honeydew3", "deeppink4", "cornsilk3",
        "azure3", "cadetblue2", "gold", "seagreen", "wheat4", "snow3", "purple3", "orchid4", "mediumslateblue", "lightsteelblue3", "lightcyan3", "ivory2", "gold2", "darkorange3", "cadetblue4", "azure1", "darkorange1",
        "paleturquoise2", "steelblue2", "tomato1", "seagreen2", "palevioletred4", "navy", "maroon1", "lightsalmon", "lavenderblush4", "honeydew2", "deeppink3", "cornsilk1", "azure2", "cadetblue3", "gold4", "seagreen1",
        "yellow1", "snow1", "purple1", "orchid2", "mediumseagreen", "lightsteelblue2", "lightcyan2", "ivory1", "gold1"
      )
  } # WGCNA ordered standardColors()

  if (is.numeric(labels)) {
    if (zeroIsGrey) {
      minLabel <- 0
    } else {
      minLabel <- 1
    }
    if (any(labels < 0, na.rm = TRUE)) minLabel <- min(c(labels), na.rm = TRUE)
    nLabels <- labels
  }
  else {
    if (commonColorCode) {
      factors <- factor(c(as.matrix(as.data.frame(labels))))
      nLabels <- as.numeric(factors)
      dim(nLabels) <- dim(labels)
    }
    else {
      labels <- as.matrix(as.data.frame(labels))
      factors <- list()
      for (c in 1:ncol(labels)) factors[[c]] <- factor(labels[, c])
      nLabels <- sapply(factors, as.numeric)
    }
  }
  if (max(nLabels, na.rm = TRUE) > length(colorSeq)) {
    nRepeats <- as.integer((max(labels) - 1) / length(colorSeq)) + 1
    warning(paste0("Number of labels exceeds number of available colors. Some colors will be repeated ", nRepeats, " times."))
    extColorSeq <- colorSeq
    for (rep in 1:nRepeats) extColorSeq <- c(extColorSeq, paste(colorSeq, ".", rep, sep = ""))
  }
  else {
    nRepeats <- 1
    extColorSeq <- colorSeq
  }
  colors <- rep("grey", length(nLabels))
  fin <- !is.na(nLabels)
  colors[!fin] <- naColor
  finLabels <- nLabels[fin]
  colors[fin][finLabels != 0] <- extColorSeq[finLabels[finLabels != 0]]
  if (!is.null(dim(labels))) dim(colors) <- dim(labels)
  colors
}
