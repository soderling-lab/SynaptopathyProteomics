#' setFont

setFont <- function(font, font_path = NULL) {

  # Use a custom font.
  suppressPackageStartupMessages({
    library(extrafont)
  })

  # MY FONT PATH.
  if (is.null(font_path)) {
    font_path <- file.path(
      "/mnt/c/Program Files/Adobe",
      "Adobe Illustrator 2020",
      "Support Files/Required/PDFL Resource",
      "Resource/Fonts/TTF/"
    )
  }

  # Check if font is in available fonts.
  if (!(font %in% fonts())) {
    suppressWarnings({
      font_import(paths = font_path, prompt = FALSE)
    })
  }

  # Check again.
  if (!(font %in% fonts())) {
    stop(paste("Font not found in path."))
  }

  # Set current font.
  invisible(choose_font(font))
}
