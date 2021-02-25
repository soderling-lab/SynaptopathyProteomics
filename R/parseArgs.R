parseArgs <- function(tissue = commandArgs(trailingOnly = TRUE)) {
  # parse the command line arguments
  msg <- c(
    "Please specify a tissue type to be analyzed:\n",
    "Choose either 'Cortex' or 'Striatum'."
  )
  # if interactive, return a random tissue
  if (interactive()) {
    return(sample(c("Cortex", "Striatum"), 1))
  } else {
    switch(tissue,
      Cortex = return(tissue),
      Striatum = return(tissue)
    )
  }
} #EOF
