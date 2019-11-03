# Split functions into seperate files.

here <- getwd()
root <- dirname(here)

all_functions <- readLines("all_fun.R")
blank_header <- readLines("blank-header.txt")

# Split at '#----'
idx <- grep("#----------", all_functions)

# Loop to collect functions into a named list.
# If a function is duplicated, it will only appear in the list once!
my_functions <- list()
for (i in seq_along(idx)) {
  if (i == 1) {
    x0 <- 1
  } else {
    x0 <- idx[i - 1]
  }
  x1 <- idx[i]
  txt <- all_functions[c(x0:x1)]
  namen <- strsplit(txt[grep(" <- function", txt)], "\\ |<-")[[1]][1]
  my_functions[[namen]] <- txt
}

# Loop to save as seperate files with blank header.
for (i in 1:length(my_functions)) {
  namen <- names(my_functions[i])
  myfile <- paste0(namen, ".R")
  myheader <- blank_header
  myheader[1] <- gsub("function_name", namen, blank_header[1])
  x <- my_functions[[i]]
  x <- c(myheader, x[2:length(x) - 1])
  writeLines(x, myfile)
}
