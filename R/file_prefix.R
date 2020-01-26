#' file_prefix
#' add a prefix to filename
file_prefix <- function(output_dir,width=3) {
  indexed_files <- list.files(output_dir,pattern=("[0-9]{2,4}_"))
  if (length(indexed_files)==0) {
    last_file <- 0
  } else {
    last_file <- max(as.numeric(sapply(strsplit(indexed_files,"_"),"[",1)))
  }
  index <- formatC(last_file+1, width, format = "d", flag = "0")
  prefix <- paste(index,Sys.Date(),sep="_")
  return(prefix)
}