#!/usr/bin/env Rscript

load_data <- function(path, method, ...) {
  # Load data from file given a path and the method by which to load it.
  switch(method,
    RData = base::readRDS(file = path, ...),
    RDS = base::readRDS(file = path, ...),
    excel = readxl::read_excel(file = path, ...),
    csv = data.table::fread(file = path, ...)
  )
}
