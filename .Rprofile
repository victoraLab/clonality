.First <- function() {
  packages <- c("tidyverse", "data.table", "ggplot2", "devtools", "roxygen2", "usethis", "testthat")
  sapply(packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Package", pkg, "not found."))
    } else {
      library(pkg, character.only = TRUE)
    }
  })
}
