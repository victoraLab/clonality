#' Parse and Validate Cell Ranger Data
#'
#' Accepts a data frame, a variable name as a character string, or a file path to a `.csv` file.
#' Returns a parsed data frame and (optionally) validates the number of columns according to Cell Ranger version.
#'
#' @param data Either a data frame, the name of a data object, or a file path to a `.csv` file.
#' @param cellranger_version Integer. Supported: 3, 7, or 9. Used to validate the expected number of columns. Default is 7.
#' @param ignore_version Logical. If `TRUE`, skip version-based column-count validation. Default is `FALSE`.
#'
#' @return A (optionally validated) data frame with structure compatible with the specified Cell Ranger version.
#'
#' @examples
#' parsed <- cellranger_version_parser(Cellranger3_TCR, cellranger_version = 3)
#' parsed <- cellranger_version_parser("Cellranger3_TCR", cellranger_version = 3)
#' parsed <- cellranger_version_parser("path/to/file.csv", ignore_version = TRUE)
#'
#' @export
cellranger_version_parser <- function(data, cellranger_version = 7, ignore_version = FALSE) {
  # If data is not a data.frame, determine its type
  if (!is.data.frame(data)) {
    if (is.character(data)) {
      ext <- tools::file_ext(data)
      if (tolower(ext) == "csv") {
        data <- read.csv(data, stringsAsFactors = FALSE)
      } else {
        if (exists(data, envir = .GlobalEnv)) {
          data <- get(data, envir = .GlobalEnv)
        } else {
          stop("Could not interpret 'data'. It is neither a valid file path nor an object name.")
        }
      }
    } else {
      # Non-data.frame, non-character input
      if (exists(deparse(substitute(data)), envir = .GlobalEnv)) {
        data <- get(deparse(substitute(data)), envir = .GlobalEnv)
      } else {
        stop("Unsupported 'data' type. Provide a data.frame, object name (character), or CSV path.")
      }
    }
  }

  # Optionally validate number of columns by version
  if (!isTRUE(ignore_version)) {
    col_count <- ncol(data)

    if (cellranger_version == 3 && col_count != 18) {
      stop("Cell Ranger v3 is expected to have 18 columns, but found ", col_count, ".")
    }

    if (cellranger_version == 7 && col_count != 31) {
      stop("Cell Ranger v7 is expected to have 31 columns, but found ", col_count, ".")
    }

    if (cellranger_version == 9 && col_count != 32) {
      stop("Cell Ranger v9 is expected to have 32 columns, but found ", col_count, ".")
    }
  }

  return(data)
}
