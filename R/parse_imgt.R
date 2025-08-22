#' Parse IMGT Results
#'
#' This function parses the IMGT output files (Summary and Junction) and returns
#' a cleaned and annotated data frame. The function reads in IMGT summary and
#' junction data, cleans column names, and merges both datasets.
#' It also reports the number and proportion of rows containing "no results".
#'
#' @param path A character string. The directory path where the IMGT output files are located.
#' It expects "1_Summary.txt" and "6_Junction.txt" files in the provided directory.
#' @return A data frame with the parsed and annotated data from IMGT.
#' @importFrom readr read_delim problems
#' @importFrom dplyr left_join rename
#' @export
#' @examples
#' \dontrun{
#'   data <- parse_imgt("/path/to/imgt/files/")
#' }

parse_imgt <- function(path) {
  summary_path <- file.path(path, "1_Summary.txt")
  junction_path <- file.path(path, "6_Junction.txt")

  # Read the files
  summary <- suppressWarnings(
    read_delim(summary_path, delim = "\t", show_col_types = FALSE, progress = FALSE)
  )
  summary_problems <- problems(summary)

  if (nrow(summary_problems) > 0) {
    bad_rows <- unique(summary_problems$row)
    total_rows <- nrow(summary)
    percent_bad <- length(bad_rows) / total_rows * 100

    message(sprintf(
      "⚠️ IMGT Summary Warning: %d rows (~%.1f%%) appear to contain 'No results' (fewer columns).",
      length(bad_rows), percent_bad
    ))
  }

  junction <-  suppressWarnings(
    read_delim(junction_path, delim = "\t", show_col_types = FALSE)
  )


  clean_column_names <- function(df) {
    colnames(df) <- gsub(" ", "_", colnames(df))
    colnames(df) <- gsub("[^A-Za-z0-9_]", "", colnames(df))
    colnames(df) <- gsub("_$", "", colnames(df))
    return(df)
  }

  summary <- clean_column_names(summary)
  junction <- clean_column_names(junction)

  if (!"Sequence_ID" %in% colnames(summary) || !"Sequence_ID" %in% colnames(junction)) {
    stop("The required 'Sequence_ID' column is missing from one or both files.")
  }

  result <- dplyr::left_join(summary, junction[, c("Sequence_ID", "JUNCTION")], by = "Sequence_ID")
  result <- dplyr::rename(result, DNA_Junction = JUNCTION)

  return(result)
}
