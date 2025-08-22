#' Parse General Table Files
#'
#' This function parses general table files such as `.xlsx`, `.csv`, or `.txt` into a tibble or data frame.
#' If the input is already a data frame, it returns it directly.
#'
#' @param data A data frame object or the full path to a `.xlsx`, `.csv`, or `.txt` file.
#' The file extension is used to determine how to read the file.
#' @return A tibble or data frame containing the parsed data from the input file.
#' @importFrom readr read_csv read_delim
#' @importFrom readxl read_excel
#' @importFrom stringr str_extract
#' @importFrom tibble as_tibble
#' @export
#' @examples
#' \dontrun{
#'   df <- parse_general("data.csv")   # To parse a CSV file
#'   df <- parse_general("data.xlsx")  # To parse an Excel file
#'   df <- parse_general("data.txt")   # To parse a tab-delimited text file
#'   df <- parse_general(existing_df)  # To return an existing data frame as is
#' }

parse_general <- function(data) {
  # If input is a data frame, return it directly
  if (is.data.frame(data)) {
    return(data)
  }

  # If input is a character, assume it's a file path
  if (is.character(data)) {

    # Extract the file extension
    ext <- stringr::str_extract(data, "\\.[a-zA-Z0-9]+$")

    # Read file based on extension
    if (ext == ".xlsx") {
      imported_df <- readxl::read_excel(data)
    } else if (ext == ".csv") {
      imported_df <- readr::read_csv(data)
    } else if (ext == ".txt") {
      imported_df <- readr::read_delim(data, delim = "\t")
    } else {
      stop("Unsupported file format. Please provide a '.csv', '.xlsx', or '.txt' file.", call. = FALSE)
    }

    return(tibble::as_tibble(imported_df))

  } else {
    stop("Input must be a data frame or a file path to a '.csv', '.xlsx', or '.txt' file.", call. = FALSE)
  }
}
