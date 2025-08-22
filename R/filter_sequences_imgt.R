#' Filter Bad Sequences
#'
#' This function filters out bad sequences from the input data frame based on the following criteria:
#' - Rows with missing CDR3 sequences (`NA`, `""`).
#' - Non-productive sequences (if a column containing "Functionality" is found, only rows with values like `TRUE` or `productive` are kept).
#' - CDR3 sequences containing invalid characters (e.g., `*`, `#`) if they are amino acid sequences (non-nucleotide).
#'
#' @param df A data frame containing repertoire data.
#' @param cdr3_col `Character`. The column name containing CDR3 sequences.
#'
#' @return A filtered data frame with the bad sequences removed and a message reporting the number of rows removed.
#'
#' @examples
#' filtered_df <- filter_sequences_imgt(df, cdr3_col = "JUNCTION")
#' @importFrom dplyr filter
#' @export
filter_sequences_imgt <- function(df, cdr3_col) {

  # Initial row count
  initial_row_count <- nrow(df)

  # Step 1: Remove rows with missing or empty CDR3 sequences
  df <- df[!is.na(df[[cdr3_col]]) & df[[cdr3_col]] != "", ]

  # Step 2: Filter out non-productive sequences based on the "Functionality" column (case-insensitive)
  functionality_col <- grep("functionality$", names(df), ignore.case = TRUE, value = TRUE)

  if (length(functionality_col) == 1) {
    productive.cols <- grep("^productive$", df[[functionality_col]], ignore.case = TRUE, value = FALSE)
    productive.cols.warnings <- grep("productive \\(see comment\\)", df[[functionality_col]], ignore.case = TRUE, value = FALSE)

    df <- df[productive.cols, ]
  }

  # Step 3: Determine if sequences are nucleotides or amino acids
  # Nucleotide sequences should only contain A, T, C, G
  nucleotide_pattern <- "^[ATCGatcg]+$"
  is_nucleotide <- grepl(nucleotide_pattern, df[[cdr3_col]])

  # Step 4: Remove rows with weird characters in CDR3 sequences if they are amino acids (non-nucleotide)
  if (!all(is_nucleotide)) {
    weird_chars <- "[*#]"
    df <- df[!grepl(weird_chars, df[[cdr3_col]]), ]
  }

  # Final row count after filtering
  final_row_count <- nrow(df)

  # Report how many rows were removed and how many rows were productive with warnings
  rows_removed <- initial_row_count - final_row_count
  message(sprintf("%d rows were removed during filtering.", rows_removed))
  message(sprintf("%d rows were productive with warnings", length(productive.cols.warnings)))

  return(df)
}
