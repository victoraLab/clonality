#' Extract and Optionally Concatenate Columns by Pattern
#'
#' Selects columns from a data frame that match a given prefix or regular expression,
#' and optionally concatenates them into a single character vector per row.
#'
#' This function is useful for collapsing multi-chain columns (e.g. `v_gene_IGH`, `v_gene_IGK`)
#' into a unified column (e.g. `"IGH_IGK"`), or simply extracting subsets of columns for downstream processing.
#'
#' @param df A data frame containing chain-specific columns (e.g. `v_gene_IGH`, `cdr3_nt_IGK`, etc.).
#' @param column_name A string pattern. If `use_regex = FALSE`, it is treated as a prefix passed to `startsWith()`.
#'                    If `use_regex = TRUE`, it is treated as a regular expression passed to `str_detect()`.
#' @param unite Logical. If `TRUE` (default), the matched columns are collapsed row-wise with `_` as separator.
#'              If `FALSE`, the matched columns are returned as-is in a tibble.
#' @param use_regex Logical. If `TRUE`, treats `column_name` as a regular expression for column matching.
#'                  Default is `FALSE`, for prefix-based selection.
#'
#' @return If `unite = TRUE`, returns a character vector of concatenated values per row.
#'         If `unite = FALSE`, returns a tibble with the matched columns.
#'         If no columns match, a warning is issued and either `character(0)` or an empty tibble is returned.
#'
#' @examples
#' df <- tibble::tibble(
#'   v_gene_1 = c("V1", "V2"),
#'   v_gene_2 = c("V3", "V4"),
#'   other = c("x", "y")
#' )
#'
#' # Concatenate v_gene columns into a single vector
#' extract_and_unite(df, "v_gene", unite = TRUE)
#'
#' # Extract v_gene columns as a tibble
#' extract_and_unite(df, "v_gene", unite = FALSE)
#'
#' # Use regex to match any column with "gene" in the name
#' extract_and_unite(df, "gene", unite = TRUE, use_regex = TRUE)
#'
#' @export
extract_and_unite <- function(df, column_name, unite = TRUE, use_regex = FALSE) {
  df <- dplyr::ungroup(df)

  # Get matching column names safely
  matched_cols <- if (use_regex) {
    names(df)[stringr::str_detect(names(df), column_name)]
  } else {
    names(df)[startsWith(names(df), column_name)]
  }

  if (length(matched_cols) == 0) {
    warning("No columns matched: ", column_name)
    return(if (unite) character(0) else tibble::tibble())
  }

  df_selected <- df %>% dplyr::select(dplyr::any_of(matched_cols))

  if (unite) {
    df_united <- tidyr::unite(df_selected, col = !!column_name, sep = "_") %>%
      dplyr::pull({{ column_name }})
    return(df_united)
  } else return(df_selected)
}
