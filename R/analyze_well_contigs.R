#' Analyze Contigs in a Well for Quality and Background Signals
#'
#' Parses cell IDs, identifies and classifies contigs from the same well as real or background.
#' Contigs with the same clonality ID are collapsed; different clonalities in the same tree are BG1;
#' unrelated trees are BG2. Also returns filtered original input for contigs tagged as REAL.
#'
#' @param df A data frame output from `clonality()` containing at least `cell_id`, `clonality`, and `CDR3L` or `CDR3`.
#' @param cell_id_col Name of the column with cell IDs (default: "cell_id").
#' @param barcode_format Character indicating format of barcode in `cell_id`.
#' Options: `"WELL_PLATE"` (e.g., A01P09), `"PLATE_WELL"` (e.g., P09A01), or `"auto"` to detect automatically. Default: `"WELL_PLATE"`.
#'
#' @return A list with:
#'   - `collapsed`: A data frame with collapsed contigs and background tags.
#'   - `real_contigs_input`: Original input filtered to only contigs classified as REAL.
#' @export
analyze_well_contigs <- function(df, cell_id_col = "cell_id", barcode_format = "WELL_PLATE") {
  ids <- df[[cell_id_col]]

  # Attempt automatic detection
  if (barcode_format == "auto") {
    match_well_plate  <- stringr::str_match(ids, "^(.*?)([A-H]\\d{2}P\\d{2})_(\\d+)-(\\d+)$")
    match_plate_well  <- stringr::str_match(ids, "^(.*?)(P\\d{2}[A-H]\\d{2})_(\\d+)-(\\d+)$")

    n1 <- sum(!is.na(match_well_plate[, 1]))
    n2 <- sum(!is.na(match_plate_well[, 1]))

    if (n1 >= n2) {
      matches <- match_well_plate
      barcode_format <- "WELL_PLATE"
    } else {
      matches <- match_plate_well
      barcode_format <- "PLATE_WELL"
    }
  } else if (barcode_format == "WELL_PLATE") {
    matches <- stringr::str_match(ids, "^(.*?)([A-H]\\d{2}P\\d{2})_(\\d+)-(\\d+)$")
  } else if (barcode_format == "PLATE_WELL") {
    matches <- stringr::str_match(ids, "^(.*?)(P\\d{2}[A-H]\\d{2})_(\\d+)-(\\d+)$")
  } else {
    stop("Invalid 'barcode_format'. Use 'WELL_PLATE', 'PLATE_WELL', or 'auto'.")
  }

  if (anyNA(matches)) stop("Failed to parse some cell IDs. Check format or barcode_format setting.")

  df$PLATE_WELL     <- matches[, 3]
  df$PLATE          <- stringr::str_extract(df$PLATE_WELL, "P\\d{2}")
  df$WELL           <- stringr::str_extract(df$PLATE_WELL, "[A-H]\\d{2}")
  df$CONTIG_NUMBER  <- as.integer(matches[, 4])
  df$CONTIG_DEPTH   <- as.integer(matches[, 5])
  df$TREE           <- sub("\\..*", "", df$clonality)

  # Collapse by clonality
  df_summary <- df %>%
    dplyr::group_by(PLATE_WELL) %>%
    dplyr::group_modify(.f = function(d, keys) {
      d %>%
        dplyr::group_by(clonality) %>%
        dplyr::summarise(
          PLATE         = dplyr::first(PLATE),
          WELL          = dplyr::first(WELL),
          CONTIG_NUMBER = dplyr::first(CONTIG_NUMBER),
          CONTIG_DEPTH  = sum(CONTIG_DEPTH, na.rm = TRUE),
          TREE          = dplyr::first(TREE),
          collapsed     = dplyr::n() > 1,
          .groups       = "drop"
        )
    }) %>%
    dplyr::ungroup()

  # Tag REAL, BG1, BG2
  df_tagged <- df_summary %>%
    dplyr::group_by(PLATE_WELL) %>%
    dplyr::mutate(
      real_idx = which.max(CONTIG_DEPTH),
      real_tree = TREE[real_idx],
      background_tag = dplyr::case_when(
        dplyr::row_number() == real_idx ~ "REAL",
        TREE == real_tree               ~ "BG1",
        TRUE                            ~ "BG2"
      )
    ) %>%
    dplyr::select(-real_idx, -real_tree) %>%
    dplyr::ungroup()

  # Build composite key to filter original df
  df_tagged$key <- paste(df_tagged$PLATE_WELL, df_tagged$CONTIG_NUMBER)
  df$key <- paste(df$PLATE_WELL, df$CONTIG_NUMBER)

  real_keys <- df_tagged$key[df_tagged$background_tag == "REAL"]
  real_contigs_input <- df[df$key %in% real_keys, , drop = FALSE]
  real_contigs_input$key <- NULL  # clean

  df_tagged$key <- NULL  # clean

  return(list(
    collapsed = df_tagged,
    real_contigs_input = real_contigs_input
  ))
}
