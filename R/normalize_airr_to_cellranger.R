#' Normalize AIRR / IgBLAST-style contigs to a CellRanger-shape table
#'
#' IgBLAST AIRR output and other long-read pipelines (e.g. nanopore + BLAZE)
#' encode information slightly differently from CellRanger:
#' - V/J/C calls may be comma-separated lists of alleles (`IGHV1-72*01,IGHV1-72*02`)
#'   while CellRanger reports a single allele.
#' - C-gene calls are usually absent and must be supplied externally
#'   (e.g. via minimap2 against constant-region references).
#' - The `productive` flag is sometimes a logical, sometimes the strings
#'   `"T"`/`"F"`, while [tenx()] expects `"true"`/`"false"`.
#'
#' This helper applies those normalizations in place. It does **not** rename
#' columns; pass the result to [tenx()] which calls `normalize_tenx_columns()`
#' to map AIRR aliases (`v_call`, `junction`, ...) to CellRanger names.
#'
#' @param data A data frame with AIRR-style columns.
#' @param c_gene_col Optional. Name of an external C-gene column (for example
#'   from minimap2). If supplied and `c_call`/`c_gene` is missing or empty,
#'   the values are copied across.
#' @param collapse_alleles Logical. If `TRUE` (default), keep only the first
#'   allele when V/J/C calls contain comma-separated lists.
#' @param productive_to_logical Logical. If `TRUE` (default), coerce the
#'   `productive` column to lowercase strings `"true"` / `"false"`.
#' @param allele_cols Character. Columns whose comma-separated allele lists
#'   should be collapsed. Default: `c("v_call", "d_call", "j_call", "c_call")`,
#'   plus already-renamed CellRanger equivalents `c("v_gene", "d_gene", "j_gene", "c_gene")`.
#'
#' @return The data frame with normalized columns.
#'
#' @examples
#' \dontrun{
#' airr <- readr::read_tsv("igblast_airr.tsv")
#' airr_norm <- normalize_airr_to_cellranger(airr, c_gene_col = "isotype_minimap2")
#' out <- tenx(airr_norm, method = "unique_paired", cell = "B")
#' }
#'
#' @export
normalize_airr_to_cellranger <- function(data,
                                         c_gene_col = NULL,
                                         collapse_alleles = TRUE,
                                         productive_to_logical = TRUE,
                                         allele_cols = c("v_call", "d_call", "j_call", "c_call",
                                                         "v_gene", "d_gene", "j_gene", "c_gene")) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  # ── 1. Collapse comma-separated alleles ──────────────────────────────────
  if (collapse_alleles) {
    for (col in intersect(allele_cols, names(data))) {
      vals <- as.character(data[[col]])
      has_comma <- grepl(",", vals, fixed = TRUE)
      if (any(has_comma, na.rm = TRUE)) {
        vals[has_comma] <- sub(",.*$", "", vals[has_comma])
        data[[col]] <- vals
      }
    }
  }

  # ── 2. Fill C-gene from an external column ───────────────────────────────
  if (!is.null(c_gene_col)) {
    if (!c_gene_col %in% names(data)) {
      stop(sprintf("`c_gene_col` = '%s' not found in data.", c_gene_col), call. = FALSE)
    }
    target <- if ("c_gene" %in% names(data)) "c_gene" else if ("c_call" %in% names(data)) "c_call" else "c_gene"
    if (!target %in% names(data)) data[[target]] <- NA_character_
    cur <- as.character(data[[target]])
    ext <- as.character(data[[c_gene_col]])
    fill <- is.na(cur) | !nzchar(cur)
    cur[fill] <- ext[fill]
    data[[target]] <- cur
  }

  # ── 3. Coerce productive to "true"/"false" ───────────────────────────────
  if (productive_to_logical && "productive" %in% names(data)) {
    p <- data[["productive"]]
    if (is.logical(p)) {
      p <- ifelse(is.na(p), NA_character_, ifelse(p, "true", "false"))
    } else {
      p <- as.character(p)
      lp <- tolower(p)
      lp[lp %in% c("t", "true", "1", "yes", "y")]   <- "true"
      lp[lp %in% c("f", "false", "0", "no", "n")]   <- "false"
      lp[!lp %in% c("true", "false")]               <- NA_character_
      p <- lp
    }
    data[["productive"]] <- p
  }

  data
}
