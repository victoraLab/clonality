#' Simplify Gene Names Based on Cell Type and Create a Summary Report
#'
#' Simplifies V, D, and J gene names (when present) in the input data frame according to the specified
#' cell type (T cell or B cell). Generates a report summarizing functionality ("F", "ORF", "P")
#' and multiple-annotation flags.
#'
#' @param data A data frame containing immune repertoire data with V/J (and optionally D) gene columns.
#' @param cell `Character`. "T" for T cells or "B" for B cells.
#' @param drop_all_na_genes Logical. If TRUE, drop genes that are entirely NA after simplification.
#' @return A list with:
#' - `clonal.df`: data frame with simplified gene names and added metadata columns.
#' - `report`: summary data frame with counts of functional types and multiple annotations.
#'
#' @importFrom stringr str_extract str_detect str_extract_all
#' @export
simplify_genenames <- function(data, cell, drop_all_na_genes = FALSE) {
  if (!cell %in% c("T", "B")) {
    stop("Invalid cell type. Please specify 'T' for T cells or 'B' for B cells.")
  }

  # If someone passes a vector, make it a data.frame (nice for DGENE_and_allele-style use)
  if (!is.data.frame(data)) {
    data <- data.frame(gene = data, stringsAsFactors = FALSE)
  }

  # Detect which columns exist (support your canonical v_genes/j_genes/d_genes names)
  has_v <- "v_genes" %in% names(data)
  has_d <- "d_genes" %in% names(data)
  has_j <- "j_genes" %in% names(data)

  # Back-compat: if single column named gene, treat it as d_genes (common for DGENE_and_allele)
  if (!has_v && !has_d && !has_j && ncol(data) == 1) {
    names(data)[1] <- "d_genes"
    has_d <- TRUE
  }

  prefix <- if (cell == "T") "TR" else "IG"

  # V/J regex include prefix in match; stop before allele (*)
  v_regex <- paste0(prefix, "[A-Z]*V[^*]+(?=\\*)")
  j_regex <- paste0(prefix, "[A-Z]*J[^*]+(?=\\*)")

  # D genes:
  # - B cell heavy D genes look like IGHD...
  # - T cell delta D genes look like TRDD...
  # Use the chosen prefix + allow "D" segment explicitly.
  # For B: IGHD... ; For T: TRDD...
  d_regex <- if (cell == "B") {
    "IGHD[^*]+(?=\\*)"
  } else {
    "TRDD[^*]+(?=\\*)"
  }

  # --- Backup originals (only if the columns exist) ---
  if (has_v && !"v_genes_original" %in% names(data)) data$v_genes_original <- data$v_genes
  if (has_d && !"d_genes_original" %in% names(data)) data$d_genes_original <- data$d_genes
  if (has_j && !"j_genes_original" %in% names(data)) data$j_genes_original <- data$j_genes

  # --- Simplify V/D/J only when allele notation exists ---
  if (has_v) {
    to_simplify_v <- str_detect(data$v_genes, "\\*\\d+")
    data$v_genes[to_simplify_v] <- str_extract(data$v_genes[to_simplify_v], v_regex)

    data$v_genes_all <- sapply(
      str_extract_all(data$v_genes_original, v_regex),
      function(x) paste(unique(x), collapse = ";")
    )

    data$v_functionality <- NA_character_
    data$v_functionality[str_detect(data$v_genes_original, "\\bF\\b")]   <- "F"
    data$v_functionality[str_detect(data$v_genes_original, "\\bORF\\b")] <- "ORF"
    data$v_functionality[str_detect(data$v_genes_original, "\\bP\\b")]   <- "P"

    data$v_multiple_annotations <- str_detect(data$v_genes_original, ",\\s*or")
  }

  if (has_d) {
    to_simplify_d <- !is.na(data$d_genes) & str_detect(data$d_genes, "\\*\\d+")
    data$d_genes[to_simplify_d] <- str_extract(data$d_genes[to_simplify_d], d_regex)

    data$d_genes_all <- sapply(
      str_extract_all(data$d_genes_original, d_regex),
      function(x) paste(unique(x), collapse = ";")
    )

    data$d_functionality <- NA_character_
    data$d_functionality[str_detect(data$d_genes_original, "\\bF\\b")]   <- "F"
    data$d_functionality[str_detect(data$d_genes_original, "\\bORF\\b")] <- "ORF"
    data$d_functionality[str_detect(data$d_genes_original, "\\bP\\b")]   <- "P"

    data$d_multiple_annotations <- str_detect(data$d_genes_original, ",\\s*or")
  }

  if (has_j) {
    to_simplify_j <- str_detect(data$j_genes, "\\*\\d+")
    data$j_genes[to_simplify_j] <- str_extract(data$j_genes[to_simplify_j], j_regex)

    data$j_genes_all <- sapply(
      str_extract_all(data$j_genes_original, j_regex),
      function(x) paste(unique(x), collapse = ";")
    )

    data$j_functionality <- NA_character_
    data$j_functionality[str_detect(data$j_genes_original, "\\bF\\b")]   <- "F"
    data$j_functionality[str_detect(data$j_genes_original, "\\bORF\\b")] <- "ORF"
    data$j_functionality[str_detect(data$j_genes_original, "\\bP\\b")]   <- "P"

    data$j_multiple_annotations <- str_detect(data$j_genes_original, ",\\s*or")
  }

  # --- Summary report (only include fields that exist) ---
  report <- data.frame(dummy = 1L, stringsAsFactors = FALSE)[, 0, drop = FALSE]

  if (has_v) {
    report$V_Functional <- sum(data$v_functionality == "F", na.rm = TRUE)
    report$V_ORF <- sum(data$v_functionality == "ORF", na.rm = TRUE)
    report$V_Pseudo <- sum(data$v_functionality == "P", na.rm = TRUE)
    report$V_Multiple_Annotations <- sum(data$v_multiple_annotations, na.rm = TRUE)
  }

  if (has_d) {
    report$D_Functional <- sum(data$d_functionality == "F", na.rm = TRUE)
    report$D_ORF <- sum(data$d_functionality == "ORF", na.rm = TRUE)
    report$D_Pseudo <- sum(data$d_functionality == "P", na.rm = TRUE)
    report$D_Multiple_Annotations <- sum(data$d_multiple_annotations, na.rm = TRUE)
  }

  if (has_j) {
    report$J_Functional <- sum(data$j_functionality == "F", na.rm = TRUE)
    report$J_ORF <- sum(data$j_functionality == "ORF", na.rm = TRUE)
    report$J_Pseudo <- sum(data$j_functionality == "P", na.rm = TRUE)
    report$J_Multiple_Annotations <- sum(data$j_multiple_annotations, na.rm = TRUE)
  }

  # Failsafe: drop rows that have NA for all simplified gene fields that exist
  simplified_cols <- c()
  if (has_v) simplified_cols <- c(simplified_cols, "v_genes")
  if (has_d) simplified_cols <- c(simplified_cols, "d_genes")
  if (has_j) simplified_cols <- c(simplified_cols, "j_genes")

  if (drop_all_na_genes && length(simplified_cols) > 0) {
    all_na <- Reduce(`&`, lapply(simplified_cols, function(cc) is.na(data[[cc]])))
    data <- data[!all_na, , drop = FALSE]
  }
  list(clonal.df = data, report = report)
}
