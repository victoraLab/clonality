#' Simplify Gene Names Based on Cell Type and Create a Summary Report
#'
#' This helper function simplifies V and J gene names in the input data frame according to the specified cell type (T cell or B cell).
#' It also generates a report summarizing the number of functional ("F"), open reading frame ("ORF"), and pseudo ("P") genes,
#' as well as sequences with multiple annotations.
#'
#' @param data A data frame containing immune repertoire data with V and J gene columns.
#' @param cell `Character`. The cell type: "T" for T cells or "B" for B cells. The function adjusts the gene name extraction pattern accordingly.
#'
#' @return A list with:
#' - `clonal.df`: A data frame with simplified gene names and added metadata columns.
#' - `report`: A summary data frame with counts of functional types and multiple annotations.
#'
#' @importFrom stringr str_extract str_detect str_extract_all
#' @export
#'
#' @examples
#' data <- data.frame(v_genes = c("TRAV12-2*01 F", "TRBV7-9*01 F"),
#'                    j_genes = c("TRAJ33*01 F", "TRBJ2-1*01 F"))
#' result <- simplify_genenames(data, cell = "T")
#' simplified <- result$clonal.df
#' report <- result$report
simplify_genenames <- function(data, cell) {
  if (!cell %in% c("T", "B")) {
    stop("Invalid cell type. Please specify 'T' for T cells or 'B' for B cells.")
  }

  prefix <- if (cell == "T") "TR" else "IG"

  # Regex that includes the prefix in the match
  v_regex <- paste0(prefix, "[A-Z]*V[^*]+(?=\\*)")
  j_regex <- paste0(prefix, "[A-Z]*J[^*]+(?=\\*)")

  # Backup originals
  if (!"v_genes_original" %in% names(data)) data$v_genes_original <- data$v_genes
  if (!"j_genes_original" %in% names(data)) data$j_genes_original <- data$j_genes

  # Only simplify those with allele notation
  to_simplify_v <- str_detect(data$v_genes, "\\*\\d+")
  to_simplify_j <- str_detect(data$j_genes, "\\*\\d+")

  data$v_genes[to_simplify_v] <- str_extract(data$v_genes[to_simplify_v], v_regex)
  data$j_genes[to_simplify_j] <- str_extract(data$j_genes[to_simplify_j], j_regex)

  # Extract all matches from original fields (include prefix already in match)
  data$v_genes_all <- sapply(
    str_extract_all(data$v_genes_original, v_regex),
    function(x) paste(unique(x), collapse = ";")
  )

  data$j_genes_all <- sapply(
    str_extract_all(data$j_genes_original, j_regex),
    function(x) paste(unique(x), collapse = ";")
  )

  # Functional flags
  data$v_functionality <- NA
  data$v_functionality[str_detect(data$v_genes_original, "\\bF\\b")]   <- "F"
  data$v_functionality[str_detect(data$v_genes_original, "\\bORF\\b")] <- "ORF"
  data$v_functionality[str_detect(data$v_genes_original, "\\bP\\b")]   <- "P"

  data$j_functionality <- NA
  data$j_functionality[str_detect(data$j_genes_original, "\\bF\\b")]   <- "F"
  data$j_functionality[str_detect(data$j_genes_original, "\\bORF\\b")] <- "ORF"
  data$j_functionality[str_detect(data$j_genes_original, "\\bP\\b")]   <- "P"

  # Multiple annotations flag
  data$v_multiple_annotations <- str_detect(data$v_genes_original, ",\\s*or")
  data$j_multiple_annotations <- str_detect(data$j_genes_original, ",\\s*or")

  # Summary
  report <- data.frame(
    V_Functional = sum(data$v_functionality == "F", na.rm = TRUE),
    V_ORF = sum(data$v_functionality == "ORF", na.rm = TRUE),
    V_Pseudo = sum(data$v_functionality == "P", na.rm = TRUE),
    J_Functional = sum(data$j_functionality == "F", na.rm = TRUE),
    J_ORF = sum(data$j_functionality == "ORF", na.rm = TRUE),
    J_Pseudo = sum(data$j_functionality == "P", na.rm = TRUE),
    V_Multiple_Annotations = sum(data$v_multiple_annotations, na.rm = TRUE),
    J_Multiple_Annotations = sum(data$j_multiple_annotations, na.rm = TRUE)
  )

  # Failsafe: remove rows with NA in both simplified V and J
  data <- data[!(is.na(data$v_genes) & is.na(data$j_genes)), ]

  return(list(clonal.df = data, report = report))
}
