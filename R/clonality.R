#' Annotate Clonality in Immune Repertoire Data
#'
#' Identifies clonal groups based on V gene, J gene, and CDR3 sequence similarity.
#' Sequences with identical V-J-CDR3 length combinations are clustered using hierarchical clustering.
#' Can return a full joined data frame or a simplified annotated output.
#'
#' @param data A data frame or a file path to a `.xlsx` file containing repertoire data.
#' @param file `Character`. Saves output as an excel file.
#' @param ident_col `Character`. The name of the column containing unique identifiers. Default: "Sequence_ID".
#' @param vgene_col `Character`. The name of the column containing V gene annotations. Default: "VGENE_and_allele".
#' @param jgene_col `Character`. The name of the column containing J gene annotations. Default: "JGENE_and_allele".
#' @param cdr3_col `Character`. The name of the column containing CDR3 sequences. Default: "JUNCTION".
#' @param cell `Character`. One of "T" or "B". Specifies T cell or B cell data. Default: "T".
#' @param output_original `Logical`. If `TRUE`, joins results to the original input. Default: `FALSE`.
#' @param mismatch `Numeric`. Percent mismatch tolerance used to subcluster sequences within the same V/J/CDR3-length bin (range: 0–100).
#'   Internally, pairwise CDR3 string distances are normalized by CDR3 length and the hierarchical tree
#'   is cut at `hclust_height = mismatch / 100`.
#'   - A value of `0` means **no mismatch is allowed**: only identical CDR3 sequences cluster together.
#'   - A value of `10` allows clustering of sequences that differ by up to about 10% of CDR3 length.
#'   - A value of `100` is the loosest setting and effectively prevents subclustering within a duplicated V/J/length bin.
#' @param suffix `Character`. A string to append to the clonality group IDs.
#' @param stringdist_method String distance method passed to `stringdist()`.
#'   Defaults to "osa" (Optimal String Alignment). Other options include:
#'   "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard",
#'   "jw" (Jaro-Winkler), and "soundex". See `?stringdist::stringdist` for details.
#' @param filter_IMGT_sequences `Logical`. Whether to filter for productive IMGT CDR3 sequences. Default: `FALSE`.
#' @param project `Character`. Optional label to assign to all rows in the resulting data frame to identify the project or dataset. If not provided, defaults to "project".
#' @param ... Additional arguments passed to `read_repertoire()`.
#' @param quiet `Logical`. Whether to print silently the result.
#' @return A data frame with a `clonality` column added, or written to disk if input is a file.
#'
#' @examples
#' results <- clonality(data = tra, vgene_col = "V_GENE_and_allele", jgene_col = "J_GENE_and_allele")
#' results <- clonality(data = trb,  vgene_col = "V_GENE_and_allele", jgene_col = "J_GENE_and_allele")
#' results <- clonality(data = IMGT_TCR_Summary_File, cdr3_col  = "DNA_Junction")
#'
#' @importFrom stringdist stringdistmatrix
#' @importFrom readxl read_excel
#' @importFrom stringr str_extract
#' @importFrom gtools mixedorder
#' @import dplyr
#' @export
clonality <- function(data = "example.xlsx",
                      ident_col = "Sequence_ID",
                      vgene_col = "VGENE_and_allele",
                      jgene_col = "JGENE_and_allele",
                      cdr3_col = "JUNCTION",
                      cell = "T",
                      output_original = FALSE,
                      mismatch = 0,
                      suffix = NULL,
                      file = NULL,
                      quiet = FALSE,
                      filter_IMGT_sequences = FALSE,
                      stringdist_method = "osa",
                      project = "project",
                      ...) {

  if (!is.numeric(mismatch) || mismatch < 0 || mismatch > 100) {
    stop("'mismatch' must be between 0 and 100 (percent mismatch tolerance).")
  }

  hclust_height <- mismatch / 100

  # Read the input data (data frame or .xlsx file)
  # If data is not a file path, assume it's a data frame from the R environment
  imported_df <- read_repertoire(data, ...)

  missing_cols <- setdiff(c(ident_col, vgene_col, jgene_col, cdr3_col), names(imported_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required column(s): %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Remove rows with missing CDR3 sequences regardless of filter_sequences option
  rm.cols <- c(NA, "")
  df <- imported_df[!imported_df[[cdr3_col]] %in% rm.cols, ]

  if(filter_IMGT_sequences == TRUE){
    df <- filter_sequences_imgt(df = df, cdr3_col = cdr3_col)
  }

  # Create a simplified data frame with relevant columns
  sub.df <- data.frame(cell_id = df[[ident_col]],
                       v_genes = df[[vgene_col]],
                       j_genes = df[[jgene_col]],
                       CDR3 = df[[cdr3_col]],
                       CDR3L = nchar(df[[cdr3_col]]),
                       stringsAsFactors = FALSE)

  # Check for missing V or J genes and stop execution if any are found
  if (any(!stats::complete.cases(sub.df))) {
    stop("Some rows contain missing V or J genes.")
  }

  # Simplify gene names using IMGT nomenclature if requested
  simple.df <- simplify_genenames(data = sub.df, cell = cell)

  # Filter by V gene, J gene, and CDR3 length to generate a unique ID for each sequence
  clonal.df <- simple.df[["clonal.df"]]
  v <- clonal.df[["v_genes"]]
  j <- clonal.df[["j_genes"]]
  l <- clonal.df[["CDR3L"]]
  clonal_id <- paste(v, j, l, sep = "_")

  # Group into clones or assign unique IDs
  if (any(duplicated(clonal_id))) {
    clonal.df <- cluster_duplicates_by_cdr3(clonal.df = clonal.df, clonal_id = clonal_id, cdr3_lengths = l, dist_method = stringdist_method, h_height = hclust_height)
  } else {
    clonal.df$clonality <- assign_unique_ids(clonal.df)
  }

  # Add project name if specified
  clonal.df$project <- project

  # Merge with original data or save a minimal version depending on user preference
  if (output_original) {
    clonal.df <- safe_full_join_coalesce(
      x = imported_df,
      y = clonal.df,
      by = stats::setNames("cell_id", ident_col)
    )
  } else {
    clonal.df <- clonal.df[mixedorder(clonal.df$clonality), ]
  }

  # Append the suffix to clonality IDs if provided
  if (!is.null(suffix)) {
    clonal.df$clonality <- paste(suffix, clonal.df$clonality, sep = "_")
  }

  if (!is.null(file)) {
    openxlsx::write.xlsx(clonal.df, file, rowNames = FALSE)
  }

  if (quiet) return(invisible(clonal.df))
  clonal.df
}
