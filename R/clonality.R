#' Annotate Clonality in Immune Repertoire Data
#'
#' Identifies clonal groups based on V gene, J gene, and CDR3 sequence similarity.
#' Sequences with identical V-J-CDR3 length combinations are clustered using hierarchical clustering.
#' Can return a full joined data frame or a simplified annotated output.
#'
#' When \code{cdr3_col} is a vector of length > 1, the function runs one independent
#' clustering pass per CDR3 column, each with its own \code{mismatch} and
#' \code{stringdist_method} (recycled to match).
#' The output contains one \code{clonality_<tag>} column per CDR3 column
#' (e.g. \code{clonality_nt}, \code{clonality_aa}).
#'
#' @param data A data frame or a file path to a `.xlsx` file containing repertoire data.
#' @param file `Character`. Saves output as an excel file.
#' @param ident_col `Character`. The name of the column containing unique identifiers. Default: "Sequence_ID".
#' @param vgene_col `Character`. The name of the column containing V gene annotations. Default: "VGENE_and_allele".
#' @param jgene_col `Character`. The name of the column containing J gene annotations. Default: "JGENE_and_allele".
#' @param cdr3_col `Character` vector. One or more column names containing CDR3 sequences.
#'   When multiple columns are provided, clustering is performed independently for each,
#'   producing separate \code{clonality_<tag>} columns. Default: "JUNCTION".
#' @param cell `Character`. One of "T" or "B". Specifies T cell or B cell data. Default: "T".
#' @param output_original `Logical`. If `TRUE`, joins results to the original input. Default: `FALSE`.
#' @param mismatch `Numeric` vector. Percent mismatch tolerance used to subcluster sequences within the same V/J/CDR3-length bin (range: 0--100).
#'   Internally, pairwise CDR3 string distances are normalized by CDR3 length and the hierarchical tree
#'   is cut at `hclust_height = mismatch / 100`.
#'   - A value of `0` means **no mismatch is allowed**: only identical CDR3 sequences cluster together.
#'   - A value of `10` allows clustering of sequences that differ by up to about 10% of CDR3 length.
#'   - A value of `100` is the loosest setting and effectively prevents subclustering within a duplicated V/J/length bin.
#'   When \code{cdr3_col} has length > 1, \code{mismatch} is recycled to match.
#' @param suffix `Character` vector. A string to append to the clonality group IDs. Recycled to match `cdr3_col`.
#' @param stringdist_method String distance method passed to `stringdist()`. Recycled to match `cdr3_col`.
#'   Defaults to "osa" (Optimal String Alignment). Other options include:
#'   "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard",
#'   "jw" (Jaro-Winkler), and "soundex". See `?stringdist::stringdist` for details.
#' @param cdr3_tags `Character` vector. Optional override for the CDR3 column tags used
#'   in output column names (e.g. `c("nt","aa")`). If `NULL` (default), tags are
#'   auto-detected as `"nt"` or `"aa"` based on sequence content.
#' @param filter_IMGT_sequences `Logical`. Whether to filter for productive IMGT CDR3 sequences. Default: `FALSE`.
#' @param project `Character`. Optional label to assign to all rows in the resulting data frame to identify the project or dataset. If not provided, defaults to "project".
#' @param ... Additional arguments passed to `read_repertoire()`.
#' @param quiet `Logical`. Whether to print silently the result.
#' @return A data frame with one or more `clonality_*` columns added, or written to disk if `file` is specified.
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
                      cdr3_tags = NULL,
                      ...) {

  # --- internal helpers -----------------------------------------------------
  recycle2 <- function(x, n) {
    if (length(x) == 0) return(rep(NA, n))
    if (length(x) == n) return(x)
    if (length(x) == 1) return(rep(x, n))
    stop(sprintf("Length mismatch: expected length 1 or %d, got %d.", n, length(x)))
  }

  detect_tag <- function(vec_chr) {
    # crude but robust: >80 % chars in A/C/G/T/N -> nt, else aa
    s <- paste(vec_chr, collapse = "")
    s <- gsub("[^A-Za-z]", "", s)
    if (nchar(s) == 0) return("nt")
    nt_frac <- nchar(gsub("[^ACGTNacgtn]", "", s)) / nchar(s)
    if (is.na(nt_frac)) return("nt")
    if (nt_frac >= 0.8) "nt" else "aa"
  }
  # --------------------------------------------------------------------------

  if (!is.numeric(mismatch) || any(mismatch < 0 | mismatch > 100)) {
    stop("'mismatch' must be between 0 and 100 (percent mismatch tolerance).")
  }

  # Read the input data (data frame or .xlsx file)
  imported_df <- read_repertoire(data, ...)

  missing_cols <- setdiff(c(ident_col, vgene_col, jgene_col, cdr3_col), names(imported_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required column(s): %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Drop rows missing ANY requested CDR3 column
  rm.cols <- c(NA, "")
  keep_idx <- Reduce(`&`, lapply(cdr3_col, function(cc) !(imported_df[[cc]] %in% rm.cols)))
  df <- imported_df[keep_idx, , drop = FALSE]

  if (isTRUE(filter_IMGT_sequences)) {
    df <- filter_sequences_imgt(df = df, cdr3_col = cdr3_col[[1]])
  }

  # Create a simplified data frame with V/J info
  sub.df <- data.frame(cell_id = df[[ident_col]],
                       v_genes = df[[vgene_col]],
                       j_genes = df[[jgene_col]],
                       stringsAsFactors = FALSE)

  # Check for missing V or J genes
  if (any(!stats::complete.cases(sub.df[, c("v_genes", "j_genes")]))) {
    stop("Some rows contain missing V or J genes.")
  }

  # Simplify gene names; attach CDR3 columns so they survive into clonal.df
  simple.df <- simplify_genenames(
    data = cbind(sub.df, df[, cdr3_col, drop = FALSE]),
    cell = cell
  )
  clonal.df <- simple.df[["clonal.df"]]

  # --- Per-CDR3 clustering --------------------------------------------------
  K <- length(cdr3_col)
  mismatch          <- recycle2(mismatch, K)
  stringdist_method <- recycle2(stringdist_method, K)
  suffix            <- recycle2(if (is.null(suffix)) NA_character_ else suffix, K)

  # Decide tags for output column names
  if (is.null(cdr3_tags)) {
    tags <- vapply(cdr3_col, function(cc) detect_tag(clonal.df[[cc]]), character(1))
  } else {
    if (length(cdr3_tags) == 1) tags <- rep(cdr3_tags, K) else tags <- cdr3_tags
  }
  tags <- ifelse(tags %in% c("nt", "aa"), tags, paste0("cdr3", seq_len(K)))

  out_list <- vector("list", K)

  for (k in seq_len(K)) {
    ccol <- cdr3_col[[k]]
    tag  <- tags[[k]]

    tmp <- clonal.df
    tmp$CDR3  <- clonal.df[[ccol]]
    tmp$CDR3L <- nchar(tmp$CDR3)

    # Unique bin by V/J/CDR3 length
    clonal_id <- paste(tmp[["v_genes"]], tmp[["j_genes"]], tmp[["CDR3L"]], sep = "_")

    if (any(duplicated(clonal_id))) {
      tmp <- cluster_duplicates_by_cdr3(
        clonal.df    = tmp,
        clonal_id    = clonal_id,
        cdr3_lengths = tmp$CDR3L,
        dist_method  = stringdist_method[[k]],
        h_height     = mismatch[[k]] / 100
      )
    } else {
      tmp$clonality <- assign_unique_ids(tmp)
    }

    # project label and suffix
    tmp$project <- project
    if (!is.na(suffix[[k]])) {
      tmp$clonality <- paste(suffix[[k]], tmp$clonality, sep = "_")
    }

    # keep only id + clonality for this CDR3 column
    out_list[[k]] <- tmp[, c("cell_id", "clonality"), drop = FALSE]
    # When a single CDR3 column is used, keep the simple "clonality" name
    # for backward compatibility with downstream functions (add_seurat_metadata, etc.)
    if (K > 1) {
      names(out_list[[k]])[2] <- paste0("clonality_", tag)
    }
  }

  # Merge per-CDR3 results together
  clonal.merge <- Reduce(function(x, y) dplyr::left_join(x, y, by = "cell_id"), out_list)

  # Add back base columns (V/J + all CDR3 columns) for context
  base_cols <- unique(c("cell_id", "v_genes", "j_genes", cdr3_col))
  clonal.merge <- dplyr::left_join(
    clonal.merge,
    clonal.df[, base_cols, drop = FALSE],
    by = "cell_id"
  )

  # Add project column
  clonal.merge$project <- project

  # Merge with original data or sort
  if (isTRUE(output_original)) {
    clonal.merge <- safe_full_join_coalesce(
      x = imported_df,
      y = clonal.merge,
      by = stats::setNames("cell_id", ident_col)
    )
  } else {
    clonality_cols <- names(clonal.merge)[grepl("^clonality", names(clonal.merge))]
    primary <- clonality_cols[1]
    clonal.merge <- clonal.merge[gtools::mixedorder(clonal.merge[[primary]]), ]
  }

  if (!is.null(file)) {
    openxlsx::write.xlsx(clonal.merge, file, rowNames = FALSE)
  }

  if (quiet) return(invisible(clonal.merge))
  clonal.merge
}
