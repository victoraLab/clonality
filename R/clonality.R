#' Annotate Clonality in Immune Repertoire Data
#'
#' Identifies clonal groups based on V gene, J gene, and CDR3 sequence similarity.
#' Sequences with identical V-J-CDR3 length combinations are clustered using hierarchical clustering.
#' Can return a full joined data frame or a simplified annotated output.
#'
#' @param data A data frame or a file path to a `.xlsx` file containing repertoire data.
#' @param output `Character`. The name of the output variable or file. Default: "output".
#' @param ident_col `Character`. The name of the column containing unique identifiers. Default: "Sequence_ID".
#' @param vgene_col `Character`. The name of the column containing V gene annotations. Default: "VGENE_and_allele".
#' @param jgene_col `Character`. The name of the column containing J gene annotations. Default: "JGENE_and_allele".
#' @param cdr3_col `Character`. The name of the column containing CDR3 sequences. Default: "JUNCTION".
#' @param cell `Character`. One of "T" or "B". Specifies T cell or B cell data. Default: "T".
#' @param output_original `Logical`. If `TRUE`, joins results to the original input. Default: `FALSE`.
#' @param mismatch `Numeric`. Percent similarity threshold used to define clonality clusters (range: 0–100).
#'   This defines the maximum allowable difference (as a percentage of CDR3 sequence length) within a clone.
#'   - A value of `0` means **no differences allowed** (strict clustering: only identical sequences cluster).
#'   - A value of `100` means **any difference is allowed** (loose clustering: sequences are not subclustered).
#'   Internally used to calculate clustering height: `hclust_height = mismatch / 100`.
#' @param suffix `Character`. A string to append to the clonality group IDs.
#' @param stringdist_method String distance method passed to `stringdist()`.
#'   Defaults to "osa" (Optimal String Alignment). Other options include:
#'   "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard",
#'   "jw" (Jaro-Winkler), and "soundex". See `?stringdist::stringdist` for details.
#' @param filter_IMGT_sequences `Logical`. Whether to filter for productive IMGT CDR3 sequences. Default: `FALSE`.
#' @param project `Character`. Optional label to assign to all rows in the resulting data frame to identify the project or dataset. If not provided, defaults to "project".
#' @param ... Additional arguments passed to `read_repertoire()`.
#'
#' @return A data frame with a `clonality` column added, or written to disk if input is a file.
#'
#' @examples
#' clonality(data = tra, vgene_col = "V_GENE_and_allele", jgene_col = "J_GENE_and_allele")
#' clonality(data = trb, output = "clones", vgene_col = "V_GENE_and_allele", jgene_col = "J_GENE_and_allele")
#' clonality(data = IMGT_TCR_Summary_File, cdr3_col  = "DNA_Junction")
#'
#' @importFrom stringdist stringdistmatrix
#' @importFrom readxl read_excel
#' @importFrom stringr str_extract
#' @importFrom safejoin safe_full_join
#' @importFrom gtools mixedorder
#' @import dplyr
#' @export
clonality <- function(data = "example.xlsx",
                      output = "output",
                      ident_col = "Sequence_ID",
                      vgene_col = "VGENE_and_allele",
                      jgene_col = "JGENE_and_allele",
                      cdr3_col = "JUNCTION",                 # now accepts length 1 or 2
                      cell = "T",
                      output_original = FALSE,
                      mismatch = 0,
                      suffix = NULL,
                      filter_IMGT_sequences = FALSE,
                      stringdist_method = "osa",
                      project = "project",
                      cdr3_tags = NULL,                      # optional override for column tags, e.g. c("nt","aa")
                      ...) {

  # --- helpers ---------------------------------------------------------------
  recycle2 <- function(x, n) {
    if (length(x) == 0) return(rep(NA, n))
    if (length(x) == n) return(x)
    if (length(x) == 1) return(rep(x, n))
    stop(sprintf("Length mismatch: expected length 1 or %d, got %d.", n, length(x)))
  }

  detect_tag <- function(vec_chr) {
    # crude but robust: >80% chars in A/C/G/T/N -> nt, else aa
    s <- paste(vec_chr, collapse = "")
    s <- gsub("[^A-Za-z]", "", s)
    if (nchar(s) == 0) return("nt")
    nt_frac <- nchar(gsub("[^ACGTNacgtn]", "", s)) / nchar(s)
    if (is.na(nt_frac)) return("nt")
    if (nt_frac >= 0.8) "nt" else "aa"
  }

  # --------------------------------------------------------------------------
  if (!is.numeric(mismatch) || any(mismatch < 0 | mismatch > 100)) {
    stop("'mismatch' must be between 0 and 100 (percent similarity threshold).")
  }

  # Read input
  imported_df <- read_repertoire(data, ...)
  missing_cols <- setdiff(c(ident_col, vgene_col, jgene_col, cdr3_col), names(imported_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required column(s): %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Always drop empty CDR3 rows (across any requested CDR3 col)
  rm.cols <- c(NA, "")
  keep_idx <- Reduce(`&`, lapply(cdr3_col, function(cc) !(imported_df[[cc]] %in% rm.cols)))
  df <- imported_df[keep_idx, , drop = FALSE]

  if (isTRUE(filter_IMGT_sequences)) {
    # you may keep your existing filter; it will run on the primary column
    df <- filter_sequences_imgt(df = df, cdr3_col = cdr3_col[[1]])
  }

  # reduced frame
  sub.df <- data.frame(
    cell_id = df[[ident_col]],
    v_genes = df[[vgene_col]],
    j_genes = df[[jgene_col]],
    stringsAsFactors = FALSE
  )

  # sanity on V/J
  if (any(!complete.cases(sub.df[, c("v_genes","j_genes")]))) {
    stop("Some rows contain missing V or J genes.")
  }

  # simplify names
  simple.df <- simplify_genenames(data = cbind(sub.df, df[, cdr3_col, drop = FALSE]), cell = cell)
  clonal.df  <- simple.df[["clonal.df"]]

  # Prepare per-CDR3 configs
  K <- length(cdr3_col)
  mismatch          <- recycle2(mismatch, K)
  stringdist_method <- recycle2(stringdist_method, K)
  suffix            <- recycle2(suffix %||% NA_character_, K)

  # decide tags
  if (is.null(cdr3_tags)) {
    tags <- vapply(cdr3_col, function(cc) detect_tag(clonal.df[[cc]]), character(1))
  } else {
    if (length(cdr3_tags) == 1) tags <- rep(cdr3_tags, K) else tags <- cdr3_tags
  }
  tags <- ifelse(tags %in% c("nt","aa"), tags, paste0("cdr3", seq_len(K)))

  # container for outputs
  out_list <- vector("list", K)

  # run once per CDR3 column
  for (k in seq_len(K)) {
    ccol <- cdr3_col[[k]]
    tag  <- tags[[k]]

    tmp <- clonal.df
    tmp$CDR3  <- clonal.df[[ccol]]
    tmp$CDR3L <- nchar(tmp$CDR3)

    # Unique bin by V/J/length
    clonal_id <- paste(tmp[["v_genes"]], tmp[["j_genes"]], tmp[["CDR3L"]], sep = "_")

    if (any(duplicated(clonal_id))) {
      tmp <- cluster_duplicates_by_cdr3(
        clonal.df      = tmp,
        clonal_id      = clonal_id,
        cdr3_lengths   = tmp$CDR3L,
        dist_method    = stringdist_method[[k]],
        h_height       = mismatch[[k]] / 100
      )
    } else {
      tmp$clonality <- assign_unique_ids(tmp)
    }

    # project, suffixing per CDR3 if requested
    tmp$project <- project
    if (!is.na(suffix[[k]])) {
      tmp$clonality <- paste(suffix[[k]], tmp$clonality, sep = "_")
    }

    # keep only id + clonality
    out_list[[k]] <- tmp[, c("cell_id","clonality"), drop = FALSE]
    names(out_list[[k]])[2] <- paste0("clonality_", tag)
  }

  # merge per-CDR3 results together
  clonal.merge <- Reduce(function(x, y) dplyr::left_join(x, y, by = "cell_id"), out_list)

  # add back basics (V/J + the provided CDR3 columns) for context
  base_cols <- c("cell_id","v_genes","j_genes", cdr3_col)
  clonal.merge <- dplyr::left_join(
    clonal.merge,
    clonal.df[, unique(base_cols), drop = FALSE],
    by = "cell_id"
  )

  # merge with original or order
  if (isTRUE(output_original)) {
    clonal.merge <- safe_full_join(
      x = imported_df, y = clonal.merge,
      by = setNames("cell_id", ident_col), conflict = coalesce
    )
  } else {
    primary <- names(clonal.merge)[grepl("^clonality_", names(clonal.merge))][1]
    clonal.merge <- clonal.merge[gtools::mixedorder(clonal.merge[[primary]]), ]
  }

  # write or assign
  if (is.character(data)) {
    openxlsx::write.xlsx(clonal.merge, sprintf("%s.xlsx", output), row.names = FALSE)
  } else {
    assign(output, clonal.merge, envir = .GlobalEnv)
  }
}
