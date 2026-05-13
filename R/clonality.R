#' Annotate clonality in immune-repertoire data
#'
#' Groups sequences into clones based on V/J gene identity and CDR3 similarity.
#' Sequences sharing V, J, and CDR3 length are clustered by CDR3 string distance
#' and the tree is cut at `mismatch / 100`.
#'
#' `cdr3_col` and `mismatch` accept vectors. Whichever is longer determines how
#' many independent clustering passes are run; the shorter argument is recycled
#' (must be length 1). Each pass produces a `clonality_<tag>` column.
#'
#' - `cdr3_col = "Junction", mismatch = 0` → one `clonality` column (legacy)
#' - `cdr3_col = c("nt","aa"), mismatch = c(0, 10)` → `clonality_nt`, `clonality_aa`
#' - `cdr3_col = "Junction", mismatch = c(0, 10, 20)` → `clonality_mm0`, `clonality_mm10`, `clonality_mm20`
#'
#' @param data Data frame or path to `.xlsx`.
#' @param ident_col,vgene_col,jgene_col Column names for the cell ID, V gene, J gene.
#' @param cdr3_col Character vector of one or more CDR3 column names.
#' @param mismatch Numeric vector of percent mismatch tolerances (0–100). Cut height = `mismatch / 100`. `0` = identical CDR3 only; `10` ≈ 10% divergence allowed.
#' @param cell `"T"` or `"B"`.
#' @param output_original If `TRUE`, joins the clonality column(s) back onto the full input.
#' @param file Optional `.xlsx` path to save the result.
#' @param suffix,stringdist_method Recycled to match the number of passes.
#' @param cdr3_tags Optional override for the per-pass tags used in `clonality_<tag>` column names. Auto-detected as `"nt"`/`"aa"` from sequence content; when the same `cdr3_col` is reused across passes, tags fall back to `"mm<value>"`.
#' @param extra_bin_cols Optional character vector of column names whose values are
#'   appended to the V/J/length bin key. Cells must match on all of these columns
#'   to be considered for sub-clustering by CDR3 distance. Used by [assigntenx()]
#'   to enforce per-chain CDR3 length matching (so that paired cells with
#'   compensating indels across chains do not share a bin).
#' @param filter_IMGT_sequences Filter for productive IMGT CDR3s before clustering.
#' @param project Optional label column added to the result.
#' @param quiet Return invisibly.
#' @param ... Passed to `read_repertoire()`.
#'
#' @return Data frame with one or more `clonality_*` columns (or a single
#'   `clonality` column in the legacy single-pass case).
#'
#' @examples
#' results <- clonality(tra, vgene_col = "V_GENE_and_allele", jgene_col = "J_GENE_and_allele")
#' results <- clonality(IMGT_TCR_Summary_File, cdr3_col = "DNA_Junction")
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
                      extra_bin_cols = NULL,
                      ...) {

  # ── internal helpers ────────────────────────────────────────────────────
  recycle_to <- function(x, n, what) {
    if (length(x) == n) return(x)
    if (length(x) == 1) return(rep(x, n))
    stop(sprintf("'%s' must be length 1 or %d, got %d.", what, n, length(x)),
         call. = FALSE)
  }

  detect_tag <- function(vec_chr) {
    s <- paste(vec_chr, collapse = "")
    s <- gsub("[^A-Za-z]", "", s)
    if (nchar(s) == 0) return("nt")
    nt_frac <- nchar(gsub("[^ACGTNacgtn]", "", s)) / nchar(s)
    if (is.na(nt_frac)) return("nt")
    if (nt_frac >= 0.8) "nt" else "aa"
  }

  if (!is.numeric(mismatch) || any(mismatch < 0 | mismatch > 100)) {
    stop("'mismatch' must be between 0 and 100 (percent mismatch tolerance).",
         call. = FALSE)
  }

  # ── Read + validate ─────────────────────────────────────────────────────
  imported_df <- read_repertoire(data, ...)

  required <- c(ident_col, vgene_col, jgene_col, cdr3_col, extra_bin_cols)
  missing_cols <- setdiff(required, names(imported_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required column(s): %s",
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Drop rows missing ANY requested CDR3 column
  rm.cols  <- c(NA, "")
  keep_idx <- Reduce(`&`,
                     lapply(cdr3_col, function(cc) !(imported_df[[cc]] %in% rm.cols)))
  df <- imported_df[keep_idx, , drop = FALSE]

  if (isTRUE(filter_IMGT_sequences)) {
    df <- filter_sequences_imgt(df = df, cdr3_col = cdr3_col[[1]])
  }

  sub.df <- data.frame(
    cell_id = df[[ident_col]],
    v_genes = df[[vgene_col]],
    j_genes = df[[jgene_col]],
    stringsAsFactors = FALSE
  )

  if (any(!stats::complete.cases(sub.df[, c("v_genes", "j_genes")]))) {
    stop("Some rows contain missing V or J genes.", call. = FALSE)
  }

  carry_cols <- unique(c(cdr3_col, extra_bin_cols))
  simple.df <- simplify_genenames(
    data = cbind(sub.df, df[, carry_cols, drop = FALSE]),
    cell = cell
  )
  clonal.df <- simple.df[["clonal.df"]]

  # ── Determine pass count K = max(|cdr3_col|, |mismatch|) ────────────────
  # Each pass clusters a (cdr3_col, mismatch) pair. Length-1 args recycle.
  K <- max(length(cdr3_col), length(mismatch))
  cdr3_col          <- recycle_to(cdr3_col,          K, "cdr3_col")
  mismatch          <- recycle_to(mismatch,          K, "mismatch")
  stringdist_method <- recycle_to(stringdist_method, K, "stringdist_method")
  suffix            <- recycle_to(
    if (is.null(suffix)) NA_character_ else suffix,
    K, "suffix"
  )

  # ── Decide tags for output column names ─────────────────────────────────
  # Strategy:
  #   1. Auto-detect "nt"/"aa" from CDR3 content
  #   2. If duplicate tags result (e.g. same CDR3 col scanned at multiple
  #      thresholds), fall back to "mm<mismatch>"
  #   3. User-supplied cdr3_tags overrides everything
  if (is.null(cdr3_tags)) {
    tags <- vapply(cdr3_col, function(cc) detect_tag(clonal.df[[cc]]),
                   character(1))
    if (anyDuplicated(tags)) {
      tags <- paste0("mm", mismatch)
    }
  } else {
    tags <- recycle_to(cdr3_tags, K, "cdr3_tags")
  }
  # Ensure tags are unique; if still not, append index
  if (anyDuplicated(tags)) {
    tags <- make.unique(tags, sep = "_")
  }

  # ── Per-pass clustering ─────────────────────────────────────────────────
  out_list <- vector("list", K)

  for (k in seq_len(K)) {
    ccol <- cdr3_col[[k]]
    tmp  <- clonal.df
    tmp$CDR3  <- clonal.df[[ccol]]
    tmp$CDR3L <- nchar(tmp$CDR3)

    # Bin key: V + J + total CDR3 length, plus any user-supplied extra columns
    # (e.g. per-chain length string from assigntenx() to prevent paired cells
    # with compensating cross-chain indels from sharing a bin).
    bin_parts <- list(tmp[["v_genes"]], tmp[["j_genes"]], tmp[["CDR3L"]])
    if (!is.null(extra_bin_cols)) {
      for (cc in extra_bin_cols) {
        bin_parts[[length(bin_parts) + 1L]] <- as.character(tmp[[cc]])
      }
    }
    clonal_id <- do.call(paste, c(bin_parts, list(sep = "_")))

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

    tmp$project <- project
    if (!is.na(suffix[[k]])) {
      tmp$clonality <- paste(suffix[[k]], tmp$clonality, sep = "_")
    }

    out_list[[k]] <- tmp[, c("cell_id", "clonality"), drop = FALSE]
    # Single-pass case keeps the plain "clonality" name for back-compat.
    if (K > 1) {
      names(out_list[[k]])[2] <- paste0("clonality_", tags[[k]])
    }
  }

  # ── Merge passes + base context columns ─────────────────────────────────
  clonal.merge <- Reduce(function(x, y) dplyr::left_join(x, y, by = "cell_id"),
                         out_list)

  base_cols <- unique(c("cell_id", "v_genes", "j_genes", cdr3_col))
  clonal.merge <- dplyr::left_join(
    clonal.merge,
    clonal.df[, base_cols, drop = FALSE],
    by = "cell_id"
  )
  clonal.merge$project <- project

  if (isTRUE(output_original)) {
    clonal.merge <- safe_full_join_coalesce(
      x = imported_df,
      y = clonal.merge,
      by = stats::setNames("cell_id", ident_col)
    )
  } else {
    primary <- names(clonal.merge)[grepl("^clonality", names(clonal.merge))][1]
    clonal.merge <- clonal.merge[gtools::mixedorder(clonal.merge[[primary]]), ]
  }

  if (!is.null(file)) {
    openxlsx::write.xlsx(clonal.merge, file, rowNames = FALSE)
  }

  if (quiet) return(invisible(clonal.merge))
  clonal.merge
}
