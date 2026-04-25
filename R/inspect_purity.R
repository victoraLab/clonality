#' Summarize sticky-ends rescue purity per chain
#'
#' Prints a per-chain summary of single-cell rescue results after a sticky-ends pass
#' (see [add_seurat_metadata()]). For each chain it reports the number of singles,
#' how many were matched against the paired pool, the median dominant-purity score,
#' and how many were rescued at typical purity thresholds.
#'
#' Use this to decide whether the rescue is biologically trustworthy on your dataset:
#' chains with a median purity well below 1.0 are dominated by ambiguous matches
#' (e.g. light-chain CDR3s that recur across many unrelated paired clones). The
#' default chain-aware purity policy in [add_seurat_metadata()] already skips non-D
#' chains, so this helper is mostly diagnostic.
#'
#' @param x A Seurat object with `Res_*` columns added by [add_seurat_metadata()],
#'   or a data frame containing those columns.
#' @param chain_col Character. Column whose values identify the source single-chain
#'   pool. Defaults to `"Res_mapping_source"` for plain data frames, but the function
#'   reconstructs the chain from the rescue annotation when present.
#' @param purity_col Character. Column with the per-cell dominant-purity score.
#'   Default `"Res_dominant_purity"`.
#' @param correction_col Character. Column with the rescue outcome label.
#'   Default `"Res_clonality.correction"`.
#' @param thresholds Numeric. Purity thresholds at which to count rescues.
#'   Default `c(0.8, 0.9, 1.0)`.
#'
#' @return A `data.frame` with one row per chain, returned invisibly. Also printed
#'   as a tidy table.
#'
#' @seealso [add_seurat_metadata()]
#'
#' @export
inspect_purity <- function(x,
                           chain_col      = "Res_mapping_source",
                           purity_col     = "Res_dominant_purity",
                           correction_col = "Res_clonality.correction",
                           thresholds     = c(0.8, 0.9, 1.0)) {

  meta <- if (inherits(x, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat is required to inspect purity directly from a Seurat object.",
           call. = FALSE)
    }
    x@meta.data
  } else if (is.data.frame(x)) {
    x
  } else {
    stop("`x` must be a Seurat object or a data.frame.", call. = FALSE)
  }

  if (!purity_col %in% names(meta) || !correction_col %in% names(meta)) {
    stop("Required columns not found. Run add_seurat_metadata(sticky = TRUE) first.",
         call. = FALSE)
  }

  # Reconstruct chain identity from the embedded source where possible.
  # add_seurat_metadata() preserves the .clonal_source-derived suffixes implicitly
  # via the "Res_clonal_source" or "Res_mapping_source" columns; if neither is
  # available, fall back to the chain inferred from a CDR3 column name pattern.
  chain <- if ("Res_clonal_source" %in% names(meta)) {
    sub("^Clonal\\.output\\.10x", "", as.character(meta[["Res_clonal_source"]]))
  } else if ("Res_chain" %in% names(meta)) {
    as.character(meta[["Res_chain"]])
  } else if (chain_col %in% names(meta)) {
    as.character(meta[[chain_col]])
  } else {
    rep(NA_character_, nrow(meta))
  }

  # Restrict to single-chain rows (no underscore = single chain class)
  is_single <- !is.na(chain) & !grepl("_", chain)
  if (!any(is_single)) {
    message("inspect_purity: no single-chain rows detected -- nothing to report.")
    return(invisible(data.frame()))
  }

  meta <- meta[is_single, , drop = FALSE]
  chain <- chain[is_single]

  groups <- split(meta, chain)

  rows <- lapply(names(groups), function(ch) {
    g <- groups[[ch]]
    purity_vals <- as.numeric(g[[purity_col]])
    matched     <- !is.na(purity_vals)
    n           <- nrow(g)
    n_matched   <- sum(matched)
    median_p    <- if (n_matched > 0) stats::median(purity_vals[matched]) else NA_real_
    n_rescued   <- sum(g[[correction_col]] %in% c("Normal", "Corrected"))
    rescued_at  <- vapply(thresholds, function(t) {
      sum(matched & purity_vals >= t)
    }, integer(1))

    out <- data.frame(
      chain         = ch,
      singles       = n,
      matched       = n_matched,
      median_purity = round(median_p, 3),
      rescued       = n_rescued,
      stringsAsFactors = FALSE
    )
    for (k in seq_along(thresholds)) {
      out[[paste0("at_", format(thresholds[k]))]] <- rescued_at[k]
    }
    out
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL

  cat("Sticky-ends rescue summary\n")
  print(result, row.names = FALSE)

  # Heuristic warning when blanket rescue would be dangerous
  low_purity_chains <- result$chain[!is.na(result$median_purity) & result$median_purity < 0.5]
  if (length(low_purity_chains) > 0) {
    message(
      "  -> ", paste(low_purity_chains, collapse = ", "),
      " median purity below 0.5: blanket rescue of these chains is not recommended.\n",
      "     The default chain-aware policy in add_seurat_metadata() already skips non-D chains."
    )
  }

  invisible(result)
}
