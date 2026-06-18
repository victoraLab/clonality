#' Cluster duplicated clonal IDs by CDR3 similarity
#'
#' Clusters sequences sharing the same V–J–CDR3 length identifier into subclones
#' based on CDR3 sequence similarity. Clone bins are generated based on same V
#' J and CDR3 length.
#' When `h_height = 0`, subclones are defined by exact CDR3 identity.
#'
#' Clone bins are processed in decreasing order of size and assigned hierarchical
#' labels of the form `<bin>.<subcluster>`. Rows not belonging to duplicated bins
#' receive unique `U*` identifiers.
#'
#' @param clonal.df A simplified data frame containing at least a `CDR3` column.
#' @param clonal_id A character vector of V–J–CDR3 length identifiers.
#' @param cdr3_lengths A numeric vector of CDR3 lengths, same length as `clonal_id`.
#' @param dist_method String distance method passed to `stringdist::stringdistmatrix()`.
#' @param h_height Clustering height threshold (0–1); typically `mismatch / 100`. Applied to the
#'   length-normalized distance (`stringdist / CDR3 length`). Ignored when `max_edits` is set.
#' @param max_edits Optional absolute edit-distance cutoff. When non-`NULL`, the tree is cut on the
#'   RAW (un-normalized) `stringdist` at this height instead of the length-normalized `h_height` —
#'   e.g. `max_edits = 1` with `dist_method = "hamming"` is a true absolute Hamming-1 (one mismatch
#'   regardless of CDR3 length). `max_edits = 0` is exact-identity (same as `h_height = 0`).
#' @param per_chain_cdr3_cols Optional character vector of column names holding the
#'   individual per-chain CDR3 sequences (e.g. `c("cdr3_nt_IGH", "cdr3_nt_IGK")`).
#'   When provided together with `max_edits`, the distance between two cells is
#'   defined as the **maximum Hamming distance across all chains independently** —
#'   i.e. two cells cluster together only when every chain's CDR3 is within
#'   `max_edits` edits of its counterpart. This is a stricter, biologically more
#'   accurate criterion than applying `max_edits` to the concatenated CDR3 string.
#'   Ignored when `max_edits` is `NULL`.
#'
#' @return A data frame with an added `clonality` column.
#'
#' @export
cluster_duplicates_by_cdr3 <- function(clonal.df, clonal_id, cdr3_lengths,
                                          dist_method = "osa", h_height = 0,
                                          max_edits = NULL,
                                          per_chain_cdr3_cols = NULL) {
  absolute  <- !is.null(max_edits)
  per_chain <- absolute && !is.null(per_chain_cdr3_cols) && length(per_chain_cdr3_cols) > 0L

  # Identify bin IDs that occur more than once (i.e., bins needing within-bin subclustering).
  duplicated_ids <- unique(clonal_id[duplicated(clonal_id)])

  # Initialize output column (explicit character to avoid factor surprises).
  clonal.df$clonality <- NA_character_

  # Process the largest bins first (mirrors original behavior and helps debugging/profiling).
  sizes <- vapply(duplicated_ids, function(id) sum(clonal_id == id), numeric(1))
  duplicated_ids <- duplicated_ids[order(sizes, decreasing = TRUE)]

  clone_counter <- 1L

  # Iterate over each duplicated bin (each corresponds to a fixed V/J/length combination).
  for (id in duplicated_ids) {
    idx  <- which(clonal_id == id)      # row indices belonging to this bin

    # ── Per-chain independent Hamming mode ────────────────────────────────────
    # Requires max_edits (absolute cutoff) + per_chain_cdr3_cols.
    # Distance between two cells = max(hamming(chain_k_a, chain_k_b) for each chain k).
    # Two cells cluster together iff every chain is within max_edits of its counterpart.
    if (per_chain) {
      seqs_list <- lapply(per_chain_cdr3_cols, function(col) {
        s <- as.character(clonal.df[[col]][idx])
        s[is.na(s)] <- ""
        s
      })

      # Fast path: exact identity across all chains
      if (max_edits == 0L) {
        combined <- do.call(paste, c(seqs_list, list(sep = "|")))
        grp <- match(combined, unique(combined))
        clonal.df$clonality[idx] <- sprintf("%d.%d", clone_counter, grp)
        clone_counter <- clone_counter + 1L
        next
      }

      # Single row in bin: trivially one cluster
      if (length(idx) == 1L) {
        clonal.df$clonality[idx] <- sprintf("%d.1", clone_counter)
        clone_counter <- clone_counter + 1L
        next
      }

      # For each chain: pairwise Hamming on unique sequences, then expand to all rows.
      # pmax across chains gives the per-pair "worst chain" distance.
      chain_dist_mats <- lapply(seqs_list, function(seqs) {
        uniq_seqs <- unique(seqs)
        d_uniq <- stringdist::stringdistmatrix(
          uniq_seqs, uniq_seqs, method = dist_method, useNames = FALSE
        )
        k <- match(seqs, uniq_seqs)
        d_uniq[k, k, drop = FALSE]
      })

      dist_obj <- stats::as.dist(Reduce(pmax, chain_dist_mats))

      # Complete linkage: within each cluster, every pair has worst-chain distance <= max_edits.
      hc <- stats::hclust(dist_obj, method = "complete")
      ct <- stats::cutree(hc, h = max_edits)

      clonal.df$clonality[idx] <- sprintf("%d.%d", clone_counter, ct)
      clone_counter <- clone_counter + 1L
      next
    }

    # ── Original single-CDR3 logic ────────────────────────────────────────────
    seqs <- clonal.df$CDR3[idx]         # original per-row CDR3 strings for this bin

    # Collapse to unique CDR3 sequences to avoid N^2 work on repeated identical reads.
    uniq_seqs <- unique(seqs)

    # FAST PATH: a zero threshold (percent h_height == 0, or absolute max_edits == 0) =>
    # exact identity grouping. Each unique CDR3 within this bin gets a unique subcluster label.
    if (if (absolute) max_edits == 0 else h_height == 0) {
      # For each row, map its CDR3 to the index of the corresponding unique CDR3.
      # Example: if 90 rows share the same CDR3, they all map to the same k.
      k <- match(seqs, uniq_seqs)  # integers 1..length(uniq_seqs)
      clonal.df$clonality[idx] <- sprintf("%d.%d", clone_counter, k)
      clone_counter <- clone_counter + 1L
      next
    }

    # If the bin contains only one unique CDR3, all rows must belong to the same subcluster.
    # No distance computation is needed.
    if (length(uniq_seqs) == 1L) {
      clonal.df$clonality[idx] <- sprintf("%d.%d", clone_counter, 1L)
      clone_counter <- clone_counter + 1L
      next
    }

    # mismatch > 0 and >1 unique CDR3:
    # Compute pairwise distances only across unique CDR3 sequences (size = U x U).
    dist_mat <- stringdist::stringdistmatrix(
      uniq_seqs, uniq_seqs,
      method = dist_method,
      useNames = FALSE
    )

    # Absolute mode: cut on the raw stringdist (e.g. Hamming-1 = one mismatch, any length).
    # Percent mode (default): normalize by CDR3 length (constant within this bin) and cut at
    # h_height = mismatch/100.
    if (absolute) {
      dist_obj <- stats::as.dist(dist_mat)
      cut_h <- max_edits
    } else {
      L <- cdr3_lengths[idx][1]
      dist_obj <- stats::as.dist(dist_mat / L)
      cut_h <- h_height
    }

    # Hierarchical clustering + cut at the requested height.
    hc <- stats::hclust(dist_obj)
    ct <- stats::cutree(hc, h = cut_h)  # cluster label for each entry of uniq_seqs

    # Map per-unique cluster labels back to every original row.
    m <- match(seqs, uniq_seqs)
    clonal.df$clonality[idx] <- sprintf("%d.%d", clone_counter, ct[m])

    clone_counter <- clone_counter + 1L
  }

  # Rows that never belonged to a duplicated bin get unique IDs (same behavior as original).
  unique_bool <- is.na(clonal.df$clonality)
  if (any(unique_bool)) {
    clonal.df$clonality[unique_bool] <- paste0("U", seq_len(sum(unique_bool)))
  }

  clonal.df
}
#' Assign Unique Clone Identifiers
#'
#' Assigns a unique ID to each row (e.g., `"U1"`, `"U2"`, ...) when no duplicates are found.
#'
#' @param df A data frame whose rows will receive unique clonality labels.
#'
#' @return A character vector of unique clone IDs.
#'
#' @export
assign_unique_ids <- function(df) {
  sprintf("U%s", seq_len(nrow(df)))
}
