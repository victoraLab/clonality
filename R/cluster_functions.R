#' Cluster Clonal IDs Based on CDR3 Similarity
#'
#' Clusters sequences with identical V-J-CDR3L identifiers into subclones based on CDR3 sequence similarity.
#' Uses hierarchical clustering with a specified height threshold. Clone groups are sorted by decreasing size.
#'
#' @inheritParams clonality
#' @param clonal.df A simplified data frame with at least CDR3, `cell_id`, and gene columns.
#' @param clonal_id A character vector of V-J-CDR3L combinations.
#' @param cdr3_lengths A numeric vector of CDR3 lengths, same length as `clonal_id`.
#' @param h_height The height threshold (between 0 and 1) for clustering.
#'
#' @return A data frame with a `clonality` column assigning clone IDs.
#'
#' @export
cluster_duplicates_by_cdr3 <- function(clonal.df, clonal_id, cdr3_lengths, dist_method, h_height) {
  duplicated_ids <- unique(clonal_id[duplicated(clonal_id)])

  clonal.df$clonality <- NA
  groups <- vector("list", length = length(duplicated_ids))
  sizes <- numeric(length(groups))

  # Step 1: build group info (distance matrices and indices)
  for (i in seq_along(duplicated_ids)) {
    id <- duplicated_ids[i]
    idx <- which(clonal_id == id)
    dist_mat <- stringdist::stringdistmatrix(clonal.df[idx, "CDR3"], method = dist_method, useNames = FALSE)
    dist_mat <- dist_mat / cdr3_lengths[idx][1]
    groups[[i]] <- list(idx = idx, dist = dist_mat)
    sizes[i] <- length(idx)
  }

  # Step 2: sort groups by descending size
  ord <- order(sizes, decreasing = TRUE)
  groups <- groups[ord]

  # Step 3: assign clonality
  clone_counter <- 1
  for (g in groups) {
    hclust_result <- hclust(g$dist)
    cutree_result <- cutree(hclust_result, h = h_height)

    for (k in unique(cutree_result)) {
      group_idx <- g$idx[cutree_result == k]
      clonal.df$clonality[group_idx] <- sprintf("%s.%s", clone_counter, k)
    }

    clone_counter <- clone_counter + 1
  }

  # Step 4: assign unique IDs to leftovers
  unique_bool <- is.na(clonal.df$clonality)
  if (any(unique_bool)) {
    clonal.df$clonality[unique_bool] <- sprintf("U%s", seq_len(sum(unique_bool)))
  }

  return(clonal.df)
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
