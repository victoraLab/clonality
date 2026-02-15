safe_full_join_coalesce <- function(x, y, by, relationship = "one-to-one") {
  stopifnot(is.data.frame(x), is.data.frame(y))
  stopifnot(is.character(by), length(by) >= 1)

  out <- dplyr::full_join(
    x, y,
    by = by,
    suffix = c(".x", ".y"),
    relationship = relationship
  )

  # columns that exist from both sides (after suffixing)
  x_cols <- sub("\\.x$", "", names(out)[grepl("\\.x$", names(out))])
  y_cols <- sub("\\.y$", "", names(out)[grepl("\\.y$", names(out))])
  dup <- intersect(x_cols, y_cols)

  # coalesce duplicates: prefer x, fill from y
  for (nm in dup) {
    out[[nm]] <- dplyr::coalesce(out[[paste0(nm, ".x")]],
                                 out[[paste0(nm, ".y")]])
  }

  # drop the suffixed columns
  out <- out[, !grepl("\\.(x|y)$", names(out)), drop = FALSE]

  out
}
