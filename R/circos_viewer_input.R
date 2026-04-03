#' Convert a wide clonotype table into Circos viewer input
#'
#' Turns a wide table with one clonotype column and one column per group into
#' the long format expected by the bundled D3 Circos viewer.
#'
#' @param data A data frame in wide format.
#' @param clone_col Column name or position containing clonotype identifiers.
#' @param drop_zero `Logical`. If `TRUE`, rows with `size <= 0` are removed.
#'
#' @return A data frame with columns `clonotype`, `size`, and `group`.
#'   Group order follows the original column order in `data`.
#'
#' @examples
#' wide_df <- data.frame(
#'   Clone_ID = c("1.1", "2.1", "6.1"),
#'   M1LP_PC = c(2, 3, 3),
#'   M1PP_GC = c(11, 0, 0),
#'   M2LP_PC = c(0, 0, 1)
#' )
#'
#' circos_viewer_wide_to_long(wide_df)
#' @export
circos_viewer_wide_to_long <- function(data, clone_col = "Clone_ID", drop_zero = TRUE) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  clone_col <- resolve_circos_col(data, clone_col)
  value_cols <- setdiff(names(data), clone_col)
  group_order <- value_cols

  if (!length(value_cols)) {
    stop("`data` must contain at least one group column besides `clone_col`.", call. = FALSE)
  }

  out <- tidyr::pivot_longer(
    data = data,
    cols = dplyr::all_of(value_cols),
    names_to = "group",
    values_to = "size"
  )

  names(out)[names(out) == clone_col] <- "clonotype"
  out$clonotype <- as.character(out$clonotype)
  out$group <- as.character(out$group)
  out$size <- suppressWarnings(as.numeric(out$size))

  if (anyNA(out$clonotype) || anyNA(out$group) || anyNA(out$size)) {
    stop("Converted `clonotype`, `group`, and `size` values cannot contain missing data.", call. = FALSE)
  }

  if (isTRUE(drop_zero)) {
    out <- out[out$size > 0, , drop = FALSE]
  }

  out <- out[order(match(out$group, group_order), -out$size, out$clonotype), c("clonotype", "size", "group")]
  rownames(out) <- NULL
  out
}

#' Summarize clonotypes into Circos viewer input
#'
#' Converts cell-level or contig-level tables into the long format expected by
#' the bundled D3 Circos viewer. If `size_col` is not supplied, the function
#' counts rows per clonotype and group.
#'
#' @param data A data frame containing clonotype and group information.
#' @param clonotype_col Column name or position with clonotype identifiers.
#' @param group_col Column name or position with the group or condition.
#' @param size_col Optional column name or position with numeric abundances.
#'   If omitted, rows are counted.
#' @param drop_zero `Logical`. If `TRUE`, rows with `size <= 0` are removed.
#'
#' @return A data frame with columns `clonotype`, `size`, and `group`.
#'   If `group_col` is a factor, ideograms follow its level order.
#'
#' @examples
#' example_df <- data.frame(
#'   clonality = c("1.1", "1.1", "2.1", "6.1", "6.1"),
#'   condition = c("M1LP_PC", "M1PP_GC", "M1LP_PC", "M1LP_PC", "M2LP_PC")
#' )
#'
#' circos_viewer_input(example_df, group_col = "condition")
#' @export
circos_viewer_input <- function(data,
                                clonotype_col = "clonality",
                                group_col,
                                size_col = NULL,
                                drop_zero = TRUE) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  clonotype_col <- resolve_circos_col(data, clonotype_col)
  group_col <- resolve_circos_col(data, group_col)
  group_values <- data[[group_col]]
  group_order <- if (is.factor(group_values)) levels(group_values) else unique(as.character(group_values))

  if (is.null(size_col)) {
    out <- data |>
      dplyr::group_by(.data[[clonotype_col]], .data[[group_col]]) |>
      dplyr::summarise(size = dplyr::n(), .groups = "drop")

    names(out) <- c("clonotype", "group", "size")
  } else {
    size_col <- resolve_circos_col(data, size_col)

    out <- data |>
      dplyr::group_by(.data[[clonotype_col]], .data[[group_col]]) |>
      dplyr::summarise(size = sum(as.numeric(.data[[size_col]]), na.rm = TRUE), .groups = "drop")

    names(out) <- c("clonotype", "group", "size")
  }

  out$clonotype <- as.character(out$clonotype)
  out$group <- as.character(out$group)
  out$size <- as.numeric(out$size)

  if (anyNA(out$clonotype) || anyNA(out$group) || anyNA(out$size)) {
    stop("Summarized `clonotype`, `group`, and `size` values cannot contain missing data.", call. = FALSE)
  }

  if (isTRUE(drop_zero)) {
    out <- out[out$size > 0, , drop = FALSE]
  }

  out <- out[order(match(out$group, group_order), -out$size, out$clonotype), c("clonotype", "size", "group")]
  rownames(out) <- NULL
  out
}

resolve_circos_col <- function(data, col) {
  if (is.numeric(col)) {
    if (length(col) != 1L || is.na(col) || col < 1L || col > ncol(data)) {
      stop("Numeric column indices must point to an existing column.", call. = FALSE)
    }
    return(names(data)[[col]])
  }

  col <- as.character(col)
  if (!col %in% names(data)) {
    stop(sprintf("Column `%s` was not found in `data`.", col), call. = FALSE)
  }

  col
}
