#' Extract Plate Metadata from Cell IDs
#'
#' Extracts plate, well, contig, and read number from `cell_id` strings typically produced by the `clonality()` function.
#' Supports three barcode formats: where WELL precedes PLATE (e.g. "A01P09"), PLATE precedes WELL (e.g. "P09A01"), or auto-detected.
#'
#' @param df A data frame containing a `cell_id` column with formatted strings.
#' @param cell_id_col Name of the column containing cell IDs. Default is `"cell_id"`.
#' @param barcode_format One of `"WELL_PLATE"`, `"PLATE_WELL"`, or `"auto"` (default: `"WELL_PLATE"`).
#' @return A data frame with columns: `ProjectID`, `PLATE_WELL`, `PLATE`, `WELL`, `CONTIG_NUMBER`, `CONTIG_DEPTH`, and `ggplate_well`.
#' @export
extract_plate_metadata <- function(df, cell_id_col = "cell_id", barcode_format = "WELL_PLATE") {
  ids <- df[[cell_id_col]]

  # Auto-detect barcode format
  if (barcode_format == "auto") {
    match_wp <- stringr::str_match(ids, "^(.*?)([A-H]\\d{2}P\\d{2})_(\\d+)-(\\d+)$")
    match_pw <- stringr::str_match(ids, "^(.*?)(P\\d{2}[A-H]\\d{2})_(\\d+)-(\\d+)$")

    n_wp <- sum(!is.na(match_wp[, 1]))
    n_pw <- sum(!is.na(match_pw[, 1]))

    if (n_wp >= n_pw) {
      matches <- match_wp
      format <- "WELL_PLATE"
    } else {
      matches <- match_pw
      format <- "PLATE_WELL"
    }
  } else {
    format <- match.arg(barcode_format, choices = c("WELL_PLATE", "PLATE_WELL"))
    regex <- switch(format,
                    WELL_PLATE = "^(.*?)([A-H]\\d{2}P\\d{2})_(\\d+)-(\\d+)$",
                    PLATE_WELL = "^(.*?)(P\\d{2}[A-H]\\d{2})_(\\d+)-(\\d+)$"
    )
    matches <- stringr::str_match(ids, regex)
  }

  if (anyNA(matches)) stop("Failed to parse some cell IDs. Check barcode_format or cell ID formatting.")

  PLATE_WELL <- matches[, 3]
  plate <- stringr::str_extract(PLATE_WELL, "P\\d{2}")
  well  <- stringr::str_extract(PLATE_WELL, "[A-H]\\d{2}")

  # Format well for ggplate
  well_fmt <- sub("^([A-H])0?([1-9]|1[0-9]|2[0-4])$", "\\1\\2", well)

  result <- data.frame(
    ProjectID      = matches[, 2],
    PLATE_WELL     = PLATE_WELL,
    PLATE          = plate,
    WELL           = well,
    CONTIG_NUMBER  = as.integer(matches[, 4]),
    CONTIG_DEPTH   = as.integer(matches[, 5]),
    ggplate_well   = well_fmt,
    stringsAsFactors = FALSE
  )

  return(result)
}
#' Summarize and Plot Contig Depth per Well for a Specific Plate
#'
#' Aggregates contig depth per well within a specified plate and plots it using `ggplate`.
#'
#' @param df A data frame from `plot_plate_metadata()`, containing `PLATE`, `ggplate_well`, and `CONTIG_DEPTH`.
#' @param plate_id A character string specifying which plate to plot (e.g., `"P01"`).
#' @param stat Summary statistic to apply for wells with multiple contigs. One of `"mean"`, `"median"`, or a custom function. Default is `"mean"`.
#' @param plate_size Number of wells in the plate (e.g., 96, 384). Default is 96.
#' @param plate_type Well type: `"round"` or `"square"`. Default is `"round"`.
#' @param label Logical. Whether to label each well with its value. Default is `TRUE`.
#' @import ggplate
#' @return A `ggplot` object for the requested plate.
#' @export
plot_plate_contig_depth <- function(df,
                                    plate_id,
                                    stat = "mean",
                                    plate_size = 96,
                                    plate_type = "round",
                                    label = TRUE) {
  if (!requireNamespace("ggplate", quietly = TRUE)) {
    stop("Package 'ggplate' is required. Install it with install.packages('ggplate')")
  }

  if (!plate_id %in% df$PLATE) {
    stop(sprintf("Plate '%s' not found in data.", plate_id))
  }

  stat_fun <- switch(stat,
                     mean = mean,
                     median = median,
                     stop("Unsupported 'stat' value. Use 'mean', 'median', or provide a function.")
  )

  require(ggplate)

  df_plate <- dplyr::filter(df, PLATE == plate_id)

  summary_df <- df_plate %>%
    dplyr::group_by(ggplate_well) %>%
    dplyr::summarise(Value = round(stat_fun(CONTIG_DEPTH, na.rm = TRUE), 2), .groups = "drop")

  plate_plot(
    data = summary_df,
    position = ggplate_well,
    value = Value,
    title = sprintf("%s", plate_id),
    label = if (label) Value else NULL,
    plate_size = plate_size,
    plate_type = plate_type
  )
}
#' Plot Multiple Plates Side-by-Side Using ggplate and patchwork
#'
#' Plots contig depth per well for multiple plates using `plot_plate_contig_depth()` and arranges them with `patchwork`.
#'
#' @param df A data frame from `extract_plate_metadata()`.
#' @param plate_ids A character vector of plate IDs to include (e.g., c("P01", "P02")).
#' @param stat Summary statistic to apply. One of `"mean"`, `"median"`, or a custom function. Default is `"mean"`.
#' @param plate_size Number of wells per plate. Default is 96.
#' @param plate_type Well type. `"round"` or `"square"`. Default is `"round"`.
#' @param label Logical. Show values in wells. Default is `TRUE`.
#' @param ncol Number of columns in patchwork layout. Default is `2`.
#'
#' @return A combined patchwork plot object.
#' @export
plot_multiple_plates <- function(df,
                                 plate_ids,
                                 stat = "mean",
                                 plate_size = 96,
                                 plate_type = "round",
                                 label = TRUE,
                                 ncol = 2) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install it with install.packages('patchwork')")
  }

  plots <- lapply(plate_ids, function(pid) {
    p <- plot_plate_contig_depth(df, plate_id = pid, stat = stat, plate_size = plate_size, plate_type = plate_type, label = label)
    p + ggplot2::ggtitle(sprintf("Plate %s", pid))
  })

  patchwork::wrap_plots(plots, ncol = ncol)
}
