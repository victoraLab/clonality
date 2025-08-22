#' Evaluate Plate Quality with Negative Binomial Fit and Efficiency Metrics
#'
#' Fits a Negative Binomial distribution to real contig depths per plate and calculates
#' efficiency metrics such as the number of wells detected, collapsed contigs, and success rate.
#'
#' @param analysis_output A list returned by `analyze_well_contigs()` with elements `collapsed` and `real_contigs_input`.
#' @param plate_size Expected number of wells per plate (default: 96).
#' @return A tibble with one row per plate and the following columns:
#' \describe{
#'   \item{PLATE}{Plate ID}
#'   \item{mu}{Mean of fitted NB model}
#'   \item{size}{Dispersion parameter of NB model}
#'   \item{n_wells_detected}{Number of unique wells with successful contigs}
#'   \item{n_successful_contigs}{Total number of real contigs}
#'   \item{n_fully_collapsed}{Number of contigs marked as collapsed}
#'   \item{efficiency_score}{Percent of expected wells (default 96) that were detected}
#'   \item{collapse_quality}{Proportion of collapsed contigs (consistency measure)}
#' }
#' @export
evaluate_plate_quality <- function(analysis_output, plate_size = 96) {
  if (!all(c("collapsed", "real_contigs_input") %in% names(analysis_output))) {
    stop("Input must be a list returned by `analyze_well_contigs()`.")
  }

  df_all <- analysis_output$collapsed
  df_real <- dplyr::filter(df_all, background_tag == "REAL")


  # Aggregate real contig depths per PLATE_WELL for NB fitting
  per_well_depth <- df_real %>%
    dplyr::group_by(PLATE, PLATE_WELL) %>%
    dplyr::summarise(CONTIG_DEPTH = sum(CONTIG_DEPTH, na.rm = TRUE), .groups = "drop")

  # Fit NB model per plate
  nb_stats <- per_well_depth %>%
    group_by(PLATE) %>%
    group_modify(~ {
      data <- .x
      n_missing <- plate_size - nrow(data)
      if (n_missing > 0) {
        missing_rows <- tibble(
          PLATE = unique(data$PLATE),
          PLATE_WELL = paste0("MISSING_", seq_len(n_missing)),
          CONTIG_DEPTH = 0
        )
        data <- bind_rows(data, missing_rows)
      }

      fit <- tryCatch(
        MASS::fitdistr(data$CONTIG_DEPTH, "Negative Binomial"),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        return(tibble(mu = NA_real_, size = NA_real_))
      }

      tibble(mu = fit$estimate["mu"], size = fit$estimate["size"])
    }) %>%
    ungroup()

# Efficiency metrics using full collapsed data filtered to REAL
  efficiency <- df_all %>%
    dplyr::filter(background_tag == "REAL") %>%
    group_by(PLATE) %>%
    summarise(
      n_wells_detected = n_distinct(WELL),
      n_successful_contigs = n(),
      n_fully_collapsed = sum(collapsed, na.rm = TRUE),
      efficiency_score = round((n_wells_detected / plate_size) * 100, 2),
      collapse_quality = round(n_fully_collapsed / n_successful_contigs, 2),
      .groups = "drop"
    )


  # Merge results
  merged_summary <- left_join(nb_stats, efficiency, by = "PLATE")
  merged_summary$score <- merged_summary$mu * merged_summary$size
  return(merged_summary)
}
