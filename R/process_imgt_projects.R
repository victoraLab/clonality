#' Process Multiple IMGT Datasets in Parallel
#'
#' Processes multiple IMGT summary folders using parallel execution.
#' Tags each dataset with its project name and returns merged results.
#'
#' @param project_paths Vector of IMGT project folder paths.
#' @param cell_type "B" or "T" for B/T cell analysis.
#' @param mismatch Numeric value (0–100) for CDR3 mismatch tolerance.
#' @param exclude_index Optional vector of indices to exclude.
#' @param workers Number of parallel workers (cores). Default = available cores - 1.
#'
#' @return A list with combined `metadata`, `quality`, and `real_contigs`.
#' @export
#' @import dplyr
#' @importFrom furrr future_map
#' @importFrom stringr str_extract
process_imgt_projects <- function(project_paths,
                                  cell_type = "B",
                                  mismatch = 0,
                                  exclude_index = NULL,
                                  workers = parallel::detectCores() - 1) {
  if (!is.null(exclude_index)) {
    project_paths <- project_paths[-exclude_index]
  }

  # Load parallel backend
  future::plan(future::multisession, workers = workers)

  results <- furrr::future_map(project_paths, function(project_path) {
    project_name <- basename(project_path)
    imgt_data <- parse_imgt(project_path)

    clonality(data = imgt_data,
              cdr3_col = "DNA_Junction",
              cell = cell_type,
              mismatch = mismatch,
              project = project_name)

    analysed <- analyze_well_contigs(output, barcode_format = "auto")
    metadata <- extract_plate_metadata(analysed$real_contigs_input, barcode_format = "auto")
    quality  <- evaluate_plate_quality(analysed)

    metadata$project <- quality$project <- analysed$real_contigs_input$project <- project_name

    list(meta = metadata, qual = quality, real = analysed$real_contigs_input)
  })

  # Collapse list into data.frames
  combined_metadata     <- dplyr::bind_rows(purrr::map(results, "meta"))
  combined_quality      <- dplyr::bind_rows(purrr::map(results, "qual"))
  combined_real_contigs <- dplyr::bind_rows(purrr::map(results, "real"))

  combined_quality$Author <- stringr::str_extract(combined_quality$project, "^..")

  return(list(
    metadata = combined_metadata,
    quality = combined_quality,
    real_contigs = combined_real_contigs
  ))
}
