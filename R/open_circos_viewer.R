#' Open the bundled Circos D3 viewer
#'
#' Opens the package's standalone HTML viewer in the user's default browser.
#' The viewer is installed from `inst/circos/index.html`. If `data` is
#' supplied, the viewer is preloaded and rendered automatically from a
#' clonality-style data frame without writing a manual CSV file.
#'
#' @param data Optional data frame to preload into the viewer.
#' @param ideogram_by Optional column name or position defining the ideograms.
#'   When `data` is provided, this argument is required.
#' @param clonality_col Column name or position containing clonotype IDs.
#' @param size_col Optional column name or position with per-row abundances.
#'   If omitted, rows are counted within each `clonality_col` x `ideogram_by`
#'   combination.
#' @param band_labels `Logical`. If `TRUE`, show clone IDs as band labels when
#'   preloading data into the viewer.
#' @param open `Logical`. If `TRUE`, open the viewer in the default browser.
#'   If `FALSE`, only return the installed file path.
#'
#' @return Invisibly returns the path to the installed HTML file.
#'
#' @examples
#' \dontrun{
#' open_circos_viewer()
#'
#' example_df <- data.frame(
#'   clonality = c("1.1", "1.1", "2.1", "6.1", "6.1"),
#'   tissue = factor(
#'     c("M1LP_PC", "M1PP_GC", "M1LP_PC", "M1LP_PC", "M2LP_PC"),
#'     levels = c("M1LP_PC", "M1PP_GC", "M2LP_PC")
#'   )
#' )
#'
#' open_circos_viewer(example_df, ideogram_by = "tissue", open = FALSE)
#' }
#' @export
open_circos_viewer <- function(data = NULL,
                               ideogram_by = NULL,
                               clonality_col = "clonality",
                               size_col = NULL,
                               band_labels = TRUE,
                               open = TRUE) {
  viewer_path <- system.file("circos", "index.html", package = "clonality")

  if (!nzchar(viewer_path)) {
    stop("The Circos viewer could not be found in the installed clonality package.", call. = FALSE)
  }

  viewer_path <- normalizePath(viewer_path, winslash = "/", mustWork = TRUE)

  if (!is.null(data)) {
    if (missing(ideogram_by) || is.null(ideogram_by)) {
      stop("When `data` is provided, `ideogram_by` must also be provided.", call. = FALSE)
    }

    viewer_data <- circos_viewer_input(
      data = data,
      clonotype_col = clonality_col,
      group_col = ideogram_by,
      size_col = size_col
    )

    viewer_path <- build_circos_viewer_html(
      template_path = viewer_path,
      viewer_data = viewer_data,
      band_labels = band_labels
    )
  }

  if (isTRUE(open)) {
    utils::browseURL(viewer_path)
  }

  invisible(viewer_path)
}

build_circos_viewer_html <- function(template_path, viewer_data, band_labels = TRUE) {
  payload <- paste(
    "<script>",
    sprintf(
      "window.__CLONALITY_CIRCOS_PRELOAD__ = { mode: \"clonotypes\", data: \"%s\", autoRender: true, showBandLabels: %s };",
      encode_js_string(circos_viewer_to_tsv(viewer_data)),
      tolower(as.character(isTRUE(band_labels)))
    ),
    "</script>",
    sep = "\n"
  )

  template_lines <- readLines(template_path, warn = FALSE)
  marker <- '<script src="https://d3js.org/d3.v7.min.js"></script>'
  match_idx <- match(marker, template_lines)

  if (is.na(match_idx)) {
    stop("Could not inject preloaded data into the Circos viewer template.", call. = FALSE)
  }

  template_lines <- append(template_lines, payload, after = match_idx - 1L)

  temp_path <- tempfile(pattern = "clonality-circos-", fileext = ".html")
  writeLines(template_lines, temp_path, useBytes = TRUE)
  normalizePath(temp_path, winslash = "/", mustWork = TRUE)
}

circos_viewer_to_tsv <- function(viewer_data) {
  header <- paste(names(viewer_data), collapse = "\t")
  rows <- apply(viewer_data, 1, function(row) paste(as.character(row), collapse = "\t"))
  paste(c(header, rows), collapse = "\n")
}

encode_js_string <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"", "\\\\\"", x)
  x <- gsub("\r", "", x, fixed = TRUE)
  x <- gsub("\n", "\\\\n", x, fixed = TRUE)
  gsub("</script>", "<\\\\/script>", x, fixed = TRUE)
}
