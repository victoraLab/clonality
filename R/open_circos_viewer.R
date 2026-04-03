#' Open the bundled Circos D3 viewer
#'
#' Opens the package's standalone HTML viewer in the user's default browser.
#' The viewer is installed from `inst/circos/index.html`.
#'
#' @param open `Logical`. If `TRUE`, open the viewer in the default browser.
#'   If `FALSE`, only return the installed file path.
#'
#' @return Invisibly returns the path to the installed HTML file.
#'
#' @examples
#' \dontrun{
#' open_circos_viewer()
#' }
#' @export
open_circos_viewer <- function(open = TRUE) {
  viewer_path <- system.file("circos", "index.html", package = "clonality")

  if (!nzchar(viewer_path)) {
    stop("The Circos viewer could not be found in the installed clonality package.", call. = FALSE)
  }

  viewer_path <- normalizePath(viewer_path, winslash = "/", mustWork = TRUE)

  if (isTRUE(open)) {
    utils::browseURL(viewer_path)
  }

  invisible(viewer_path)
}
