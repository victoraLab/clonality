#' Annotate 10x clones with clonality information
#'
#' Loads and processes 10x Genomics V(D)J output (either as a data frame, file path, or object name), filters and groups by cell barcode, and assigns clonal definitions.
#'
#' Supports BCR, TCR, and γδ TCR data from Cell Ranger versions 3 (18 columns) and 7 (31 columns).
#'
#' @param data Data frame, character string (object name), or full path to a `filtered_contig_annotations.csv` or `all_contig_annotations.csv` file.
#' @param method Character. Clonality inference method. One of: `"unique_paired"`, `"sticky_ends"`, `"unique_all"`. Default is `"unique_paired"`.
#' @param only_productive Logical. Whether to filter non-productive chains. Default is `TRUE`.
#' @param only_true_cells Logical. Whether to filter out low-quality/non-true cells (i.e., `is_cell != TRUE`). Default is `TRUE`.
#' @param clonality_input Named vector. Parameters passed to the internal `clonality()` function (e.g., mismatch thresholds).
#' @param cell Character. Cell type to process. One of: `"B"` (B cells), `"T"` (T cells), `"Tgd"` (γδ T cells). Default is `"T"`.
#' @param col_res Character. Either `"full"` to retain all chain columns, or `"reduced"` to include only paired chains. Default is `"full"`.
#' @param add_columns Character vector. Additional metadata columns to include in the output. Optional.
#' @param cellranger_version Integer. Expected Cell Ranger version (used for validation): `3` (18 columns) or `7` (31 columns) or `9` (32 columns). Default is `7`.

#'
#' @return Annotated clonality matrix with clonal group information, invisibly returned after internal call to `assigntenx()`.
#'
#' @examples
#' tenx(data = Cellranger3_TCR, method = "unique_paired", cell = "T", clonality_input = c(mismatch = 20), cellranger_version = 3, )
#' tenx(data = Cellranger7_BCR1, method = "sticky_ends", cell = "B", only_productive = TRUE)
#'
#' @import dplyr
#' @importFrom stringr str_extract
#' @importFrom openxlsx write.xlsx
#' @importFrom data.table rbindlist
#' @export
tenx <- function(data = NULL,
                 method = "unique_paired",
                 only_productive = TRUE,
                 only_true_cells = TRUE,
                 clonality_input = NULL,
                 cell = "B",
                 col_res = c("full"),
                 add_columns = NULL,
                 cellranger_version = 7) {

  # Validate method
  if (!method %in% c("unique_paired", "sticky_ends", "unique_all")) {
    stop("Method chosen not valid. Choose one of: unique_paired, sticky_ends, unique_all.")
  }

  # Use example fallback if NULL
  if (is.null(data)) {
    warning("No input. Running on example data.")
    data <- Cellranger7_BCR1
  }

  # Validate and load data using the parser
  data <- cellranger_version_parser(data, cellranger_version = cellranger_version)

  # Remove empty or NA CDR3s
  data <- data[data$cdr3_nt != "" & !is.na(data$cdr3_nt), ]

  # Apply shared chain filters (productive / true cells) for all cell types
  if (only_productive) {
    data <- dplyr::filter(data, grepl("true", productive, ignore.case = TRUE))
  }

  if (only_true_cells) {
    data <- dplyr::filter(data, grepl("true", is_cell, ignore.case = TRUE))
  }

  # Special handling for gamma-delta T cells
  if (cell == "Tgd") {
    data <- data %>%
      filter(chain != "None", !grepl("\\*", cdr3))
  }

  # Split by barcode and organize chains
  data.list <- split(data, data[["barcode"]])

  sort.list <- function(x) x[order(x[["chain"]]), ]
  group.chains <- function(x) paste(x[["chain"]], collapse = "_")

  # Sort chains to sync all data.frames
  data.list <- lapply(data.list, sort.list)

  # Extract chain classes per barcode
  pairs <- data.frame(Chains = unlist(lapply(data.list, group.chains)))

  #Define empty variable
  list.pairs <- list()

  # Group chains from the same chain class together
  for (i in unique(pairs$Chains)) {
    list.pairs[[i]] <- data.list[pairs$Chains == i]
  }

  # Run clonality assignment per chain class
  assigntenx(
    list.pairs = list.pairs,
    method = method,
    clonality_input = clonality_input,
    cell = cell,
    col_res = col_res,
    add_columns = add_columns
  )
}
