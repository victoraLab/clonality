#' Annotate 10x clones with clonality information
#'
#' Loads and processes 10x Genomics V(D)J output (either as a data frame, file path, or object name), filters and groups by cell barcode, and assigns clonal definitions.
#'
#' Supports BCR, TCR, and γδ TCR data from Cell Ranger versions 3 (18 columns) and 7 (31 columns).
#'
#' @param data Data frame, character string (object name), or full path to a `filtered_contig_annotations.csv` or `all_contig_annotations.csv` file.
#' @param experiment Character. Experiment name ID `optional`.
#' @param method Character. Clonality inference method. One of: `"unique_paired"`, `"sticky_ends"`, `"unique_all"`. Default is `"unique_paired"`.
#' @param only_productive Logical. Whether to filter non-productive chains. Default is `TRUE`.
#' @param only_true_cells Logical. Whether to filter out low-quality/non-true cells (i.e., `is_cell != TRUE`). Default is `TRUE`.
#' @param clonality_input Named vector. Parameters passed to the internal `clonality()` function (e.g., mismatch thresholds).
#' @param cell Character. Cell type to process. One of: `"B"` (B cells), `"T"` (T cells), `"Tgd"` (γδ T cells). Default is `"T"`.
#' @param col_res Character. Either `"full"` to retain all chain columns, or `"reduced"` to include only paired chains. Default is `"full"`.
#' @param add_columns Character vector. Additional metadata columns to include in the output. Optional.
#' @param airr_normalize Character. One of `"auto"` (default), `"always"`, `"never"`.
#'   Controls whether AIRR / IgBLAST quirks (comma-separated alleles, logical-form
#'   `productive`) are normalized via [normalize_airr_to_cellranger()] before
#'   column renaming. `"auto"` triggers normalization only when AIRR-style
#'   patterns are detected.
#' @inheritParams cellranger_version_parser
#' @return Annotated clonality matrix with clonal group information, invisibly returned after internal call to `assigntenx()`.
#'
#' @examples
#' \dontrun{
#' Cellranger7_BCR1 <- tenx(data = Cellranger7_BCR1, method = "sticky_ends", cell = "B")
#'}
#' @import dplyr
#' @importFrom stringr str_extract
#' @importFrom openxlsx write.xlsx
#' @export
tenx <- function(data = NULL,
                 experiment = "experiment1",
                 method = "unique_paired",
                 only_productive = TRUE,
                 only_true_cells = TRUE,
                 clonality_input = NULL,
                 cell = "B",
                 col_res = c("full"),
                 add_columns = NULL,
                 cellranger_version = NULL,
                 ignore_version = FALSE,
                 airr_normalize = c("auto", "always", "never")) {

  # Validate method
  if (!method %in% c("unique_paired", "sticky_ends", "unique_all")) {
    stop("Method chosen not valid. Choose one of: unique_paired, sticky_ends, unique_all.")
  }

  airr_normalize <- match.arg(airr_normalize)

  # Use example fallback if NULL
  if (is.null(data)) {
    warning("No input. Running on example data.")
    data <- Cellranger7_BCR1
  }

  # Validate and load data using the parser
  data <- cellranger_version_parser(data,
                                    cellranger_version = cellranger_version,
                                    ignore_version = ignore_version)

  # Optionally normalize AIRR / IgBLAST quirks (comma-separated alleles,
  # logical-form productive flag) before column renaming.
  if (airr_normalize == "always" ||
      (airr_normalize == "auto" && .looks_like_airr(data))) {
    data <- normalize_airr_to_cellranger(data)
  }

  data <- normalize_tenx_columns(data)

  if (!"cdr3_nt" %in% names(data) && "cdr3" %in% names(data)) {
    data$cdr3_nt <- data$cdr3
    warning(
      "Input does not contain `cdr3_nt`; falling back to amino-acid CDR3 values for clonotype assignment.",
      call. = FALSE
    )
  }

  required_cols <- c("barcode", "chain", "cdr3_nt", "v_gene", "j_gene")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Input is missing required 10x/AIRR columns after normalization: ",
      paste(missing_cols, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # Remove empty or NA CDR3s
  data <- data[data$cdr3_nt != "" & !is.na(data$cdr3_nt), ]

  if (!nrow(data)) {
    stop(
      "No rows remain after preprocessing. Check that your input includes non-missing CDR3 nucleotide data ",
      "(`cdr3_nt` or an alias such as `junction`) and that any requested filters are not removing all rows.",
      call. = FALSE
    )
  }

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
  pairs <- data.frame(
    Chains = as.character(unlist(lapply(data.list, group.chains))),
    stringsAsFactors = FALSE
  )

  #Define empty variable
  list.pairs <- list()

  # Group chains from the same chain class together
  for (i in unique(pairs$Chains)) {
    class_name <- as.character(i)
    list.pairs[[class_name]] <- data.list[pairs$Chains == class_name]
  }

  # Run clonality assignment per chain class
  out_list <- assigntenx(
              list.pairs = list.pairs,
              experiment = experiment,
              method = method,
              clonality_input = clonality_input,
              cell = cell,
              col_res = col_res,
              add_columns = add_columns
              )

  obj_names <- paste0("Clonal.output.10x", names(out_list))
  list2env(stats::setNames(out_list, obj_names), envir = .GlobalEnv)
  return(out_list)
}

# Heuristic: looks like AIRR/IgBLAST output if it has v_call/j_call columns or
# if comma-separated alleles appear in any V/J/C call column.
.looks_like_airr <- function(data) {
  airr_cols <- c("v_call", "j_call", "c_call", "d_call")
  has_airr_cols <- any(airr_cols %in% names(data))
  comma_cols <- intersect(c("v_gene", "d_gene", "j_gene", "c_gene", airr_cols), names(data))
  has_commas <- any(vapply(comma_cols, function(col) {
    any(grepl(",", as.character(data[[col]]), fixed = TRUE), na.rm = TRUE)
  }, logical(1)))
  has_logical_productive <- "productive" %in% names(data) &&
    is.logical(data[["productive"]])
  has_airr_cols || has_commas || has_logical_productive
}

normalize_tenx_columns <- function(data) {
  aliases <- list(
    barcode = c("barcode", "barcodes", "cell_id", "cell", "cellid"),
    chain = c("chain", "locus"),
    cdr3_nt = c("cdr3_nt", "junction"),
    cdr3 = c("cdr3", "cdr3_aa", "junction_aa"),
    v_gene = c("v_gene", "v_call", "v_call_genotyped"),
    j_gene = c("j_gene", "j_call"),
    c_gene = c("c_gene", "c_call", "c_call_genotyped"),
    raw_clonotype_id = c("raw_clonotype_id", "clone_id", "clonotype_id"),
    productive = c("productive"),
    is_cell = c("is_cell")
  )

  for (target in names(aliases)) {
    if (target %in% names(data)) {
      next
    }

    match <- aliases[[target]][aliases[[target]] %in% names(data)]
    if (length(match)) {
      data[[target]] <- data[[match[[1]]]]
    }
  }

  if (!"cdr3" %in% names(data) && "cdr3_nt" %in% names(data)) {
    data$cdr3 <- data$cdr3_nt
  }

  if (!"raw_clonotype_id" %in% names(data)) {
    data$raw_clonotype_id <- NA_character_
  }

  data
}
