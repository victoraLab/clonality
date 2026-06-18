#' Annotate 10x clones with clonality information
#'
#' Loads and processes 10x Genomics V(D)J output (either as a data frame, file
#' path, or object name), filters and groups by cell barcode, and assigns clonal
#' definitions.
#'
#' When `data` is a **list or character vector of multiple inputs**, each element
#' is loaded and processed independently; barcodes are suffixed with the
#' corresponding `project_names` entry (`paste0(barcode, "_", project_name)`) to
#' ensure uniqueness across samples before the combined data are grouped and
#' clonotyped together. `project_names` is required in this mode.
#'
#' Supports BCR, TCR, and γδ TCR data from Cell Ranger versions 3 (18 columns)
#' and 7 (31 columns).
#'
#' @param data Data frame, character string (object name or path to a
#'   `filtered_contig_annotations.csv`), **or** a list / character vector of
#'   multiple such inputs for multi-sample ingestion.
#' @param project_names Character vector of unique sample labels, one per
#'   element of `data`. Required when `data` contains multiple inputs. Each
#'   barcode is suffixed with its label so that barcodes remain unique across
#'   samples. Ignored for single-input calls.
#' @param experiment Character. Experiment name ID `optional`.
#' @param method Character. Clonality inference method. One of:
#'   `"unique_paired"`, `"sticky_ends"`, `"unique_all"`. Default is
#'   `"unique_paired"`.
#' @param only_productive Logical. Whether to filter non-productive chains.
#'   Default is `TRUE`.
#' @param only_true_cells Logical. Whether to filter out low-quality/non-true
#'   cells (i.e., `is_cell != TRUE`). Default is `TRUE`.
#' @param clonality_input **Named list** of overrides for the internal call
#'   to [clonality()] that defines clones once cells have been filtered and
#'   grouped by chain class. `tenx()` builds a per-cell metadata table and
#'   then calls `clonality()` once per chain class with a set of sensible
#'   defaults; `clonality_input` lets you change any of those defaults.
#'
#'   It **must be a `list()`**, not a `c()` vector. A vector silently flattens
#'   any element that is itself a vector — so
#'   `clonality_input = c(cdr3_col = c("CDR3", "cdr3_col2"))` collapses into a
#'   single 2-element character vector named `cdr3_col1` / `cdr3_col2`, which
#'   is **not** what you want. Use `list(cdr3_col = c("CDR3", "cdr3_col2"))`
#'   instead. For scalar overrides (`list(mismatch = 10)`), either form works,
#'   but `list()` is recommended for consistency.
#'
#'   The slots that are commonly overridden are:
#'
#'   - `mismatch` — percent CDR3 mismatch tolerance (0--100). Default `0` =
#'     identical CDR3s only. Pass a vector when also passing a multi-element
#'     `cdr3_col` (see below).
#'   - `cdr3_col` — column(s) in the per-cell metadata to use as the CDR3 for
#'     clonality. Defaults to `"CDR3"` (concatenated nucleotide CDR3 across
#'     chains). Use `c("CDR3", "cdr3_col2")` to cluster on both nucleotide and
#'     amino-acid CDR3 in a single pass; the output then contains both
#'     `clonality_nt` and `clonality_aa` columns.
#'   - `stringdist_method` — string-distance method passed to
#'     `stringdist::stringdistmatrix()`. Default `"osa"`. Recycled to match
#'     `cdr3_col`.
#'   - `output_original` — `TRUE` (default) joins the clonality assignment back
#'     to the per-cell metadata so the output keeps barcodes, V/J genes,
#'     concatenated CDR3s, etc. Set `FALSE` for a minimal output.
#'   - `file` — path to an `.xlsx` file to save the per-chain-class result.
#'   - `project` — character label tagged onto each row's `project` column.
#'
#'   The following slots are managed internally by `tenx()` and overriding
#'   them is rarely useful (and may break the pipeline):
#'   `vgene_col`, `jgene_col`, `ident_col`, `cell`, `suffix`. They default to
#'   the column names that `tenx()` constructs in the per-cell metadata.
#'
#'   Examples:
#'   ```r
#'   # Default clones, but allow 10% CDR3 mismatch
#'   tenx(data, clonality_input = list(mismatch = 10))
#'
#'   # Cluster on both nucleotide and amino-acid CDR3, with their own
#'   # mismatch thresholds:
#'   tenx(data,
#'        clonality_input = list(cdr3_col = c("CDR3", "cdr3_col2"),
#'                               mismatch = c(0, 10)))
#'
#'   # Save the per-chain-class clonality table to disk:
#'   tenx(data, clonality_input = list(file = "clones.xlsx"))
#'   ```
#' @param cell Character. Cell type to process. One of: `"B"` (B cells), `"T"`
#'   (T cells), `"Tgd"` (γδ T cells). Default is `"T"`.
#' @param col_res Character. Either `"full"` to retain all chain columns, or
#'   `"reduced"` to include only paired chains. Default is `"full"`.
#' @param add_columns Character vector. Additional metadata columns to include
#'   in the output. Optional.
#' @param airr_normalize Character. One of `"auto"` (default), `"always"`,
#'   `"never"`. Controls whether AIRR / IgBLAST quirks (comma-separated
#'   alleles, logical-form `productive`) are normalized via
#'   [normalize_airr_to_cellranger()] before column renaming. `"auto"` triggers
#'   normalization only when AIRR-style patterns are detected.
#' @inheritParams cellranger_version_parser
#' @return Annotated clonality matrix with clonal group information, invisibly
#'   returned after internal call to `assigntenx()`.
#'
#' @examples
#' \dontrun{
#' # Single sample
#' out <- tenx(data = Cellranger7_BCR1, method = "sticky_ends", cell = "B")
#'
#' # Multi-sample: barcodes are suffixed with the project name
#' out <- tenx(
#'   data         = list("sample1.csv", "sample2.csv"),
#'   project_names = c("S1", "S2"),
#'   method       = "sticky_ends",
#'   cell         = "B"
#' )
#' }
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
                 airr_normalize = c("auto", "always", "never"),
                 project_names = NULL) {

  # Validate method
  if (!method %in% c("unique_paired", "sticky_ends", "unique_all")) {
    stop("Method chosen not valid. Choose one of: unique_paired, sticky_ends, unique_all.")
  }

  airr_normalize <- match.arg(airr_normalize)

  # clonality_input must be a list (or NULL). A bare named vector would
  # silently flatten any vector-valued slot like cdr3_col = c("CDR3", "cdr3_col2").
  if (!is.null(clonality_input) && !is.list(clonality_input)) {
    if (is.atomic(clonality_input) && !is.null(names(clonality_input))) {
      stop(
        "`clonality_input` must be a list(), not a c() vector. A vector ",
        "flattens nested values: e.g. c(cdr3_col = c(\"CDR3\", \"cdr3_col2\")) ",
        "becomes a single character vector, not a multi-CDR3 spec. ",
        "Use list(cdr3_col = c(\"CDR3\", \"cdr3_col2\"), mismatch = c(0, 10)) instead.",
        call. = FALSE
      )
    }
    stop("`clonality_input` must be a named list (or NULL).", call. = FALSE)
  }

  # ── Multi-sample support ───────────────────────────────────────────────────
  # Detect multiple inputs: a list with > 1 element, or a character vector with
  # > 1 element (multiple file paths). A single data frame is never multi.
  is_multi <- !is.data.frame(data) &&
              (is.list(data) || (is.character(data) && length(data) > 1)) &&
              length(data) > 1

  if (is_multi) {
    if (is.null(project_names)) {
      stop(
        "`project_names` must be provided when `data` contains multiple inputs.",
        call. = FALSE
      )
    }
    if (length(project_names) != length(data)) {
      stop(
        "`project_names` must be the same length as `data` (",
        length(data), " inputs, ", length(project_names), " names provided).",
        call. = FALSE
      )
    }
    if (anyDuplicated(project_names)) {
      stop("`project_names` must be unique.", call. = FALSE)
    }

    samples <- lapply(seq_along(data), function(i) {
      message("Loading sample ", i, "/", length(data), ": ", project_names[[i]])
      .tenx_prepare_sample(
        data               = data[[i]],
        cellranger_version = cellranger_version,
        ignore_version     = ignore_version,
        airr_normalize     = airr_normalize,
        only_productive    = only_productive,
        only_true_cells    = only_true_cells,
        cell               = cell,
        project_name       = project_names[[i]]
      )
    })
    data <- dplyr::bind_rows(samples)

  } else {
    # Single sample (existing behavior)
    if (is.null(data)) {
      warning("No input. Running on example data.")
      data <- Cellranger7_BCR1
    }

    data <- .tenx_prepare_sample(
      data               = data,
      cellranger_version = cellranger_version,
      ignore_version     = ignore_version,
      airr_normalize     = airr_normalize,
      only_productive    = only_productive,
      only_true_cells    = only_true_cells,
      cell               = cell,
      project_name       = NULL
    )
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

# Load, normalize, filter, and optionally suffix barcodes for one sample.
# Used by tenx() for both single- and multi-sample paths.
.tenx_prepare_sample <- function(data,
                                  cellranger_version,
                                  ignore_version,
                                  airr_normalize,
                                  only_productive,
                                  only_true_cells,
                                  cell,
                                  project_name = NULL) {

  data <- cellranger_version_parser(data,
                                    cellranger_version = cellranger_version,
                                    ignore_version = ignore_version)

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

  # Suffix barcodes with project name so they remain unique across samples
  if (!is.null(project_name)) {
    data$barcode <- paste0(data$barcode, "_", project_name)
  }

  data
}
