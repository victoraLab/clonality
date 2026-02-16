#' Define 10x clonotypes
#'
#' Takes a list of paired chains and defines clonotypes.
#'
#' @param seu Seurat. A seurat object to embed with clonality data.
#' @param clonal.list Character. The clonalities to be inserted into the seurat metadata. Default:  `ls(pattern = "^Clonal")`.
#' @param cdr3_map Character. The columun name to be used to map single chains to paired chains. CDR3 for nt ; cdr3_col2 for aa. Default: "cdr3_col2".
#' @param purity Numeric. Minimum fraction of the dominant mapped paired to be used in a sticky assignment.
#' @param sticky Logical. If the script should merge the clonality of single chain cells with paired cells.
#' @param stick_only Character vector. Possible values: `TRA`, `TRB`, `TRD`, `TRG`, `IGH`, `IGK`, `IGL`, `all`. Default: `all`
#' @param cell Character. Possible values: `B` Bcells, `T` Tcells, `Tgd` Tcells GamaDelta.
#' @param overwrite Logical. Whether to overwrite clonality columns from seurat.

#' @importFrom stringr str_count
#' @export
#'
#'
add_seurat_metadata <- function(seu = NULL,
                                clonal.list = ls(pattern = "^Clonal", envir = .GlobalEnv),
                                cdr3_map = "cdr3_col2",     # now accepts c("cdr3_col2","CDR3")
                                cell = "T",
                                purity = 0.8,
                                stick_only = "all",
                                sticky = FALSE,
                                overwrite = TRUE){

  if (is.null(seu)) stop("Please add Seurat object to merge clonality.")

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("This function requires the 'Seurat' package. Install it first.", call. = FALSE)
  }

  # ensure vector
  cdr3_map <- as.character(cdr3_map)

  # ... (unchanged preamble) ...

  # mapping helper now takes a map column name
  map_singlets <- function(x, map_col) {
    mapped <- c2[grepl(x[[map_col]], c2[[map_col]], fixed = TRUE), ]
    if (nrow(mapped) == 0) return(mapped)
    # ensure V/J/length consistency
    if (length(table(mapped[["v_genes"]])) == 1 &&
        length(table(mapped[["j_genes"]])) == 1 &&
        length(table(mapped[["cdr3_length"]])) == 1) {
      return(mapped)
    }
    mapped[0, ] # empty
  }

  get_mapped.singles <- function(x, purity) {
    if (nrow(x) == 0) return(FALSE)
    tb <- table(x[["clonality"]])
    if (length(tb) == 1) return(TRUE)
    max(tb) / sum(tb) > purity
  }

  get_fixed_clonality <- function(x) names(which.max(table(x[["clonality"]])))

  #This function maps each singlet to a paired
  map_singlets <- function(x) {
    mapped <- c2[grepl(x[[cdr3_map]], c2[[cdr3_map]], fixed = T),]
    if(nrow(mapped) == 0){}
    if(length(table(mapped[["v_genes"]])) == 1 ){
      if(length(table(mapped[["j_genes"]])) == 1 ){
        if(length(table(mapped[["v_genes"]])) == 1 ){
          if(length(table(mapped[["cdr3_length"]])) == 1 ){
            return(mapped)
          }else{}
        }else{}
      }else{}
    }else{}
    }

  #Get single chain cells to map
  singles <- c(clonal.list[grep("_", clonal.list, invert = T)])

  #Filter single chains selected by the user
  if(stick_only != "all"){
    singles <- singles[grepl(stick_only, singles)]
    if(length(singles) == 0){stop(sprintf("No detected %s", stick_only))}
  }

  pair_match <- switch(cell, "B" =   c("Clonal.output.10xIGH_IGK", "Clonal.output.10xIGH_IGL"),
                             "T" =   c("Clonal.output.10xTRA_TRB"),
                             "Tgd" = c("Clonal.output.10xTRD_TRG"))

  pairs   <- c(clonal.list[grep("^^[^_]*_[^_]*$", clonal.list)])
  pairs   <- pairs[grepl(pair_match, pairs)]

  if(length(pairs) == 0){stop(print("No detected %s", pairs))}

  rest <- clonal.list[!clonal.list %in% c(singles,pairs)]

  c1 <- do.call(bind_rows, mget(singles, envir = .GlobalEnv))
  c2 <- do.call(bind_rows, mget(pairs,   envir = .GlobalEnv))
  c3 <- do.call(bind_rows, mget(rest,    envir = .GlobalEnv))

  if (sticky) {
    # try maps in priority order
    tried <- list()
    success <- FALSE
    for (mc in cdr3_map) {
      c1.mapped <- apply(c1, MARGIN = 1, map_singlets, map_col = mc)
      names(c1.mapped) <- c1$barcodes

      ok <- vapply(c1.mapped, get_mapped.singles, logical(1), purity = purity)
      c1.correct <- c1.mapped[ok]

      if (length(c1.correct) > 0) {
        # success with this map
        df <- data.frame(names = names(c1.correct))
        df$clonality.corrected <- unlist(lapply(c1.correct, get_fixed_clonality))

        c1 <- dplyr::left_join(c1, df, by = c("barcodes" = "names"))
        c1$clonality.correction <- c("Corrected", "Normal")[
          as.numeric(coalesce(c1$clonality.corrected, c1$clonality) == c1$clonality) + 1
        ]
        c1$clonality <- coalesce(c1$clonality.corrected, c1$clonality)
        c1$clonality.corrected <- NULL
        c1$mapping_source <- mc
        success <- TRUE
        break
      } else {
        tried[[mc]] <- TRUE
      }
    }
    if (!success) {
      c1$mapping_source <- NA_character_
    }
  }

  # merge all clonotype tables
  c0 <- dplyr::bind_rows(c1, c2, c3)
  c0 <- as.data.frame(c0)
  colnames(c0) <- paste0("Res_", colnames(c0))
  rownames(c0) <- c0[["Res_barcodes"]]

  if (sticky) {
    c0$Res_clonality.correction[is.na(c0$Res_clonality.correction)] <- "Not tested"
    # carry mapping source if present
    if (!"Res_mapping_source" %in% colnames(c0) && "Res_mapping_source" %in% names(c0)) {
      # no-op; naming already handled by prefix above
    }
  }

  if (overwrite) {
    seu@meta.data <- seu@meta.data %>% dplyr::select(-dplyr::contains("Res_"))
  }

  seu <- Seurat::AddMetaData(seu, metadata = c0)
  return(seu)
}
