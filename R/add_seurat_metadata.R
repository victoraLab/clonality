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
                                cdr3_map = "cdr3_col2",
                                cell = "T",
                                purity = 0.8,
                                stick_only = "all",
                                sticky = FALSE,
                                overwrite = T){

  #Seu parameter is mandatory
  if(is.null(seu)){stop("Please add Seurat object to merge clonality.")}

  #Function to select cells to map
  get_mapped.singles <- function(x){
    if(length(x) == 0){return(FALSE)}
    if(nrow(x) == 0){
      return(FALSE)
    } else {
      if(nrow(x) > 0){
        if(length(table(x[["clonality"]])) > 1){
          if(max(table(x[["clonality"]]))/nrow(x) > purity){
            return(TRUE)
          } else {
            return(FALSE)
          }
        }
      }

      if(length(table(x[["clonality"]])) == 1){
        return(TRUE)
      }

    }

  }

  #This function gives the top matched paired clone
  get_fixed_clonality <- function(x){
    return(names(which.max(table(x[["clonality"]]))))
  }

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

  if(sticky == T){

  #Which single chain, matches the paired cdr3 sequences
  c1.mapped <- apply(c1, MARGIN = 1 , map_singlets)
  names(c1.mapped) <- c1$barcodes

  #Filter out non mapping clones
  c1.filter <- lapply(c1.mapped, FUN = get_mapped.singles)
  c1.correct <- c1.mapped[unlist(c1.filter)]

  #Get the top matching clone
  df <- data.frame(names = names(c1.correct))
  df$clonality.corrected <- NA
  df$clonality.corrected <- unlist(lapply(c1.correct, FUN = get_fixed_clonality))

  #Replace old IDs
  c1 <- left_join(c1, df, by = c("barcodes" = "names"))
  c1$clonality.correction <- c("Corrected", "Normal")[as.numeric(coalesce(c1$clonality.corrected, c1$clonality) == c1$clonality) + 1]

  c1$clonality <- coalesce(c1$clonality.corrected, c1$clonality)
  c1$clonality.corrected <- NULL

}

  #Merge all clonotype tables
  c0 <- bind_rows(c1,c2,c3)
  colnames(c0) <- paste0("Res_", colnames(c0))
  rownames(c0) <- c0[["Res_barcodes"]]

  if(sticky == T){
  c0$Res_clonality.correction[is.na(c0$Res_clonality.correction)] <- "Not tested"
  }

  if(overwrite == T){
    seu@meta.data <- seu@meta.data %>% select(-contains("Res_"))
  }

  seu <- AddMetaData(seu, metadata = c0)

  return(seu)

}
