#' Define 10x clonotypes
#'
#' Takes a list of paired chains and defines clonotypes.
#'
#' @param seu Seurat. A seurat object to embed with clonality data.
#' @param clonal.list Character. The clonalities to be inserted into the seurat metadata. Default:  `ls(pattern = "^Clonal")`.
#' @param cdr3_map Character. The columun name to be used to map single chains to paired chains. CDR3 for nt ; cdr3_col2 for aa. Default: "cdr3_col2".
#' @param purity Numeric. Minimum fraction of the dominant mapped paired to be assigned to single.
#'

#' @importFrom stringr str_count
#' @export
#'
add_seurat_metadata <- function(seu = NULL, clonal.list = ls(pattern = "^Clonal", envir = .GlobalEnv), cdr3_map = "cdr3_col2", purity = 0.8){

  if(length(seu) == 0){stop("Please add Seurat object to merge clonality.")}

  get_mapped.singles <- function(x){
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


  get_fixed_clonality <- function(x){

    return(names(which.max(table(x[["clonality"]]))))
  }

  singles <- c(clonal.list[grep("_", clonal.list, invert = T)])
  pairs   <- c(clonal.list[grep("^^[^_]*_[^_]*$", clonal.list)])
  multi <- clonal.list[!clonal.list %in% c(singles,pairs)]
 c1 <- do.call(bind_rows, mget(singles, envir = .GlobalEnv))
  c2 <- do.call(bind_rows, mget(pairs, envir = .GlobalEnv))
  c4 <- do.call(bind_rows, mget(multi, envir = .GlobalEnv))

  # which single chain, matches the paired cdr3
  c1.mapped <- apply(c1, 1 , function(x) c2[grepl(x[[cdr3_map]], c2[[cdr3_map]]),])

  names(c1.mapped) <- c1$barcodes
  c1.filter <- lapply(c1.mapped, FUN = get_mapped.singles)
  c1.correct <- c1.mapped[unlist(c1.filter)]
  df <- data.frame(names = names(c1.correct))
  df$clonality.corrected <- NA
  df$clonality.corrected <- unlist(lapply(c1.correct, FUN = get_fixed_clonality))

  c1 <- left_join(c1, df, by = c("barcodes" = "names"))

  c1$clonality <- coalesce(c1$clonality.corrected, c1$clonality)
  c1$clonality.corrected <- NULL


  c3 <- do.call(bind_rows, mget(setdiff(ls(pattern = "^Clonal"), clonal.list)))

  c0 <- bind_rows(c1,c2,c3,c4)

  meta <- seu@meta.data
  meta$barcodes <- rownames(meta)
  meta <- left_join(meta, c0,  by = "barcodes")
  rownames(meta) <- meta$barcodes


  seu@meta.data <- meta

  return(seu)


}
