
#' @importFrom stringr str_count
#' @export
#'
add_seurat_metadata <- function(seu = NULL, clonal.list = ls(pattern = "^Clonal"), cell = "T", cdr3_map = "cdr3_col2"){

  get_mapped.singles <- function(x){
    if(nrow(x) == 0){
      return(FALSE)
    }
    else if(nrow(x) > 0){
      if(length(table(x[["clonality"]])) > 1){
        return(FALSE)
      }
      else if(table(x[["clonality"]]) == nrow(x)){
        return(TRUE)
      }
    }
  }

  get_fixed_clonality <- function(x){
    return(unique(x[["clonality"]]))
  }

  if(cell == "T"){
   cnt_a <- str_count(clonal.list, pattern = "TRA")
   cnt_b <- str_count(clonal.list, pattern = "TRB")

   singles <- c(clonal.list[(cnt_a == 0 | cnt_b == 0)])
   pairs   <- c(clonal.list[cnt_a > 0 & cnt_b > 0])

   c1 <- do.call(bind_rows, mget(singles, envir = .GlobalEnv))
   c2 <- do.call(bind_rows, mget(pairs, envir = .GlobalEnv))

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
   meta <- seu@meta.data
   meta$barcodes <- rownames(meta)

   c3 <- do.call(bind_rows, mget(setdiff(ls(pattern = "^Clonal"), clonal.list)))

   c0 <- bind_rows(c1,c2,c3)
   c0$barcodes <- str_extract(c0$barcodes, pattern = "[ATCG]{16}")
   meta <- left_join(meta, c0,  by = "barcodes")
   rownames(meta) <- meta$barcodes


   seu@meta.data <- meta

   return(seu)

  }

  if(cell == "B"){
    cnt_h <- str_count(clonal.list, pattern = "IGH")
    cnt_k <- str_count(clonal.list, pattern = "IGK")
    cnt_l <- str_count(clonal.list, pattern = "IGL")

    singles <- c(clonal.list[(cnt_h == 0 | cnt_k == 0)])
    pairs   <- c(clonal.list[cnt_h > 0 & cnt_k > 0])

    c1 <- do.call(bind_rows, mget(singles, envir = .GlobalEnv))
    c2 <- do.call(bind_rows, mget(pairs, envir = .GlobalEnv))

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
    meta <- seu@meta.data
    meta$barcodes <- rownames(meta)

    c3 <- do.call(bind_rows, mget(setdiff(ls(pattern = "^Clonal"), clonal.list)))

    c0 <- bind_rows(c1,c2,c3)
    c0$barcodes <- str_extract(c0$barcodes, pattern = "[ATCG]{16}")
    meta <- left_join(meta, c0,  by = "barcodes")
    rownames(meta) <- meta$barcodes


    seu@meta.data <- meta

    return(seu)
  }


}
