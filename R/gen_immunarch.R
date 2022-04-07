#' Generates an immunarch compatible input
#'
#' Returns a dataset compatible with the immunarch package.
#'
#' @param seu Seurat with clonality generated columns.
#' @param proportion Character. Possible values: `all` against all cells, `rm.na` remove NAs prior to proportion.
#' @param sample_column Character. Column to be used as sample to cut the data.
#' @param metadata_columns Character. Columns to split the data.
#' @examples
#' gen_immunarch(seu, proportion = "rm.na",  sample_column = "Conditions")
#' @import dplyr
#' @importFrom gtools mixedsort
#' @importFrom tidyr separate
#' @exports
#'

gen_immunarch <- function(seu = seu, proportion = "rm.na",  sample_column = NULL, metadata_columns = NULL){

  if(is.null(metadata_columns)){metadata_columns <- sample_column}
  if(is.null(sample_column)){stop("Define a column to be used as sample")}

  if(any(grepl("^Res", colnames(seu@meta.data))) == FALSE){
    stop("No clonality columns found. Run tenx and add_seurat_metadata first.", call. = FALSE)
  }

  cond<- sample_column
  if(sample_column != metadata_columns){
   newcond <- paste0(seu@meta.data[[sample_column]], "_", seu@meta.data[[metadata_columns]])

   seu@meta.data[[  "newcond"]] <- newcond
   cond <- "newcond"
   }

  meta.data <- seu@meta.data[, grepl("^Res", colnames(seu@meta.data))]

  meta <- split(meta.data, f = seu@meta.data[[cond]])

  prop1 <- function(x){
    x %>%
      filter(!is.na(Res_clonality)) %>%
      count(Res_clonality) %>%
      rename("Clones" = n) %>%
      mutate(Proportion = Clones/sum(Clones)) %>%
      arrange(desc(Clones)) %>%
      select(Clones, Proportion, Res_clonality)
  }

  prop2 <- function(x){
    x %>%
      count(Res_clonality) %>%
      rename("Clones" = n) %>%
      mutate(Proportion = Clones/sum(Clones)) %>%
      arrange(desc(Clones)) %>%
      select(Clones, Proportion, Res_clonality)
  }


  clono.counts <- switch (proportion,
                          rm.na = lapply(meta, prop1),
                          all = lapply(meta, prop2)
                          )

  meta_pairs <- meta.data %>% filter(Res_clonality.correction != "Corrected")


  gen_emoonarc <- function(x){

    clone.idx <- match(x[["Res_clonality"]], meta_pairs$Res_clonality)

    clono.meta <- meta_pairs[clone.idx,]

    clono.meta <- as.data.frame(apply(clono.meta, MARGIN = 2, FUN = function(x){gsub(";.*","",x)}))
    clono.meta <- as.data.frame(apply(clono.meta, MARGIN = 2, FUN = function(x){gsub("_",";", x)}))

    x <- x %>% mutate(CDR3.nt = clono.meta$Res_CDR3)
    x <- x %>% mutate(CDR3.aa = clono.meta$Res_cdr3_col2)
    x <- x %>% mutate(CDR3.aa = clono.meta$Res_cdr3_col2)
    x <- x %>% mutate(V.name  = clono.meta$Res_v_genes)
    x <- x %>% mutate(D.name  = NA)
    x <- x %>% mutate(V.end   = NA)
    x <- x %>% mutate(D.start = NA)
    x <- x %>% mutate(J.start = NA)
    x <- x %>% mutate(VJ.ins  = NA)
    x <- x %>% mutate(VD.ins  = NA)
    x <- x %>% mutate(DJ.ins  = NA)
    x <- x %>% mutate(Sequence  = clono.meta$Res_CDR3)
    x <- x %>% mutate(chain  = gsub("(.*)\\;(.*)", "\\1", clono.meta$Res_clonality))
    x <- x %>% mutate(Barcode  = NA)
    x <- x %>% mutate(raw_clonotype_id  = x$Res_clonality)
    x <- x %>% mutate(ContigID  = NA)

  }

  clono.counts <- lapply(clono.counts, FUN = gen_emoonarc)
  seu <- SetIdent(seu, value = cond)

  for(name in names(clono.counts)){
    x <- clono.counts[[name]]
    seu.subset <- subset(seu, idents = name)

    idx <- 0
    for(i in x[["Res_clonality"]]){
      idx <- idx + 1
      x[["Barcode"]][idx] <- paste(seu.subset@meta.data %>% filter(Res_clonality == i) %>% pull(Res_barcodes), collapse=";")
    }

    clono.counts[[name]] <- x

  }

  df <- data.frame(
    Sample = mixedsort(as.character(unique(seu@meta.data[[cond]])))
)

  if(sample_column != metadata_columns){
    df2 <- separate(data = df, col = Sample, into = c("Cond1", "Cond2"), sep = "_")
    df <- cbind(df,df2)
  }



immdata <- list(data = clono.counts,
     meta = df)

return(immdata)


}






