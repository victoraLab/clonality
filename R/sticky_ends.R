#' Define clonality for 10x genomics
#'
#' Assign clonality on 10x VDJ filtered contigs using sticky_ends method
#'
#' @param list.pairs
#' @param classes
#' @param pairing
#' @param clonality_input
#'
#' @export

sticky_ends <- function(list.pairs = list.pairs, classes = NULL, pairing = "_", clonality_input = clonality_input){


  if(pairing == "_"){
    classes <- classes[classes %in% c("TRA", "TRB", "TRA_TRB")]
  }

  else if(pairing == "__"){
    classes <- classes[classes %in% c("TRA", "TRB", "TRA_TRB", "TRA_TRA", "TRB_TRB", "TRA_TRA_TRB_TRB")]
  }

  list.pairs_filt <- list.pairs[names(list.pairs) %in% classes]

  res <- do.call(c, list.pairs_filt)
  res <- lapply(res, row1)

  res <- bind_rows(res, .id = "classes")

  res$classes <- gsub("\\..*", "", res$classes)

  res <- res %>% mutate_at(vars(matches("raw")), as.character) %>% mutate(sc.raw_clonotypes = coalesce(!!!select(., matches("raw"))))
  res <- res %>% mutate_at(vars(matches("bar")), as.character) %>% mutate(sc.barcodes = coalesce(!!!select(., matches("bar"))))

  res <- res %>% select(-contains("raw_clonotype_id"))
  res <- res %>% select(-contains("barcode_"))

  #res <- res %>% group_by(classes) %>% filter(n() > classes_n)

  for (i in unique(res$classes)) {

    res.sub <- res %>% filter(classes == i)

    cols <- str_split_fixed(i, "_", n = Inf)

    barcodes <- res.sub$sc.barcodes
    raw_clonotypes <- res.sub$sc.raw_clonotypes

    v_gene <- res.sub[, sprintf("v_gene_%s", cols)] %>% tidyr::unite("v_gene", sep = "_") %>% pull(v_gene)
    j_gene <- res.sub[, sprintf("j_gene_%s", cols)] %>% tidyr::unite("j_gene", sep = "_") %>% pull(j_gene)
    cdr3_col <- res.sub[, sprintf("%s_%s", "cdr3_nt", cols)] %>% tidyr::unite("cdr3_col", sep = "_") %>% pull(cdr3_col)
    cdr3_col2 <- res.sub[, sprintf("%s_%s", "cdr3", cols)] %>% tidyr::unite("cdr3_col2", sep = "_") %>% pull(cdr3_col2)
    cdr3_length <- as.data.frame(apply(res.sub[, sprintf("%s_%s", "cdr3_nt", cols)], MARGIN = 2, FUN = nchar)) %>% tidyr::unite("cdr3_length",
                                                                                                                                sep = "_")


    df1 <- data.frame(barcodes = barcodes, v_genes = v_gene,  j_genes = j_gene, CDR3 = cdr3_col,
                      cdr3_col2 = cdr3_col2, cdr3_length = cdr3_length, raw_clonotypes = raw_clonotypes)

    df1$v_genes <- gsub("\\+",";",df1$v_genes)
    df1$j_genes <- gsub("\\+",";",df1$j_genes)

    if(length(clonality_input) == 0){
      clonality_input <- c(output = "Clonal.output.10x", vgene_col = "v_genes", jgene_col = "j_genes", cdr3_col = "CDR3",
                           cell = "T", output_original = T,  id_col = "barcodes", mm = 0, search_gname = F)
    }

    else{

      input <- c(output = "Clonal.output.10x", vgene_col = "v_genes", jgene_col = "j_genes", cdr3_col = "CDR3",
                 cell = "T", output_original = T,  id_col = "barcodes", mm = 0, search_gname = F)

      input[names(clonality_input)] <- clonality_input

      clonality_input <- input
    }


    clonality(data = df1,
              output = paste0(clonality_input["output"], i),
              vgene_col = clonality_input["vgene_col"],
              jgene_col = clonality_input["jgene_col"],
              cdr3_col = clonality_input["cdr3_col"],
              cell = clonality_input["cell"],
              output_original = as.logical(clonality_input["output_original"]),
              suffix = i,
              id_col = clonality_input["id_col"],
              mm = as.numeric(clonality_input["mm"]),
              search_gname = as.logical(clonality_input["search_gname"]))



  }


}
