#' @noRd
row1 <- function(x, add_columns) {
  selected <- c("v_gene", "j_gene", "c_gene", "cdr3", "cdr3_nt", "chain", "raw_clonotype_id", "barcode")
  if(!is.null(add_columns)){
    selected <- c(selected, add_columns)
    }
  x <- dplyr::select(x, dplyr::all_of(selected))
  x[["chain"]] <- make.unique(x[["chain"]])
  x[["barcode"]] <- make.unique(x[["barcode"]])

  x %>%
    dplyr::mutate(rn = 1) %>%
    tidyr::pivot_wider(names_from = "chain", values_from = dplyr::all_of(selected)) %>%
    dplyr::select(-rn)
}
