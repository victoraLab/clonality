#'
#' @importFrom tidyr pivot_wider
#' @export
#'
row1 <- function(x) {
  selected <- c("v_gene", "j_gene", "cdr3", "cdr3_nt", "chain", "raw_clonotype_id", "barcode")
  x <- x %>% select(selected)
  x[["chain"]] <- make.unique(x[["chain"]])
  x[["barcode"]] <- make.unique(x[["barcode"]])
  x %>%
    mutate(rn = 1) %>%
    pivot_wider(names_from = "chain", values_from = selected) %>%
    select(-rn)
}
