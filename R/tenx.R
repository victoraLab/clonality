#' Define tenx
#'
#' Returns the input with an additional clonality column with clonal definitions.
#'
#' @param File Data frame object or the full path to a xlsx input containing TCR/BCR repertoire data.
#' @param method One of: unique_paired, unique_all, sticky_ends. Default: unique_paired.

#' @examples
#' clonality(File = tra)

#' @import tidyverse

#' @export


tenx <- function (File = "example.xlsx",
                  method = "unique_paired"){

#Here we  split the 10x dataframe based on each barcode.
File.list <- split( File , f = bar)

#This function sorts every barcode table for the chain, to align every cell chains.
sort.list <- function(x){
  return(x[order(x[["chain"]]),])
  }

#This function creates a tag with the paired type for each cell barcode
group.chains <- function(x){
  paste(x[["chain"]], collapse = "_")
}

row1 <- function(x){
  x <- x %>% select(c(v_gene, j_gene, cdr3, cdr3_nt, raw_clonotype_id, chain, barcode))
  x[["chain"]] <- make.unique(x[["chain"]])
  x[["barcode"]] <- make.unique(x[["barcode"]])
  x %>% mutate(rn = 1) %>% pivot_wider(names_from = 'chain', values_from = c(v_gene, j_gene, cdr3, cdr3_nt, raw_clonotype_id, barcode)) %>% select(-rn)
}


#We apply sort.list function into our input
File.list <- lapply(File.list, sort.list)
#Create the pair tag
pairs <- data.frame(Chains = unlist(lapply(File.list, group.chains)))

list.pairs <- list()
for(i in unique(pairs)[["Chains"]]){
  list.pairs[[i]] <- File.list[pairs[["Chains"]] == i]
}
#

if(method == "unique_paired"){

  list.pairs <- list.pairs[grep("_", names(list.pairs))]
  res <- do.call(c, list.pairs)
  res <- lapply(res, row1)
  res <- rbind.fill(res)

  res <- res %>%
    mutate_at(vars(matches("raw")), as.character) %>%
    mutate(sc.raw_clonotypes = coalesce(!!! select(., matches("raw"))))

  res <- res %>%
    mutate_at(vars(matches("bar")), as.character) %>%
    mutate(sc.barcodes = coalesce(!!! select(., matches("bar"))))

  library(tidyverse)
  res <- res %>% select(-contains("raw_clonotype_id"))
  res <- res %>% select(-contains("barcode_"))

  clonality(File = res,
            Vgene_Column = "v_gene_IGH",
            Jgene_Column = "j_gene_IGH",
            CDR3_Column = "cdr3_nt",
            Cell = "B",
            ID_Column = "sc.barcodes"
            )


 }

if(method == "unique_all"){


}

if(method == "sticky_ends"){

}









}
