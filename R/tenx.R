#' Define 10x genomics clonotypes
#'
#' Returns the input with an additional clonality column with clonal definitions.
#'
#' @param data Data frame object or the full path to a filtered_contig_annotations.csv input.
#' @param method One of: unique_paired, unique_all, sticky_ends. Default: unique_paired. Character.
#' @param only_productive Filter non productive chains. Logical.
#' @param classes Types of chains to be annotated. Character.
#' @param clonality_input Input parameters for the clonality function. List.
#' @param pairing Character. _ for sticky_ends only.
#' @examples
#' clonality(data = filtered_contig_annotations)
#' @import tidyverse
#' @import dplyr
#' @importFrom tidyr unite
#' @importFrom stringr str_split_fixed
#' @importFrom data.table rbindlist
#' @export

tenx <- function(data = NULL, method = "unique_paired", only_productive = T, clonality_input = NULL, cell = "T", classes = "all", pairing = "_", save.files = F) {

    if( !any(method %in% c("unique_paired", "sticky_ends", "unique_all")) ) stop('Method chosen not valid. Choose one of: unique_paired, sticky_ends, unique_all.')

    if(is.null(data) == T ){
        warning("No input. Running on example data.")
        data <- filtered_contig_annotations
    }

    if(is.character(data) == T){
        data <- read.csv(data)
    }

    if(only_productive == T){
        data <- data %>% filter(grepl("true", productive, ignore.case = TRUE), ignore.case = TRUE)
    }

    # Split the 10x dataframe based on each barcode.
    data.list <- split(data, f = data[["barcode"]])

    # Sorts every barcode table for the chain, to align every cell chains.
    sort.list <- function(x) {
        return(x[order(x[["chain"]]), ])
    }

    # Creates a tag with the paired type for each cell barcode
    group.chains <- function(x) {
        paste(x[["chain"]], collapse = "_")
    }

    # We apply sort.list function
    data.list <- lapply(data.list, sort.list)

    # Create the pair tag
    pairs <- data.frame(Chains = unlist(lapply(data.list, group.chains)))

    list.pairs <- list()
    for (i in unique(pairs)[["Chains"]]) {
        list.pairs[[i]] <- data.list[pairs[["Chains"]] == i]
    }
    #

    if(classes == "all"){

        classes <- unique(list.pairs = list.pairs, names(list.pairs))

    }
    else{

        classes <- classes

    }

    if (method == "unique_paired") {

        unique_paired(list.pairs = list.pairs, classes = classes, clonality_input = clonality_input)

    }

    else if (method == "unique_all") {

        unique_all()

    }

    else if (method == "sticky_ends") {

        sticky_ends(list.pairs = list.pairs, classes = classes, pairing = pairing, clonality_input = clonality_input)

    }





}
