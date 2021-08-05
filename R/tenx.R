#' Define 10x clonotypes
#'
#' Returns the input with an additional clonality column with clonal definitions.
#'
#' @param data Data frame object or the full path to a filtered_contig_annotations.csv input.
#' @param method One of: unique_paired, unique_all, sticky_ends. Default: unique_paired. Character.
#' @param only_productive Filter non productive chains. Logical.
#' @param clonality_input Input parameters for the clonality function. List.
#' @param cell Character. Possible values: `B` Bcells, `T` Tcells, `Tgd` Tcells GamaDelta.

#' @examples
#' tenx(data = "filtered_contig_annotations", method = "sticky_ends", only_productive = T, clonality_input = c("mm" = 0.25), cell = "T",  save.files = F)
#' @import dplyr
#' @importFrom stringr str_extract
#' @importFrom data.table rbindlist
#' @export

tenx <- function(data = NULL, method = "unique_paired", only_productive = T, clonality_input = NULL, cell = "T", save.files = F) {

    # Test parameters for correct input
    if( !any(method %in% c("unique_paired", "sticky_ends", "unique_all")) ) stop('Method chosen not valid. Choose one of: unique_paired, sticky_ends, unique_all.')

    if(is.null(data) == T ){
        warning("No input. Running on example data.")
        data <- filtered_contig_annotations
    }

    # Extract extention
    ext <- str_extract(string = data, pattern = "\\..*") == ".csv"

    if(!is.data.frame(data)){
        if(!is.na(ext)){
            data <- read.csv(data)
        }else{
            data <- get(data)
        }
    }


    if(only_productive == T){
        data <- data %>% filter(grepl("true", productive, ignore.case = TRUE), ignore.case = TRUE)
    }

    # Split the 10x dataframe based on each barcode.
    data.list <- split(data, f = data[["barcode"]])

    # Sort barcodes to align paired chains.
    sort.list <- function(x) {
        return(x[order(x[["chain"]]), ])
    }

    # Creates a tag with the paired type for each cell barcode
    group.chains <- function(x) {
        paste(x[["chain"]], collapse = "_")
    }

    # Apply sort.list function
    data.list <- lapply(data.list, sort.list)

    # Create the pair tag
    pairs <- data.frame(Chains = unlist(lapply(data.list, group.chains)))

    list.pairs <- list()
    for (i in unique(pairs)[["Chains"]]) {
        list.pairs[[i]] <- data.list[pairs[["Chains"]] == i]
    }

    assigntenx(list.pairs = list.pairs, method = method, clonality_input = clonality_input, cell = cell)



}
