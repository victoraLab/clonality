#' Annotate 10x clones
#'
#' Returns the input with clonality annotation column with clonal definitions.
#'
#' @param data Character. Data frame object or the full path to a filtered_contig_annotations.csv file.
#' @param method Character. One of: unique_paired, unique_all, sticky_ends. Default: unique_paired.
#' @param only_productive Logical. Filter non productive chains.
#' @param only_true_cells Logical. Filter low quality / non true cells. Using the default `TRUE` with the all_contig_annotations.csv file is equivalent of just running on the filtered_contig_annotations.csv.
#' @param clonality_input Named vector. Input parameters for the clonality function.
#' @param cell Character. Possible values: `B` Bcells, `T` Tcells, `Tgd` Gamma Delta T cells.
#' @param col_res Character. Possible values: `full` paired and unique chain columns, `reduced` only paired chain columns.
#' @param add_columns Character. Vector containing columns to add into the final result.
#' @param save_files Logical. Whether to save the Cl matrices as xlsx files or not.
#' @examples
#' tenx(data = "filtered_contig_annotations", method = "sticky_ends", only_productive = T, clonality_input = c("mismatch" = 0.25), cell = "T",  save.files = F)
#' @import dplyr
#' @importFrom stringr str_extract
#' @importFrom openxlsx write.xlsx
#' @importFrom data.table rbindlist
#' @export

tenx <- function(data = NULL,
                 method = "unique_paired",
                 only_productive = T,
                 only_true_cells = T,
                 clonality_input = NULL,
                 cell = "T",
                 col_res = c("full"),
                 save_files = F,
                 add_columns = NULL) {

    # Test parameters for correct input
    if( !any(method %in% c("unique_paired", "sticky_ends", "unique_all")) ){
      stop('Method chosen not valid.
           Choose one of: unique_paired, sticky_ends, unique_all.')
      }

    if(is.null(data) == T ){
        warning("No input. Running on example data.")
        data <- filtered_contig_annotations
    }



    if(!is.data.frame(data)){
        # Extract extension
        ext <- gsub("^.*\\.", "", data)

        if(!is.na(ext)){
            data <- read.csv(data)
        }else{
            data <- get(data)
        }
    }

  #Remove empty but not NA, CDR3s
  data <- data[data$cdr3_nt != "",]

  #Accept only productive chains
    if(only_productive == T){
      data <- data %>%
        filter(grepl("true", productive, ignore.case = TRUE), ignore.case = TRUE)
    }

  if(only_true_cells == T){
    data <- data %>%
      filter(grepl("true", is_cell, ignore.case = TRUE), ignore.case = TRUE)
  }


    # Split the 10x dataframe based on each barcode.
    data.list <- split(data, f = data[["barcode"]])

    # This function sorts barcodes to align paired chains.
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

    assigntenx(list.pairs = list.pairs, method = method, clonality_input = clonality_input, cell = cell, col_res = col_res, save_files = save_files, add_columns = add_columns)

}
