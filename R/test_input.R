#' Test input
#'
#' Evaluate if input data have the required artributes.
#'
#' @param test_data A data frame object or the full path to a xlsx with repertoire data.
#'

test_input <- function(test_data = NULL){

  if(!all(class(test_data) %in% "character")) {
    if("data.frame" %in% class(test_data)){
      imported_df <- test_data
    }
    else{
      stop("data is not a data.frame", call. = FALSE)
    }
  } else {
    #Read xlsx or csv files from text input
    ext <- str_extract(string = test_data, pattern = "\\..*$")

    if (any(c(".xlsx", ".csv")  %in% ext)){
      matched.ext <- match(ext, c(".xlsx", ".csv"))
      switch(matched.ext,
             imported_df <- read_excel(test_data),
             imported_df <- as_tibble(read.csv(test_data)))
    } else{
      stop("Import data in csv or xlsx", call. = FALSE)
    }
    return(imported_df)
  }
}
