read_repertoire <- function(data, data_format = "general") {
  if (is.data.frame(data)) return(data)

  if (!is.character(data)) stop("Input must be a data frame or file path.")

  switch(data_format,
         "imgt" = parse_imgt(data),
         "general" = parse_general(data),
         stop("Unsupported format: choose 'imgt', '10x', or 'general'")
  )
}
