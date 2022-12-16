#' Annotate clonality
#'
#' The \pkg{clonality} package returns the input table with the clonality annotations.
#'
#' @param data A data frame object or the full path to a xlsx with repertoire data.
#' @param output `Character`. The output name. Default: `output.`
#' @param ident_col `Character`. The column with unique sequence IDs. Default: `Sequence_ID`.
#' @param vgene_col `Character`. The column with V gene names. Default: `V_GENE_and_allele`.
#' @param jgene_col `Character`. The column with J gene names. Default: `J_GENE_and_allele`.
#' @param cdr3_col `Character`. The column with CDR3 sequences - nt/aa. Default: `JUNCTION`.
#' @param cell `Character`. Choose T for T cell or B for B cell. Default: `T`.
#' @param output_original `Logical`. if`TRUE`, output and input have same format. Default: `FALSE`.
#' @param mismatch `Numeric.` Percent of mismatches allowed in subgroups. \code{0 == most_stringent}, \code{1 == less_stringent}. Default: `0`
#' @param suffix `Character.` String to be appended to the clonality ID.
#' @param search_genename `Logical.` Accepts IMGT nomenclature and simplify gene IDs. Default: `TRUE`.
#' @param stringdist_method `Character.` Method for distance calculation. See \code{?stringdist-metrics} for more methods.
#' @examples
#' clonality(data = tra)
#' clonality(data = trb)
#' clonality(data = 'example.xlsx')
#' @importFrom stringdist stringdistmatrix
#' @importFrom readxl read_excel
#' @importFrom stringr str_extract
#' @importFrom safejoin safe_full_join
#' @importFrom gtools mixedorder
#' @import dplyr
#' @export

clonality <- function(data = "example.xlsx",
                      output = "output",
                      ident_col = "Sequence_ID",
                      vgene_col = "V_GENE_and_allele",
                      jgene_col = "J_GENE_and_allele",
                      cdr3_col = "JUNCTION",
                      cell = "T",
                      output_original = FALSE,
                      mismatch = 0,
                      suffix = NULL,
                      search_genename = TRUE,
                      stringdist_method = "osa") {


  #Read data.frame/input csv or xlsx table.

  ##ifinput is not a path to an xlsx or csv file, read from R env
  imported_df <- test_input(data)

  #Remove all NA junctions prior to running
  rm.cols <- c(NA,"")
  df <- imported_df[!imported_df[[cdr3_col]] %in% rm.cols, ]

  #Create a simpler data table

  clonal_table <- data.frame(CellId = df[[ident_col]],
                             v_genes = df[[vgene_col]],
                             j_genes = df[[jgene_col]],
                             CDR3 = df[[cdr3_col]],
                             CDR3L = nchar(df[[cdr3_col]]),
                             stringsAsFactors = F)


  if(any(complete.cases(clonal_table) == FALSE)){
    stop("Some rows might have empty V or J genes")
  }

  #Clean IMGT format

  if(search_genename == TRUE) {

    if(cell == "T") {
      clonal_table[["v_genes"]] <- str_extract(string = clonal_table[["v_genes"]], pattern = "TR[AB]V.*?(?=\\*|\\/)")
      clonal_table[["j_genes"]] <- str_extract(string = clonal_table[["j_genes"]], pattern = "TR[AB]J.*?(?=\\*|\\/)")
    }

    if(cell == "B") {
      clonal_table[["v_genes"]] <- str_extract(string = clonal_table[["v_genes"]], pattern = "IG[HK]V.*?(?=\\*|\\/)")
      clonal_table[["j_genes"]] <- str_extract(string = clonal_table[["j_genes"]], pattern = "IG[HK]J.*?(?=\\*|\\/)")
    }

  }

  #First filter, cells must match V gene, J gene and CDR3 Length
  v <- clonal_table[["v_genes"]]
  j <- clonal_table[["j_genes"]]
  l <- clonal_table[["CDR3L"]]

  #Generate a unique ID for each sequence

  id <- paste(v, j, l, sep = "_")

  #Detect identical IDs
  if(any(duplicated(id))){
    duplicated_bool <- unique(id) %in% id[duplicated(id)]
    n_duplicated <- sum(duplicated_bool)
    duplicated_ids <- unique(id)[duplicated_bool]

    freq_matrices <- as.list(c())
    clonal_names <- as.list(c())

  #for each unique sequence:
  #capture the index position of each clonal_table group
  #calculate the distance matrix ratio to the total sequence length
  #get the frequency matrix to freq_matrices
  #get the list of clonal_table names to clonal_names

    for(i in seq(1, n_duplicated)) {
      id_position <- grep(duplicated_ids[i], id)
      dist_mat <- stringdistmatrix(clonal_table[id_position, ][["CDR3"]], useNames = T, method = stringdist_method)
      freq_matrices[[i]] <- dist_mat / l[id_position][1]
      clonal_names[[i]] <- grep(duplicated_ids[i], id, value = T)
    }

  #create a new column for the result
    clonal_table$clonality <- NA
  #order freq_matrices by size
    or <- order(unlist(lapply(lapply(freq_matrices, as.matrix), nrow)), decreasing = T)
    freq_matrices <- freq_matrices[or]
  #order clonal_names
    clonal_names <- clonal_names[or]

  #Subset clones on the CDR3 level accordingly to mismatch parameter.

  #for each frequency matrix,
  #create a data frame with the labels of each sequence
  #cluster number is based on cutree
  #function, 0 = most stringent, 1 = less stringent
  #create an ID for each clonal_table group assign the ID to the original input

    for(i in seq(1, length(freq_matrices))){
      hclust_result <- hclust(freq_matrices[[i]])
      cutree_result <- cutree(hclust_result, h = mismatch)

      df <- data.frame(Seq = labels(cutree_result),
                       Clones = cutree_result,
                       Annotaded_Clone_ID = sprintf("%s.%s", i, cutree_result),
                       Old_Clone_ID = clonal_names[[i]], stringsAsFactors = F)

      clonal_table[id %in% df$Old_Clone_ID, ][["clonality"]] <- df$Annotaded_Clone_ID
    }

  #Add to the unique sequences an unique ID
    unique_bool <- is.na(clonal_table$clonality)
    if(any(unique_bool)){
      n_unique <- nrow(clonal_table[unique_bool, ])
      clonal_table[unique_bool, ]$clonality <- sprintf("U%s", seq(1, n_unique))
    }else {}


  #ifthere is no duplicated sequence id, just call every id unique.

  }else {
  #Make the unique sequences as an unique ID
    clonal_table$clonality <- sprintf("U%s", 1:nrow(clonal_table))
  }


  #Import the output to the original input or save a minimal version
  if(output_original == TRUE) {
    clonal_table <- safe_full_join(x = imported_df, y = clonal_table, by = setNames("CellId", ident_col), conflict = coalesce)
  } else {

    clonal_table <- clonal_table[mixedorder(clonal_table$clonality), ]
  }

  #Add suffix ifthere is one to add.
  if(length(suffix) != 0) {
    clonal_table$clonality <- paste(suffix, clonal_table$clonality, sep = "_")
  }

  #Save output
  if(class(data)[1] != "character") {
    assign(x = output, value = clonal_table, envir = .GlobalEnv)
  }
  if(class(data)[1] == "character") {
    openxlsx::write.xlsx(x = clonal_table, data = sprintf("%s.xlsx", output), row.names = F)
  }



}
