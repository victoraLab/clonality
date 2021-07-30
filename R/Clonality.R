#' Define clonality
#'
#' Returns table with clonality definitions.
#'
#' @param data Data frame object or the full path to a xlsx with rep data.
#' @param output The data output name. Default: output.
#' @param id_col The column with an unique ID. Default: 'Sequence_ID'.
#' @param vgene_col The column with V gene names. Default: 'V_GENE_and_allele'.
#' @param jgene_col The column with J gene names. Default: 'J_GENE_and_allele'.
#' @param cdr3_col The column with CDR3 sequences - nt/aa. Default: 'JUNCTION'.
#' @param cell Choose between T for T cell or B for B cell. Default: 'T'.
#' @param rm.na Remove NA junctions. If `FALSE`, set to unique. Default: 'TRUE'.
#' @param output_original If `TRUE`, output = input with extra clonality column. Default: 'FALSE'.
#' @param mm Percent of mismatches allowed in subgroups. Default: 0
#' @examples
#' clonality(data = tra)
#' clonality(data = trb)
#' clonality(data = 'example.xlsx')
#' @importFrom stringdist stringdistmatrix
#' @importFrom readxl read_excel
#' @importFrom stringr str_extract
#' @importFrom safejoin safe_full_join
#' @importFrom gtools mixedorder
#' @export

clonality <- function(data = "example.xlsx",
                      output = "output",
                      id_col = "Sequence_ID",
                      vgene_col = "V_GENE_and_allele",
                      jgene_col = "J_GENE_and_allele",
                      cdr3_col = "JUNCTION",
                      cell = "T",
                      rm.na = TRUE,
                      output_original = FALSE,
                      mm = 0,
                      suffix = NULL,
                      search_gname = T) {



    # Read data.frame/input xlsx table.

    if (class(data)[1] != "character") {
        df_import <- data
    } else{

        if (any(c("xlsx", "csv")  %in% gsub("^.*\\.", "", data))){
            if (gsub("^.*\\.", "", data) == "xlsx"){
                df_import <- read_excel(data)
            }

            if (gsub("^.*\\.", "", data) == "csv"){
                df_import <- read.csv2(data, sep = ",")
            }
        } else{
            stop("Data can only be csv or xlsx", call. = FALSE)
        }
    }


    # If rm.na is TRUE, remove all NA junctions prior to running

    if (rm.na == TRUE) {
        df <- df_import[!is.na(df_import[[cdr3_col]]), ]
    } else {
        df <- df_import
    }

    # Create a simple table

    clonal_tab <- data.frame(CellId = df[[id_col]],
                             v_genes = df[[vgene_col]],
                             j_genes = df[[jgene_col]],
                             CDR3 = df[[cdr3_col]],
                             CDR3L = nchar(df[[cdr3_col]]),
                             stringsAsFactors = F)

    # Clean IMGT table format

    if (search_gname == T) {

        if (cell == "T") {
            clonal_tab[["v_genes"]] <- str_extract(string = clonal_tab[["v_genes"]], pattern = "TR[AB]V.*?(?=\\*|\\/)")
        }

        if (cell == "B") {
            clonal_tab[["v_genes"]] <- str_extract(string = clonal_tab[["v_genes"]], pattern = "IG[HK]V.*?(?=\\*|\\/)")
        }

        if (cell == "T") {
            clonal_tab[["j_genes"]] <- str_extract(string = clonal_tab[["j_genes"]], pattern = "TR[AB]J.*?(?=\\*|\\/)")
        }

        if (cell == "B") {
            clonal_tab[["j_genes"]] <- str_extract(string = clonal_tab[["j_genes"]], pattern = "IG[HK]J.*?(?=\\*|\\/)")
        }

    }

    v <- clonal_tab[["v_genes"]]
    j <- clonal_tab[["j_genes"]]
    l <- clonal_tab[["CDR3L"]]

    # Generate a unique ID for each sequence

    id <- paste(v, j, l, sep = "_")

    # Detect identical IDs
    if(any(duplicated(id))){

        index_true <- grep(TRUE, id %in% id[duplicated(id)])
        pass <- as.list(c())
        vdl <- as.list(c())

        #for each unique sequence:
        #capture the index position of each clonal_tab group
        #calculate the distance matrix ratio to the total sequence length
        #pass the frequency matrix to pass
        #pass the list of clonal_tab names to vdl

        for (i in seq(1, length(unique(id[index_true]))) ) {

            p <- grep(unique(id[index_true])[i], id)

            dist_mat <- stringdistmatrix(clonal_tab[p, ][["CDR3"]], useNames = T)

            pass[[i]] <- dist_mat / l[p][1]

            vdl[[i]] <- grep(unique(id[index_true])[i], id, value = T)
        }

        # create a new column for the result
        clonal_tab$clonality <- NA
        # order
        or <- order(unlist(lapply(lapply(pass, as.matrix), nrow)), decreasing = T)

        pass <- pass[or]
        # order
        vdl <- vdl[or]

        #for each frequency matrix,
        #create a dataframe with the labels of each sequence
        #cluster number based on the cutree
        # function, 0 = most strigent, 1, less strigent
        #create an ID for each clonal_tab group assign the ID to the original input

        for (i in seq(1, length(pass))) {

            df <- data.frame(Seq = labels(cutree(hclust(pass[[i]]), h = mm)),

                             Clones = cutree(hclust(pass[[i]]), h = mm),

                             ID = sprintf("%s.%s", i, cutree(hclust(pass[[i]]), h = mm)),

                             orig.ident = vdl[[i]], stringsAsFactors = F)

            clonal_tab[id %in% df$orig.ident, ]$clonality <- df$ID
        }

        # Make the unique sequences as an unique ID
        if(all(is.na(clonal_tab$clonality))){}else{
            n <- nrow(clonal_tab[is.na(clonal_tab$clonality), ])
            clonal_tab[is.na(clonal_tab$clonality), ]$clonality <- sprintf("U%s", seq(1, n))
        }



    }

    else{
        # Make the unique sequences as an unique ID
        clonal_tab$clonality <- sprintf("U%s", 1:nrow(clonal_tab))
    }


    # Make the output as the original input or save a minimal version
    if (output_original == TRUE) {
        clonal_tab <- safejoin::safe_full_join(x = df_import, y = clonal_tab, by = setNames("CellId", id_col), conflict = coalesce)
    } else {

        clonal_tab <- clonal_tab[gtools::mixedorder(clonal_tab$clonality), ]
    }

    if (length(suffix) != 0) {
        clonal_tab$clonality <- paste(suffix, clonal_tab$clonality, sep = "_")
    }

    # Save output
    if (class(data)[1] != "character") {
        assign(x = output, value = clonal_tab, envir = .GlobalEnv)
    } else {
        openxlsx::write.xlsx(x = clonal_tab, data = sprintf("%s.xlsx", output), row.names = F)
    }



}
