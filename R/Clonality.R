#' Define clonality
#'
#' Returns the input with an additional clonality column with clonal definitions.
#'
#' @param File Path to a xlsx input containing the clonality data to be defined.
#' @param NewF The output name.
#' @param Vgene_Column The column name containing V gene names.
#' @param Jgene_Column The column name containing J gene names.
#' @param CDR3_Column The column name containing CDR3 sequences - nt or aa.
#' @param Mismatch Percent of mismatches allowed in CDR3 before subsetting a group.
#' @param rm.junc.na Remove rows with NA junctions. If `FALSE`, NA rows are considered unique.
#' @param cell Choose between T for T cell or B for Bcell.
#' @importFrom stringdist stringdistmatrix
#' @importFrom readxl read_excel
#' @importFrom stringr str_extract
#' @importFrom dplyr full_join
#' @export


Clonality <- function (File = "t.test.xlsx",
                       NewF = "output.xlsx",
                       ID_Column = "Sequence ID",
                       Vgene_Column = "V-GENE and allele",
                       Jgene_Column = "J-GENE and allele",
                       CDR3_Column = "JUNCTION",
                       cell = "T",
                       rm.junc.na = TRUE,
                       output.orig = FALSE,
                       Mismatch = 0)


{


#Read input xlsx table.
  df.import <- read_excel(File)

  if(rm.junc.na == TRUE){
    df <- df.import[!is.na(df.import[[CDR3_Column]]),]
  }

  #Create simple table
  clonal <- data.frame(CellId = df[[ID_Column]],
                       Vgene = df[[Vgene_Column]],
                       Jgene = df[[Jgene_Column]],
                       CDR3 = df[[CDR3_Column]],
                       CDR3L = nchar(df[[CDR3_Column]]), stringsAsFactors = F)
  #Clean IMGT format

  if(cell == "T"){
    clonal[["Vgene"]] <- str_extract(string = clonal[["Vgene"]], pattern = "(TR[AB]V[0-9]{1,3}[DN]{0,1}[-]{0,1}[0-9]{0,2})")
  }

  if(cell == "B"){
    clonal[["Vgene"]] <- str_extract(string = clonal[["Vgene"]], pattern = "(IG[HK]V[0-9]{1,3}[DN]{0,1}[-]{0,1}[0-9]{0,3})")
  }

  if(cell == "T"){
    clonal[["Jgene"]] <- str_extract(string = clonal[["Jgene"]], pattern = "(TR[AB]J[0-9]{1,3}[DN]{0,1}[-]{0,1}[0-9]{0,2})")
  }

  if(cell == "B"){
    clonal[["Jgene"]] <- str_extract(string = clonal[["Jgene"]], pattern = "(IG[HK]J[0-9]{1,3}[DN]{0,1}[-]{0,1}[0-9]{0,3})")
  }


  V.genes <- clonal[["Vgene"]]
  J.genes <- clonal[["Jgene"]]
  cdr3    <- clonal[["CDR3"]]
  l.cdr3  <- clonal[["CDR3L"]]

#Creates unique ID for each cell
  comp <-  paste(V.genes, J.genes, l.cdr3, sep = "_")

  index.true <- grep(TRUE, comp %in% comp[duplicated(comp)])
  pass <- as.list(c())
  V_J_L <- as.list(c())

  for (i in seq(1, length(unique(comp[index.true])))) {
    p <- grep(unique(comp[index.true])[i], comp)
    dist <- stringdistmatrix(clonal[p, ][["CDR3"]], useNames = T)
    dist.mat <- dist/l.cdr3[p][1]
    pass[[i]] <- dist.mat
    V_J_L[[i]] <- grep(unique(comp[index.true])[i], comp, value = T)
  }

  clonal$Clonality <- NA

  for (i in seq(1, length(pass))) {
    df <- data.frame(Seq = labels(cutree(hclust(pass[[i]]), h = Mismatch)),
                     Clones = cutree(hclust(pass[[i]]), h = Mismatch),
                     ID = sprintf("%s.%s", i, cutree(hclust(pass[[i]]),h = Mismatch)), stringsAsFactors = F)

    df2 <- clonal[clonal[["CDR3"]] %in% df$Seq, ]
    V.genes.2 <- df2[["Vgene"]]
    J.genes.2 <- df2[["Jgene"]]
    l.cdr.2 <- df2[["CDR3L"]]

    comp.2 <- paste(V.genes.2, J.genes.2, l.cdr.2, sep = "_")


    df2$Clonality[comp.2 %in% V_J_L[[i]]] <- df$ID
    clonal[clonal[["CDR3"]] %in% df$Seq, ]$Clonality <- df2$Clonality
  }

  n <- nrow(clonal[is.na(clonal$Clonality), ])

  clonal[is.na(clonal$Clonality), ]$Clonality <- sprintf("U%s", seq(1, n))



  if(output.orig == TRUE){
    clonal  <- full_join(df.import, clonal, by = setNames("CellId",ID_Column))
}

  openxlsx::write.xlsx(x = clonal, file = NewF, row.names = F)
}
