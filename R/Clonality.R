#' Define clonality
#'
#' Returns the input with an additional clonality column with clonal definitions.
#'
#' @param File Data frame or full path to a xlsx input containing TCR/BCR repertoire data.
#' @param NewF The file output name.
#' @param ID_Column The column containing an individual sequence ID.
#' @param Vgene_Column The column name containing V gene names.
#' @param Jgene_Column The column name containing J gene names.
#' @param CDR3_Column The column name containing CDR3 sequences - nt or aa.
#' @param cell Choose between T for T cell or B for Bcell.
#' @param rm.junc.na Remove rows with NA junctions. If `FALSE`, NA rows are considered unique sequences.
#' @param output.orig If `TRUE`, outputs the original input, if `FALSE` outputs a minimum table.
#' @param Mismatch Percent of mismatches allowed in CDR3 before subsetting a group.
#' @examples
#' clonality(File = tra)
#' clonality(File = trb)
#' clonality(File = "example.xlsx")
#' @importFrom stringdist stringdistmatrix
#' @importFrom readxl read_excel
#' @importFrom stringr str_extract
#' @importFrom dplyr full_join
#' @importFrom gtools mixedorder
#' @export

clonality <- function (File = "example.xlsx",
                       NewF = "output",
                       ID_Column = "Sequence ID",
                       Vgene_Column = "V-GENE and allele",
                       Jgene_Column = "J-GENE and allele",
                       CDR3_Column = "JUNCTION",
                       cell = "T",
                       rm.junc.na = TRUE,
                       output.orig = FALSE,
                       Mismatch = 0)


{


#Read data.frame/input xlsx table.
  if(class(File)[1] != "character"){
    df.import <- File
  }
  else{
    df.import <- read_excel(File)
    }

# if rm.junc.na is TRUE, remove all NA junctions prior to running

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

#Detect identical sequences
  index.true <- grep(TRUE, comp %in% comp[duplicated(comp)])
  pass <- as.list(c())
  V_J_L <- as.list(c())

  for (i in seq(1, length(unique(comp[index.true])))) { #for each unique sequence
    p <- grep(unique(comp[index.true])[i], comp)  #capture the index position of each clonal group
    dist <- stringdistmatrix(clonal[p, ][["CDR3"]], useNames = T) #calculate the distance matrix
    dist.mat <- dist/l.cdr3[p][1] #get the distance ratio to the total sequence length
    pass[[i]] <- dist.mat #pass the frequency matrix to pass
    V_J_L[[i]] <- grep(unique(comp[index.true])[i], comp, value = T) #pass the list of clonal names to V_J_L
  }

  clonal$Clonality <- NA #create a new column for clonal

  pass <- pass[order(unlist(lapply(lapply(pass, as.matrix), nrow)))] #order
  V_J_L <- V_J_L[order(unlist(lapply(V_J_L, length)))] #order

  for (i in seq(1, length(pass))) {  #for each frequency matrix, create a dataframe with
    df <- data.frame(Seq = labels(cutree(hclust(pass[[i]]), h = Mismatch)), #the labels of each sequence
                     Clones = cutree(hclust(pass[[i]]), h = Mismatch), #a cluster number based on the cutree function, 0 = most strigent, 1, less strigent
                     ID = sprintf("%s.%s", i, cutree(hclust(pass[[i]]),h = Mismatch)),
                     orig.ident = V_J_L[[i]], stringsAsFactors = F) #create an ID for each clonal group


    clonal[comp %in% df$orig.ident, ]$Clonality <- df$ID #assign the ID to the original input
  }
#Make the unique sequences as an unique ID
  n <- nrow(clonal[is.na(clonal$Clonality), ])
  clonal[is.na(clonal$Clonality), ]$Clonality <- sprintf("U%s", seq(1, n))

  #Make the output as the original input or save a minimal version
  if(output.orig == TRUE){
    clonal  <- full_join(df.import, clonal, by = setNames("CellId", ID_Column))
  }else{
    clonal <- clonal[gtools::mixedorder(clonal$Clonality),]
  }


  #Save output
  if(class(File)[1] != "character"){
    assign(x = NewF, value = clonal, envir = .GlobalEnv)
  }
  else{
    openxlsx::write.xlsx(x = clonal, file = sprintf("%s.xlsx", NewF), row.names = F)
  }


}
