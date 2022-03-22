#' A small IMGT TCR Beta sequences dataset.
#'
#' A dataset containing TCR Beta sequences to use as a test for clonality.
#'
#' @format A data frame with 49 rows and 9 variables:
#' \describe{
#'   \item{Sequence_number}{Sequence number, not needed}
#'   \item{Sequence_ID}{Unique sequence ID, needed}
#'   \item{V_DOMAIN_Functionality}{If the V sequence are functional, not needed}
#'   \item{V_GENE_and_allele}{V gene name, needed}
#'   \item{J_GENE_and_allele}{J gene name, needed}
#'   \item{D_GENE_and_allele}{D gene name, not needed}
#'   \item{JUNCTION_frame}{If CDR3 junction is in-frame, not needed}
#'   \item{JUNCTION}{Junction sequence}
#'   \item{sequence}{Full TCR sequence}
#'   ...
#' }
#'
"trb"

#' A small IMGT TCR Alpha sequences dataset.
#'
#' A dataset containing TCR Alpha sequences to use as a test for clonality.
#'
#' @format A data frame with 49 rows and 9 variables:
#' \describe{
#'   \item{Sequence_number}{Sequence number, not needed}
#'   \item{Sequence_ID}{Unique sequence ID, needed}
#'   \item{V_DOMAIN_Functionality}{If the V sequence are functional, not needed}
#'   \item{V_GENE_and_allele}{V gene name, needed}
#'   \item{J_GENE_and_allele}{J gene name, needed}
#'   \item{D_GENE_and_allele}{D gene name, not needed}
#'   \item{JUNCTION_frame}{If CDR3 junction is in-frame, not needed}
#'   \item{JUNCTION}{Junction sequence}
#'   \item{sequence}{Full TCR sequence}
#'   ...
#' }
#'
"tra"

#' A 10x genomics TCR sequencing dataset.
#'
#' A dataset containing TCR Alpha and Beta sequences from 10x genomics-generated single cells.
#'
#' @format A data frame with 49 rows and 9 variables:
#' \describe{
#'   \item{barcode}{Unique barcode sequences to indentify single cells}
#'   \item{contig_id}{Unique ID for contig assemble per single cell}
#'   \item{length}{chain length in bp}
#'   \item{chain}{The TCR chain name}
#'   \item{v_gene}{The V-gene name}
#'   \item{d_gene}{The D-gene name}
#'   \item{j_gene}{The J-gene name}
#'   \item{c_gene}{The constant gene name}
#'   \item{productive}{If the TCR sequence is productive}
#'   \item{cdr3}{The amino acid sequence of the CDR3}
#'   \item{cdr3_nt}{The nucleotide sequence of the CDR3}
#'   \item{reads}{The read sequence number for each barcode}
#'   \item{umis}{The number of umis for each cell}
#'   \item{raw_clonotype_id}{The clonotype id designated to each cell by the 10x genomics protocol}
#'   ...
#' }
#'
#'
"filtered_contig_annotations"
