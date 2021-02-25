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
