#' A small IMGT TCR Beta sequences dataset.
#'
#' A dataset containing TCR Beta sequences to use as a test for clonality.
#'
#' @format A data frame with 49 rows and 9 variables:
#' \describe{
#'   \item{Sequence number}{Sequence number, not needed}
#'   \item{Sequence ID}{Unique sequence ID, needed}
#'   \item{V-DOMAIN Functionality}{If the V sequence are functional, not needed}
#'   \item{V-GENE and allele}{V gene name, needed}
#'   \item{J-GENE and allele}{J gene name, needed}
#'   \item{D-GENE and allele}{D gene name, not needed}
#'   \item{JUNCTION frame}{If CDR3 junction is in-frame, not needed}
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
#'   \item{Sequence number}{Sequence number, not needed}
#'   \item{Sequence ID}{Unique sequence ID, needed}
#'   \item{V-DOMAIN Functionality}{If the V sequence are functional, not needed}
#'   \item{V-GENE and allele}{V gene name, needed}
#'   \item{J-GENE and allele}{J gene name, needed}
#'   \item{D-GENE and allele}{D gene name, not needed}
#'   \item{JUNCTION frame}{If CDR3 junction is in-frame, not needed}
#'   \item{JUNCTION}{Junction sequence}
#'   \item{sequence}{Full TCR sequence}
#'   ...
#' }
#'
"tra"
