% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cell.classes}
\alias{cell.classes}
\title{Classes to be used}
\format{
A data frame with 49 rows and 9 variables:
\describe{
\item{unique_paired}{List of vectors for paired chain sequences}
\item{sticky_ends}{List of vectors for paired and single chain sequences}
...
}
}
\usage{
cell.classes
}
\description{
List of vectors with TCR or BCR classes.
Classes are sequencing outcomes for each cell, which can be:
no receptor, one chain sequenced or a single chain, both chains sequenced - or a paired chain.
}
\keyword{datasets}
