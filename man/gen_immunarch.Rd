% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_immunarch.R
\name{gen_immunarch}
\alias{gen_immunarch}
\title{Generates an immunarch compatible input}
\usage{
gen_immunarch(
  seu = seu,
  proportion = "rm.na",
  sample_column = NULL,
  metadata_columns = NULL
)
}
\arguments{
\item{seu}{Seurat with clonality generated columns.}

\item{proportion}{Character. Possible values: \code{all} against all cells, \code{rm.na} remove NAs prior to proportion.}

\item{sample_column}{Character. Column to be used as sample to cut the data.}

\item{metadata_columns}{Character. Columns to split the data.}
}
\description{
Returns a dataset compatible with the immunarch package.
}
\examples{
gen_immunarch(seu, proportion = "rm.na",  sample_column = "Conditions")
}
