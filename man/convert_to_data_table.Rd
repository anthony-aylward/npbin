% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{convert_to_data_table}
\alias{convert_to_data_table}
\title{Convert input data frame to a data table}
\usage{
convert_to_data_table(data_frame, minimum_coverage = 5)
}
\arguments{
\item{data_frame}{A data frame containing at least the required fields}

\item{minimum_coverage}{The minimum coverage level that will be enforced}
}
\value{
A data table
}
\description{
Convert input data frame to a data table
}
\details{
A minimum coverage level will be enforced
}
\seealso{
\code{\link{npbin_preprocess_file}}
}
