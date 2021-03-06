% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{npbin_preprocess_counts}
\alias{npbin_preprocess_counts}
\title{Preprocess allele count data}
\usage{
npbin_preprocess_counts(data_frame, chr_index = 1, location_index = 2,
  total_index = 3, ref_count_index = 4)
}
\arguments{
\item{data_frame}{A data frame containing at least the four required fields}

\item{chr_index}{Index of the chromosome field}

\item{chr_index}{Index of the location field}

\item{chr_index}{Index of the total read count field}

\item{chr_index}{Index of the reference allele count field}
}
\value{
A data frame
}
\description{
Perform preprocessing steps to prepare data for NPBin analysis
}
\details{
This function extracts the required fields from the input data and returns
a data frame formatted for NPBin analysis.
}
\seealso{
\code{\link{logical_to_winning_chip}},
  \code{\link{npbin_preprocess_file}}
}
