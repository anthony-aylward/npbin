% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npbin.R
\name{estNull1}
\alias{estNull1}
\title{Null model estimation}
\usage{
estNull1(mod, pseq = (1:9999)/10000, ncores = 1)
}
\arguments{
\item{mod}{mod}

\item{pseq}{pseq}

\item{ncores}{Number of cores to use}
}
\value{
List describing the model
}
\description{
Null model estimation
}
\details{
estimate the null model, part 1
}
\seealso{
\code{\link{evBspl}}, \code{\link{evBsplDrv}}

Other estNulls: \code{\link{estNull2}},
  \code{\link{estNull}}
}
\concept{estNulls}
