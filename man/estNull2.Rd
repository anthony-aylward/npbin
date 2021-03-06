% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npbin.R
\name{estNull2}
\alias{estNull2}
\title{Null model estimation}
\usage{
estNull2(x, m, mod, prep, init = NULL, iter.max = 200,
  err.max = 1e-06, algorithm = "NLOPT_GN_DIRECT_L",
  pseq = (1:9999)/10000, lb = c(0, 0), ub = rep(log(10000), 2))
}
\arguments{
\item{x}{x}

\item{m}{m}

\item{mod}{mod}

\item{prep}{Info prepared by estNull1}

\item{init}{init}

\item{iter.max}{Maximum number of iterations for EM algorithm}

\item{err.max}{err.max}

\item{algorithm}{Numerical library to use}

\item{pseq}{pseq}

\item{lb}{Lower bounds}

\item{ub}{Upper bounds}
}
\value{
mod
}
\description{
Null model estimation
}
\details{
estimate the null model, part 2
}
\seealso{
\code{\link{dbetabinom.vec}}

Other estNulls: \code{\link{estNull1}},
  \code{\link{estNull}}
}
\concept{estNulls}
