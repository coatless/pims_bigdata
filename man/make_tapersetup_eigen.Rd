% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{make_tapersetup_eigen}
\alias{make_tapersetup_eigen}
\title{Create a taper setup for an Eigen sparse matrix}
\usage{
make_tapersetup_eigen(d, taprange)
}
\arguments{
\item{d}{A \code{matrix} that is a symmetric distance matrix.}

\item{taprange}{A \code{double} that is used to specify the gamma.}
}
\value{
A \code{list}
\itemize{
\item \code{n} - row/col of matrix
\item \code{good.dists} - distances that meet tapering requirements
\item \code{taps} - results of applying the taper function
\item \code{ja} - column index
\item \code{ia} - row index
}
}
\description{
Generate the storage information we'll need to create a sparse eigen matrix
}
\details{
This function is meant to be used to create SPAM matrices.
The function formats data so that it can be used with a compressed sparse column format.
Do not use output for Eigen matrices.
}
\examples{
data(anom1962)
d = rdist_earth1(head(loc))
make_tapersetup_eigen(d, 20)
}

