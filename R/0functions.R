
#' @title Taper Function
#' @description Taper function we're using
#' @param d A \code{vec} containing the small distance values
#' @param taprange A \code{double} indicating the taper range
#' @return A \code{vec} that has had tapering applied to d.
#' @examples
#' data(anom1962)
#' d = rdist.earth(head(loc))
#' wendland2.1(d, 20)
wendland2.1 <- function(d, taprange){
  d <- d/taprange
  return((1-d)^4*(4*d+1)*(d<1))
}

## Get the storage information we'll need to create a sparse matrix

#' @title Create a taper setup for a SPAM matrix
#' @description Get the storage information we'll need to create a sparse matrix
#' @param d A \code{matrix} that is a symmetric distance matrix.
#' @param f.tap A \code{function} that is the tapering function to use.
#' @param taprange A \code{double} that is used to specify the gamma.
#' @return A \code{list}
#' \itemize{
#' \item \code{n} - row/col of matrix
#' \item \code{good.dists} - distances that meet tapering requirements
#' \item \code{taps} - results of applying the taper function
#' \item \code{ja} - ordered column indices
#' \item \code{ia} -  pointer to the beginning of each row in the arrays entries and colindices
#' \item \code{index} - element position in matrix 
#' }
#' @details This function is meant to be used to create SPAM matrices.
#' Do not use output for Eigen matrices.
#' @examples
#' data(anom1962)
#' d = rdist.earth(head(loc))
#' make.tapersetup(d, wendland2.1, 20)
make.tapersetup <- function(d, f.tap, taprange)
{
  n <- nrow(d)
  inrange <- d < taprange
  good.dists <- d[inrange]
  taps <- f.tap(good.dists, taprange)
  names(taps) <- names(good.dists) <- NULL
  ja <- row(d)[inrange]
  ia <- as.integer(c(1, 1 + cumsum(apply(inrange, 1, sum))))
  index <- (col(d)[inrange] - 1) * n + ja
  return(list(n = n, good.dists = good.dists, taps = taps,
              ja = ja, ia = ia, index = index))
}

## Likelihood functions -- these all have timing built into them
## For general use, you may want to remove the timings

#' @title Create a taper setup for an Eigen sparse matrix
#' @description Generate the storage information we'll need to create a sparse eigen matrix
#' @param x A \code{double} that represents rho, which is the scaling parameter of the exponential covariance function
#' @return A \code{double}
#' @examples
#' data(anom1962)
#' d = rdist_earth1(head(loc))
#' nll_arma(20)
nll <- function(x){
  print(x)
  times <- rep(NA, 5)
  times[1] <- system.time(corr.matrix <- exp(-d/x))[1]
  times[2] <- system.time(Q <- chol(corr.matrix))[1]
  times[3] <- system.time(logdet <- 2 * sum(log(diag(Q))))[1]
  times[4] <- system.time(bs <- backsolve(Q, z, transpose = TRUE))[1]
  times[5] <- system.time(distval <- drop(crossprod(bs)))[1]
  y <- (n * log(2 * pi) + n * log(distval / n) + logdet + n)/2
  attr(y, "times") <- times
  return(y)
}

#' @title Create a taper setup for an Eigen sparse matrix
#' @description Generate the storage information we'll need to create a sparse eigen matrix
#' @param x A \code{double} that represents rho, which is the scaling parameter of the exponential covariance function
#' @param d A \code{matrix} that is a symmetric distance matrix.
#' @param z A \code{vec} that contains the finite number of observations
#' @param n A \code{double} that indicates the number of rows/cols of the matrix.
#' @return A \code{double}
#' @examples
#' data(anom1962)
#' d = rdist_earth1(head(loc))
#' nll_arma(20, d, z, n)
nll_closed <- function(x,d,z,n){
  print(x)
  times <- rep(NA, 5)
  times[1] <- system.time(corr.matrix <- exp(-d/x))[1]
  times[2] <- system.time(Q <- chol(corr.matrix))[1]
  times[3] <- system.time(logdet <- 2 * sum(log(diag(Q))))[1]
  times[4] <- system.time(bs <- backsolve(Q, z, transpose = TRUE))[1]
  times[5] <- system.time(distval <- drop(crossprod(bs)))[1]
  y <- (n * log(2 * pi) + n * log(distval / n) + logdet + n)/2
  attr(y, "times") <- times
  return(y)
}

## The next two functions require the spam package

#' @title Second Tapering Implementation
#' @description The function implements the second tapering method described
#' in the Covariance Tapering paper.
#' @param x A \code{double} indicating the rho, or scaling parameter
#' @param setup Results returned from make.tapersetup 
#' @param z A \code{vec} that contains the finite number of observations
#' @examples
#' data(anom1962)
#' d = rdist_earth1(loc)
#' setup.eigen = make.tapersetup(d,taprange = 50)
#' val = nll.1taper(20, setup, z)
nll.1taper <- function(x, setup, z){
  print(x)
  times <- rep(NA, 5)
  times[1] <- system.time(corr.matrix.Taper <- new("spam",
                                                   entries = exp(-setup$good.dists/x) *
                                                   setup$taps,
                                                   colindices = setup$ja,
                                                   rowpointers = setup$ia,
                                                   dimension = as.integer(rep(n, 2))))[1]
  times[2] <- system.time(Q <- chol.spam(corr.matrix.Taper))[1]
  times[3] <- system.time(logdet <- 2 * as.numeric(determinant(Q, log = TRUE)$modulus))[1]
  times[4] <- system.time(bs <- forwardsolve.spam(Q, z))[1] # Same as usual backsolve(, transpose=TRUE)!
  times[5] <- system.time(distval <- drop(crossprod(bs)))[1]
  y <- (setup$n * log(2 * pi) + setup$n * log(distval / setup$n) + logdet + setup$n)/2
  attr(y, "times") <- times
  return(y)
}

#' @title Second Tapering Implementation
#' @description The function implements the second tapering method described
#' in the Covariance Tapering paper.
#' @param x A \code{double} indicating the rho, or scaling parameter
#' @param setup Results returned from make.tapersetup 
#' @param z A \code{vec} that contains the finite number of observations
#' @param rescol An \code{int} that contains how many non-zero elements should be in a column. 
#' @return A \code{double}
#' @examples
#' data(anom1962)
#' d = rdist_earth1(loc)
#' setup.eigen = make_tapersetup_eigen(d,taprange = 50)
#' val = nll_2taper(20, setup.eigen$n,
nll.2taper <- function(x, setup, z){
  print(x)
  times <- rep(NA, 6)
  times[1] <- system.time(corr.matrix.Taper <- new("spam",
                                                   entries = exp(-setup$good.dists/x) * setup$taps,
                                                   colindices = setup$ja,
                                                   rowpointers = setup$ia,
                                                   dimension = as.integer(rep(n, 2))))[1]
  times[2] <- system.time(Q <- chol(corr.matrix.Taper))[1]
  times[3] <- system.time(logdet <- 2 * as.numeric(determinant(Q, log = TRUE)$modulus))[1]
  times[4] <- system.time(inv <- backsolve(Q, forwardsolve(Q, diag(n))))[1]
  times[5] <- system.time(slot(corr.matrix.Taper, "entries") <- inv[setup$index] *
                          setup$taps)[1]
  times[6] <- system.time(distval <- drop(t(z) %*% corr.matrix.Taper %*% z))[1]
  y <- (n * log(2 * pi) + n * log(distval / n) + logdet + n)/2
  attr(y, "times") <- times
  return(y)
}
