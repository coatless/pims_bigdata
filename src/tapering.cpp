#include <RcppArmadillo.h>
#include <omp.h>
#include "arma_mat_index.h"
using namespace Rcpp;

// -----------------------------------------------------
// -----------------------------------------------------
//
// Tapering Functions
//
// -----------------------------------------------------
// -----------------------------------------------------

//' @title Taper Function
//' @description Taper function we're using
//' @param d A \code{vec} containing the small distance values
//' @param taprange A \code{double} indicating the taper range
//' @return A \code{vec} that has had tapering applied to d.
//' @examples
//' data(anom1962)
//' d = rdist_earth1(head(loc))
//' wendland2_1(d, 20)
// [[Rcpp::export]]
arma::vec wendland2_1(arma::vec d, double taprange){
  d = d/taprange;
  
  return pow(1.0-d,4)%(4*d+1)%(d<1);
}

//' @title Create a taper setup for a SPAM matrix
//' @description Get the storage information we'll need to create a sparse matrix
//' @param d A \code{matrix} that is a symmetric distance matrix.
//' @param taprange A \code{double} that is used to specify the gamma.
//' @return A \code{list}
//' \itemize{
//' \item \code{n} - row/col of matrix
//' \item \code{good.dists} - distances that meet tapering requirements
//' \item \code{taps} - results of applying the taper function
//' \item \code{ja} - ordered column indices
//' \item \code{ia} -  pointer to the beginning of each row in the arrays entries and colindices
//' \item \code{index} - element position in matrix 
//' }
//' @details This function is meant to be used to create SPAM matrices.
//' Do not use output for Eigen matrices.
//' @examples
//' data(anom1962)
//' d = rdist_earth1(head(loc))
//' make_tapersetup_R(d, 20)
// [[Rcpp::export]]
Rcpp::List make_tapersetup_R(arma::mat d, double taprange){
  unsigned int n = d.n_rows;
  
  // Find elements that are less than (logical matrix e.g. 0/1's)
  arma::umat inrange = d < taprange;
  
  // Same idea, but return element ids.
  arma::uvec inrange_id = find(inrange);
  
  // Convert matrix to allow operations
  arma::mat temp = arma::conv_to<arma::mat>::from(inrange);

  // Obtain the small distance elements
  arma::vec good_dists = d.elem(inrange_id);
  
  // Apply taper function (this cannot be made dynamic)
  arma::vec taps = wendland2_1(good_dists, taprange);
  
  // Create ordered column indices of the nonzero values
  arma::vec ja = row_arma(d).elem(inrange_id);
  
  // Calculate pointer to the beginning of each row in the arrays' entries and colindices
  arma::vec ia(inrange.n_rows + 1);
  ia(0) = 1;
  ia.rows(1,inrange.n_rows) = 1 + arma::cumsum(arma::sum(temp,1));
  
  // Create a column index matrix and isolate ideal elements
  arma::vec index = (col_arma(d).elem(inrange_id) - 1) * n + ja;
  
  // Export results as an R list object
  return Rcpp::List::create(
                          Rcpp::Named("n") = n,
                          Rcpp::Named("good_dists") = good_dists,
                          Rcpp::Named("taps") = taps,
                          Rcpp::Named("ja") = ja,
                          Rcpp::Named("ia") = ia,
                          Rcpp::Named("index") = index);
}


//' @title Create a taper setup for an Eigen sparse matrix
//' @description Generate the storage information we'll need to create a sparse eigen matrix
//' @param d A \code{matrix} that is a symmetric distance matrix.
//' @param taprange A \code{double} that is used to specify the gamma.
//' @return A \code{list}
//' \itemize{
//' \item \code{n} - row/col of matrix
//' \item \code{good.dists} - distances that meet tapering requirements
//' \item \code{taps} - results of applying the taper function
//' \item \code{ja} - column index
//' \item \code{ia} - row index
//' }
//' @details This function is meant to be used to create SPAM matrices.
//' The function formats data so that it can be used with a compressed sparse column format.
//' Do not use output for Eigen matrices.
//' @examples
//' data(anom1962)
//' d = rdist_earth1(head(loc))
//' make_tapersetup_eigen(d, 20)
// [[Rcpp::export]]
Rcpp::List make_tapersetup_eigen(arma::mat d, double taprange){
  double n = d.n_rows;
  
  // Find elements that are less than taprange and return element ids.
  arma::uvec inrange_id = find(d < taprange);
  
  // Approximate how many zero entries exist per columnss
  int rescol = floor(double(inrange_id.n_elem) / n); 

  // Obtain the small distance elements
  arma::vec good_dists = d.elem(inrange_id);
  
  // Apply taper function
  arma::vec taps = wendland2_1(good_dists, taprange);
  
  // Create a row index matrix and extract out ranges
  arma::vec ja = row_arma(d).elem(inrange_id) - 1;
  
  // Create a col index matrix and extract the ranges
  arma::vec ia = col_arma(d).elem(inrange_id) - 1;
    
  // Export results as an R list object
  return Rcpp::List::create(
                          Rcpp::Named("n") = n,
                          Rcpp::Named("good.dists") = good_dists,
                          Rcpp::Named("taps") = taps,
                          Rcpp::Named("ja") = ja,
                          Rcpp::Named("ia") = ia,
                          Rcpp::Named("rescol") = rescol);
}


// -----------------------------------------------------
// -----------------------------------------------------
//
// Likelihood functions
//
// -----------------------------------------------------
// -----------------------------------------------------


//' @title Create a taper setup for an Eigen sparse matrix
//' @description Generate the storage information we'll need to create a sparse eigen matrix
//' @param x A \code{double} that represents rho, which is the scaling parameter of the exponential covariance function
//' @param d A \code{matrix} that is a symmetric distance matrix.
//' @param z A \code{vec} that contains the finite number of observations
//' @param n A \code{double} that indicates the number of rows/cols of the matrix.
//' @return A \code{double}
//' @examples
//' data(anom1962)
//' d = rdist_earth1(head(loc))
//' nll_arma(20, d, z, n)
// [[Rcpp::export]]
double nll_arma(double x,
                arma::mat d,
                arma::vec z,
                double n){
  
  arma::mat corr_matrix = exp(-d/x);
  arma::mat Q = arma::chol(corr_matrix);
  double logdet = arma::as_scalar( 2.0 * sum(log(arma::diagvec(Q))) );
  arma::mat bs = arma::solve( trimatl( trans(Q) ), z);
  double distval = arma::as_scalar(trans(bs)*bs);
  return (n * log(2.0 * arma::datum::pi) + n * log(distval / n) + logdet + n)/2.0;
}

//' @describeIn nll_arma
//' @param cores A \code{integer} representing how many cores should be loaded
// [[Rcpp::export]]
arma::vec nll_parallel(arma::vec x, arma::mat d, arma::vec z, double n, int cores=1){
  unsigned int i;
  
  #pragma omp parallel for num_threads(cores) private(i)
  for(i = 0; i < x.n_elem; i++){
    x(i) = nll_arma(x(i), d, z, n);
  }
  
  return x;
}

/*
// arma lacks sparse operations for chol, det
//
// To see an implementation of the taper functions
// See the eigen_nll file.
//
// When armadillo does implement the sparse options,
// It'll probably look like the following
// [[Rcpp::export]]
double nll_1taper(double x, double n, 
                  const arma::vec& good_dists, 
                  const arma::vec& taps,
                  const arma::vec& ia,
                  const arma::vec& ja,
                  const arma::vec& z){

  arma::sp_mat corr_matrix_Tape(
                                join_rows(ia, ja), // row, col location
                                exp(-good_dists/x) % taps, //values
                                n, n // dimension of matrix
                               )
  
  // Hypothetical... Pending future armadillo update!
  // Find the Cholesky decomposition
  arma::mat Q = chol(corr_matrix_Taper);
  
  // Take the deterimant
  double logdet = 2 * log(det(Q));
  
  // Solve system
  arma::vec bs = solve(Q, z);
  
  // Take inner product
  double distval = dot(bs,bs);

  // Return y
  return (n * log(2 * pi) + n * log(distval / n) + logdet + n)/2;
}

// [[Rcpp::export]]
double nll_2taper(double x, double n, 
                  const arma::vec& good_dists, 
                  const arma::vec& taps,
                  const arma::vec& ia,
                  const arma::vec& ja,
                  const arma::vec& z)
{
  arma::sp_mat corr_matrix_Tape = arma::sp_mat(
                                join_rows(ia, ja), // row, col location
                                exp(-good_dists/x) % taps, //values
                                n, n // dimension of matrix
                               );
                          
  arma::mat Q = chol(corr_matrix_Taper);
  
  double logdet = 2 * log(det(Q));
  
  arma::mat inv = backsolve_arma(Q, 
                          forwardsolve_arma(Q, eye(n))
                          );
  
                               
  for(unsigned int k = 0; k < good_dists.n_elem; k++){
    unsigned int i = ia(k), j = ja(k);
    corr_matrix_Taper(i,j) = inv(i,j) * taps(k);
  }
  
  double distval = as_scalar(trans(z) * corr_matrix_Taper * z);
  return  (n * log(2.0 * M_PI) + n * log(distval / n) + logdet + n)/2.0;
}*/
