#include <RcppArmadillo.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
arma::mat nll_arma_corr(double x, arma::mat d){
  return exp(-d/x);
}

// [[Rcpp::export]]
arma::mat nll_arma_chol(arma::mat corr_matrix){
  return arma::chol(corr_matrix);
}

// [[Rcpp::export]]
double nll_arma_det(arma::mat Q){
  return arma::as_scalar( 2.0 * sum(log(arma::diagvec(Q))) );
}

// [[Rcpp::export]]
arma::mat nll_arma_bs(arma::mat Q, arma::vec z){
  return arma::solve( trimatl( trans(Q) ), z);
} 

// [[Rcpp::export]]
double nll_arma_dist(arma::mat bs){
  return arma::as_scalar(trans(bs)*bs);
}

// [[Rcpp::export]]
double nll_arma_all(double n, double distval, double logdet){
  return (n * log(2.0 * arma::datum::pi) + n * log(distval / n) + logdet + n)/2.0;
}
