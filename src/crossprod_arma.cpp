#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
#include "crossprod_arma.h"

// [[Rcpp::export]]
arma::mat self_crossprod_arma(const arma::mat& x){
  return trans(x) * x;
}

// [[Rcpp::export]]
arma::mat crossprod_arma(const arma::mat& x, const arma::mat& y){
  return trans(x) * y;
}

// [[Rcpp::export]]
arma::mat self_tcrossprod_arma(const arma::mat& x){
  return x * trans(x);
}

// [[Rcpp::export]]
arma::mat tcrossprod_arma(const arma::mat& x, const arma::mat& y){
  return x * trans(y);
}