#include <RcppArmadillo.h>
#include <omp.h>
#include "arma_mat_index.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
arma::mat row_arma(const arma::mat& x, int cores){
  unsigned int i;
  
  arma::mat m(x.n_rows, x.n_cols);
  
  #pragma omp parallel for
  for(i = 1; i <= x.n_rows; i++){
    m.row(i-1).fill(i);
  }
  
  return m;
}

// [[Rcpp::export]]
arma::mat col_arma(const arma::mat& x, int cores) { 
  unsigned int i;
  arma::mat m(x.n_rows, x.n_cols);
  omp_set_num_threads(cores);
  
  #pragma omp parallel for
  for(i = 1; i <= x.n_rows; i++){
      m.col(i-1).fill(i);
  }
  
  return m;
}