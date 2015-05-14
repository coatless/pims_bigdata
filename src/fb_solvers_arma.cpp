#include <RcppArmadillo.h>
#include "fb_solvers_arma.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


// Solves R * x = b
// [[Rcpp::export]]
arma::mat backsolve_arma(arma::mat r, const arma::vec& x, 
                    unsigned int k = 0,
                    bool transpose = false){
    //submat first_row, first_col, last_row, last_col
    if(k != 0 && k != r.n_rows){
      r = r.submat(0,0, k, k);
    }
    
    if(!transpose){
      r = solve(trimatu(r), x); 
    }else{
      r = solve(trans(r),x);
    }
    return r;
}
             

// solve L * x = b
// [[Rcpp::export]]
arma::mat forwardsolve_arma(arma::mat l, const arma::vec& x, 
                       unsigned int k = 0,
                       bool transpose = false){
                         
    //submat first_row, first_col, last_row, last_col
    if(k != 0 && k != l.n_rows){
      l = l.submat(0, 0, k, k);
    }
    
    if(!transpose){
      l = solve(trimatl(l), x); 
    }else{
      l = solve(trans(l),x);
    }
    
    return l;
}

// ---------------------------------------------------
