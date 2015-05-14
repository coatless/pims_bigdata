#include <RcppArmadillo.h>
#include "rdist_earth_arma.h"
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]

bool sign(const double& x){
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

/*// [[Rcpp::export]]
mat ifelse_arma(mat a) {
  unsigned int i,j;
  
  for (i=0; i < a.n_rows; i++) {
    for(j=0; j < a.n_cols; j++){
      double pt = a(i, j);
      if(abs(pt) > 1){
        a(i,j) = (pt > 0) - (pt <0);
      }
    }
  }
  return a; 
}*/

// [[Rcpp::export]]
arma::mat ifelse_arma(arma::mat a){
  arma::uvec scan = find(abs(a) > 1);
  arma::vec out = a.elem(scan);
  for(unsigned int i = 0; i < out.n_elem; i++){
    double pt = out(i);
    out(i) = (pt > 0) - (pt < 0); 
  }
  a.elem(scan) = out;
  return a;
}


// [[Rcpp::export]]
arma::mat rdist_earth1(arma::mat x1, 
                       bool miles = true, double R = 0) 
{
    if (R == 0) { // null check
        if (miles) { R = 3963.34; }
        else{ R = 6378.388; }
    }
    arma::vec coslat1 = cos((x1.col(1) * arma::datum::pi)/180.0);

    
    arma::mat pp(x1.n_rows,3);
    
    // coslat1 % coslon1
    pp.col(0) = coslat1 % cos((x1.col(0) * arma::datum::pi)/180.0);
    // coslat1 % sinlon1
    pp.col(1) = coslat1 % sin((x1.col(0) * arma::datum::pi)/180.0);
    // sinlat1
    pp.col(2) = sin((x1.col(1) * arma::datum::pi)/180.0);
    
    pp = pp * trans(pp);
    return R * acos(ifelse_arma(pp));
    //return R * acos( ifelse_arma(pp) );
}

arma::mat rdist_earth2(arma::mat x1, arma::mat x2,
                       bool miles = true, double R = 0) 
{
    if (R == 0) { // null check
        if (miles) { R = 3963.34; }
        else{ R = 6378.388; }
    }
    arma::vec coslat1 = cos((x1.col(1) * arma::datum::pi)/180.0);
    
    arma::mat pp(x1.n_rows,3);
    
    pp.col(0) = coslat1 % cos((x1.col(0) * arma::datum::pi)/180.0);
    pp.col(1) = coslat1 % sin((x1.col(0) * arma::datum::pi)/180.0);
    pp.col(2) = sin((x1.col(1) * arma::datum::pi)/180.0);
    
    arma::vec coslat2 = cos((x2.col(1) * arma::datum::pi)/180.0);
    
    arma::mat temp(x1.n_rows, 3);
    temp.col(0) = coslat2 % cos((x2.col(0) * arma::datum::pi)/180.0);
    temp.col(1) = coslat2 % sin((x2.col(0) * arma::datum::pi)/180.0);
    temp.col(2) = sin((x2.col(1) * arma::datum::pi)/180.0);
    return R * acos(ifelse_arma(pp * trans(temp)));
    
}
