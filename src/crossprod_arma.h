#ifndef CROSSPROD_ARMA_R
#define CROSSPROD_ARMA_R

arma::mat self_crossprod_arma(const arma::mat& x);

arma::mat crossprod_arma(const arma::mat& x, const arma::mat& y);

arma::mat self_tcrossprod_arma(const arma::mat& x);

arma::mat tcrossprod_arma(const arma::mat& x, const arma::mat& y);

#endif
