#ifndef SOLVERS
#define SOLVERS

arma::mat backsolve_arma(arma::mat r, const arma::vec& x, 
                    unsigned int k, 
                    bool transpose);

arma::mat forwardsolve_arma(arma::mat l, const arma::vec& x, 
                       unsigned int k, 
                       bool transpose);

#endif
