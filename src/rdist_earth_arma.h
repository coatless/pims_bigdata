#ifndef RDIST_EARTH_ARMA
#define RDIST_EARTH_ARMA

arma::mat rdist_earth1(arma::mat x1, 
                bool miles, double R);
                
arma::mat rdist_earth2(arma::mat x1, arma::mat x2,
                bool miles, double R);

#endif
