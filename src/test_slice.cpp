#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int runit(int k, arma::cube Yarray, Rcpp::NumericVector nrint, Rcpp::NumericVector ntint){
  arma::mat Ymat = Yarray.slice(k);
  arma::mat Ymatcl = Ymat(arma::span(0, nrint[k]-1), 
                           arma::span(0, ntint[k]-1));
  Rcpp::NumericMatrix Y = Rcpp::wrap(Ymatcl);
  return(1);
}

/*** R

*/
