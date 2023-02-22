#include <RcppArmadillo.h>
#include <RcppNumerical.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppNumerical)]]
using namespace Numer;

//' logdmultinomCPP
//' 
//' A function to calculate log density
//' of multinomial.
//'
//' @param x A matrix
//' @param size A double
//' @param prob A Matrix
//' @export
// [[Rcpp::export]]
double logdmultinomCPP(arma::mat x,double size, arma::mat prob  ) {
  int nrow_M = x.n_rows;
  int ncol_M = x.n_cols;
  // arma::mat out_mat = x * log(prob ) - lgamma(x + 1);
  double sum_right = 0;
  double left_side = lgamma(size + 1);
  for(int i = 0; i < nrow_M; ++i){
    for(int j = 0; j < ncol_M; ++j){
      double tmp = x(i,j) * log(prob(i,j) ) - lgamma(x(i,j) + 1);
      sum_right += tmp;
    }
  }
    return  left_side + sum_right;

}


// # Calculate CDF and p
class f_d_CPP: public Func
{
private: 
  double phi_k;
  double tau_k;
  double tmax;
public:
  f_d_CPP( double phi_k_, double tau_k_, double tmax_) : phi_k(phi_k_), tau_k(tau_k_), tmax(tmax_) {}
  
  double operator()(const double& dmax) const
  
  {
    return 2*M_PI*dmax *(1-exp(-phi_k*tmax*exp(-pow(dmax,2) / pow(tau_k,2) ) ));
  }
};





Rcpp::List  calc_p_mat(const double tau_k, const double phi_k,
                       // arma::mat Y_mat_slice,
                       int nrint_k, int ntint_k,
                       arma::rowvec tarray_k,
                       arma::rowvec rarray_k,
                       double max_r_k
){

  Rcpp::NumericMatrix CDF_binned(nrint_k, ntint_k);//nrint[k],ntint[k]);

  
  for(int j = 0; j < ntint_k; ++j){
    // std::cout << j << "--";
    double tmax = (tarray_k(j)); //max Not sure why R has max of a single value?
    for(int i = 0; i < nrint_k; ++i){
      double upper_r_t = rarray_k(i);
      if(upper_r_t == R_PosInf) upper_r_t = max_r_k;
      const double upper_r = upper_r_t;
     
      f_d_CPP f(phi_k,tau_k,tmax);
      const double lower_limit = 0.01;
      double err_est;
      int err_code;
      const double res = integrate(f,lower_limit,
                                   upper_r, err_est, err_code, 500);//,

      
      CDF_binned(i,j) = res;
    }
  }
  
  // # Difference across distance bins
  Rcpp::NumericMatrix tmp1(Rcpp::clone(CDF_binned));
  if (tmp1.nrow()>1){
    for (int i=1; i < (nrint_k);++i){
      for(int j = 0; j < tmp1.ncol(); ++j){
        tmp1(i,j) = CDF_binned(i,j) - CDF_binned(i-1,j);
      }
    }
  }
  // # Difference across time bins
  Rcpp::NumericMatrix p_matrix(Rcpp::clone(tmp1));
  if (p_matrix.ncol()>1){
    for (int j=1; j < (ntint_k); ++j){
      for(int u=0; u< p_matrix.nrow(); ++u){
        p_matrix(u,j) = tmp1(u,j) - tmp1(u,(j-1));
        
      }
    }
    
  }
  
  double sum_p_matrix =0;
  for (int j=0; j < p_matrix.ncol(); ++j){
    for(int u=0; u< p_matrix.nrow(); ++u){
      sum_p_matrix += p_matrix(u,j);
    }
  }
  
  return(Rcpp::List::create(
      Rcpp::Named("p_matrix")=p_matrix,
      Rcpp::Named("sum_p_matrix")=sum_p_matrix,
      Rcpp::Named("max_r_k")=max_r_k//,
      // Rcpp::Named("Y") = Y
  ));
    
}


arma::mat calculateY(int k, arma::cube Yarray, Rcpp::NumericVector nrint, Rcpp::NumericVector ntint){
  arma::mat Y_( nrint[k],ntint[k] );
  for(int i = 0; i < nrint[k]; ++i){
    for(int j = 0; j < ntint[k]; ++j){
      int z = Yarray(i,j,k);
      Y_(i,j) = z;
      
    }
    
  }
  return Y_;
}


//' Negative log likelihood function
//' 
//' A function to calculate negative log likelihood.
//'
//' @param params
//' @param X1
//' @param X2
//' @param tau_params
//' @param nsurvey
//' @param Yarray
//' @param tarray
//' @param rarray
//' @param nrint
//' @param ntint
//' @param max_r
//' @param Ysum
//' @param nlimit
//' @export
// [[Rcpp::export]]
double nll_fun(Rcpp::NumericVector params, arma::mat X1, arma::mat X2,
             Rcpp::StringVector  tau_params, int nsurvey,
             arma::cube Yarray,
             arma::mat tarray,
             arma::mat rarray,
             Rcpp::NumericVector nrint,
             Rcpp::NumericVector ntint,
             Rcpp::NumericVector max_r,
             Rcpp::NumericVector Ysum,
             Rcpp::NumericVector nlimit
             ){
  Rcpp::NumericVector subset = params[Rcpp::Range(0, tau_params.size()-1)];
  arma::vec sub_v = Rcpp::as<arma::vec>(subset);
  Rcpp::NumericVector subset_phi = params[Rcpp::Range(tau_params.size(), params.size()-1)];
  arma::vec sub_v_phi = Rcpp::as<arma::vec>(subset_phi);

  // std::cout<< arma::size(Yarray) << std::endl;
  arma::vec tau = exp(X1 * sub_v );
  arma::vec phi = exp(X2 * sub_v_phi );

  Rcpp::NumericVector nll(nsurvey);
  // std::cout<< "Starting Loop" <<std::endl;
  for(int k = 0; k < nsurvey; ++k){
        // std::cout<< k <<std::endl;
        // Call to duplicate calls
        const double tau_k=tau(k);
        // std::cout<< "tau ";
        const double phi_k=phi(k);
        // std::cout<< "phi ";
        int nrint_k=nrint(k);
        // std::cout<< "nrint ";
        int ntint_k=ntint(k);
        // std::cout<< "ntint ";
        arma::rowvec tarray_k=tarray.row(k);
        // std::cout<< "tarray ";
        arma::rowvec rarray_k=rarray.row(k);
        // std::cout<< "rarray ";
        double max_r_k=max_r(k);
        // std::cout<< "max_r ";
        Rcpp::List pmat_out = calc_p_mat(tau_k, 
                                         phi_k,
                                        nrint_k, 
                                         ntint_k,
                                         tarray_k,
                                         rarray_k,
                                         max_r_k
        );
        
          Rcpp::NumericMatrix p_matrix = pmat_out["p_matrix"];
          double sum_p_matrix = pmat_out["sum_p_matrix"];
          
          // Calculate Y 
          arma::mat Y = calculateY(k,Yarray, nrint, ntint);
    

        
        
        Rcpp::NumericMatrix p_matrix_norm(clone(p_matrix));
        // # Normalize
        for(int i = 0; i < p_matrix.nrow(); ++i ){
          for(int j = 0; j < p_matrix.ncol(); ++j ){
            p_matrix_norm(i,j) = p_matrix(i,j)/sum_p_matrix;
          }
        }
        
        
        nll[k] = logdmultinomCPP(Y, Ysum[k], Rcpp::as<arma::mat>(p_matrix_norm));
        // Y.reset();
        // delete Y;
    }// # close loop on k
      double nll_sum =0;// <- -sum(nll)
  for(int i = 0; i< nsurvey; ++i){
    nll_sum -= nll[i];
  }
  
  if ( R_IsNA(nll_sum) ) return(nlimit[2]); 
  else if (  !arma::is_finite(nll_sum) ) return(nlimit[2]);
  else return(nll_sum); 
     
    }

