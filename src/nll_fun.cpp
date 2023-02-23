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
double logdmultinomCPP(arma::imat x,double size, arma::mat prob  ) {
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

arma::imat calculateY(const int& k, const arma::cube& Yarray, 
                     const int& nrint_k, const int ntint_k){
  arma::imat Y_( nrint_k,ntint_k );
  for(int i = 0; i < nrint_k; ++i){
    for(int j = 0; j < ntint_k; ++j){
      int z = Yarray(i,j,k);
      Y_(i,j) = z
      ;
      
    }
    }
  return Y_;
  }



Rcpp::List  calc_p_mat(const double& tau_k, 
                       const double& phi_k,
                       // arma::mat Y_mat_slice,
                      const int& nrint_k, 
                      const int& ntint_k,
                      const arma::rowvec& tarray_k,
                      const arma::rowvec& rarray_k,
                      const double& max_r_k
){

  arma::mat CDF_binned(nrint_k, ntint_k);//nrint[k],ntint[k]);

  
  for(int j = 0; j < ntint_k; ++j){
    // std::cout << j << "--";
    double tmax = (tarray_k(j)); //max Not sure why R has max of a single value?
    for(int i = 0; i < nrint_k; ++i){
      double upper_r_t = rarray_k(i);
      if(upper_r_t == R_PosInf) upper_r_t = max_r_k;
      double upper_r = upper_r_t;
     
      f_d_CPP f(phi_k,tau_k,tmax);
       double lower_limit = 0.01;
      double err_est;
      int err_code;
      double res = integrate(f,lower_limit,
                                   upper_r, err_est, err_code, 500);//,

      
      CDF_binned(i,j) = res;
    }
  }
  
  // # Difference across distance bins
  arma::mat tmp1=CDF_binned;//(Rcpp::clone(CDF_binned));
  if (tmp1.n_rows>1){
    for (int i=1; i < (nrint_k);++i){
      for(int j = 0; j < tmp1.n_cols; ++j){
        tmp1(i,j) = CDF_binned(i,j) - CDF_binned(i-1,j);
      }
    }
  }
  // # Difference across time bins
  arma::mat p_matrix=tmp1;//(Rcpp::clone(tmp1));
  if (p_matrix.n_cols>1){
    for (int j=1; j < (ntint_k); ++j){
      for(int u=0; u< p_matrix.n_rows; ++u){
        p_matrix(u,j) = tmp1(u,j) - tmp1(u,(j-1));
        
      }
    }
    
  }
  
  // double sum_p_matrix =0;
  // for (int j=0; j < p_matrix.n_cols; ++j){
  //   for(int u=0; u< p_matrix.n_rows; ++u){
  //     sum_p_matrix += p_matrix(u,j);
  //   }
  // }
  
  
  return(Rcpp::List::create(
      Rcpp::Named("p_matrix")=p_matrix,
      // Rcpp::Named("sum_p_matrix")=sum_p_matrix,
      Rcpp::Named("max_r_k")=max_r_k//,
      // Rcpp::Named("Y") = Y
  ));
    
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
double nll_fun(arma::vec& params, const arma::mat& X1,
               const arma::mat& X2,
             const Rcpp::StringVector&  tau_params, 
             const int& nsurvey,
            const arma::cube& Yarray,
            const arma::mat& tarray,
            const arma::mat& rarray,
            const arma::ivec& nrint,
            const arma::ivec& ntint,
            const arma::vec& max_r,
            const arma::ivec& Ysum,
            const arma::vec& nlimit
             ){
  arma::vec sub_v = params(arma::span(0, tau_params.size()-1));
  // arma::vec sub_v = Rcpp::as<arma::vec>(subset);
  arma::vec sub_v_phi = params(arma::span(tau_params.size(), params.size()-1));
  // arma::vec sub_v_phi = Rcpp::as<arma::vec>(subset_phi);

  arma::vec tau = exp(X1 * sub_v );
  arma::vec phi = exp(X2 * sub_v_phi );

  
  
  arma::vec nll(nsurvey);
  // std::cout<< "Starting Loop" <<std::endl;
  for(int k = 0; k < nsurvey; ++k){
        Rcpp::List pmat_out = calc_p_mat(tau(k), 
                                         phi(k),
                                        nrint(k), 
                                         ntint(k),
                                         tarray.row(k),
                                         rarray.row(k),
                                         max_r(k)
        );
        
          arma::mat p_matrix = pmat_out["p_matrix"];
          double sum_p_matrix = arma::accu(p_matrix);//["sum_p_matrix"];
          
          // Calculate Y 
          // std::cout << "a:";
          arma::imat Y = calculateY(k,Yarray, nrint(k), ntint(k));
          arma::mat p_matrix_norm=p_matrix;//(clone(p_matrix));
        // # Normalize
        for(int i = 0; i < p_matrix.n_rows; ++i ){
          for(int j = 0; j < p_matrix.n_cols; ++j ){
            p_matrix_norm(i,j) = p_matrix(i,j)/sum_p_matrix;
          }
        }
        nll(k) = logdmultinomCPP(Y, Ysum(k), p_matrix_norm);//Rcpp::as<arma::mat>(
        
    }// # close loop on k
  //     double nll_sum =0;// <- -sum(nll)
  // for(int i = 0; i< nsurvey; ++i){
  //   nll_sum -= nll[i];
  // }
  double nll_sum =  arma::accu(nll);
  
  if ( R_IsNA(nll_sum) ) return(nlimit[2]); 
  else if (  !arma::is_finite(nll_sum) ) return(nlimit[2]);
  else return(-1 * nll_sum); 
     
    }


// [[Rcpp::export]]
Rcpp::List optim_rcpp(arma::vec& params, const arma::mat& X1,
                     const arma::mat& X2,
                     const Rcpp::StringVector&  tau_params, 
                     const int& nsurvey,
                     const arma::cube& Yarray,
                     const arma::mat& tarray,
                     const arma::mat& rarray,
                     const arma::ivec& nrint,
                     const arma::ivec& ntint,
                     const arma::vec& max_r,
                     const arma::ivec& Ysum,
                     const arma::vec& nlimit,
                     const Rcpp::String method){
  // Extract R's optim function
  Rcpp::Environment stats("package:stats");
  Rcpp::Function optim = stats["optim"];

  // Call the optim function from R in C++
  Rcpp::List opt_results = optim(Rcpp::_["par"]    = params,
                                 // Make sure this function is not exported!
                                 Rcpp::_["fn"]     = Rcpp::InternalFunction(&nll_fun),
                                 Rcpp::_["method"] = method,
                                 Rcpp::_["hessian"] = true,
                                 // Pass in the other parameters as everything
                                 // is scoped environmentally
                                 Rcpp::_["X1"]=X1,
                                 Rcpp::_["X2"]=X2, 
                                 Rcpp::_["tau_params"]=tau_params,
                                 Rcpp::_["nsurvey"]=nsurvey, 
                                 Rcpp::_["Yarray"] = Yarray,
                                 Rcpp::_["tarray"]= tarray,
                                 Rcpp::_["rarray"]=rarray, 
                                 Rcpp::_["nrint"]= nrint, 
                                 Rcpp::_["ntint"]= ntint, 
                                 Rcpp::_["max_r"]= max_r,
                                 Rcpp::_["Ysum"]=Ysum,
                                 Rcpp::_["nlimit"]= nlimit);

  // Extract out the estimated parameter values
  // arma::vec out = Rcpp::as<arma::vec>(opt_results[0]);

  // Return estimated values
  return opt_results;
}
