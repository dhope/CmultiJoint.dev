#' Calculate Offsets
#'
#' @param fit 
#' @param rarray 
#' @param tarray 
#' @param X1 
#' @param X2 
#'
#' @return
#' @export
#'
#' @examples
calculate_offsets <- function (fit,
                               rarray = rarray,
                               tarray = tarray,
                               X1 = NULL,
                               X2 = NULL) {
  
  # Data used for fitting models
  rarray_fit = fit$input_data$rarray
  tarray_fit = fit$input_data$tarray
  X1_fit = fit$input_data$X1
  X2_fit = fit$input_data$X2
  maxdistint = fit$input_data$maxdistint
  
  nsurvey <- dim(rarray)[1] # Number of surveys
  nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
  ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  
  
  if (!is.null(X1)){
    tau_params <- colnames(X1)
  } else {
    X1 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_params <- colnames(X1)[1] <- "log_tau"
  }
  
  if (!is.null(X2)){
    phi_params <- colnames(X2)
  } else {
    X2 <- matrix(1,nrow = nsurvey,ncol=1)
    phi_params <- colnames(X2)[1] <- "log_phi"
  }
  
  # Calculate maximum distance for integration at each point count
  max_r <- apply(rarray,1,max,na.rm = TRUE)
  max_r[max_r == Inf] <- maxdistint
  
  # Tau and phi
  tau_params <- fit$coefficients[1:length(tau_params)]
  phi_params <- fit$coefficients[(length(tau_params)+1):length(fit$coefficients)]
  
  tau <- poisson("log")$linkinv(drop(X1 %*% tau_params))
  phi <- poisson("log")$linkinv(drop(X2 %*% phi_params))
  
  # Calculate offsets for each survey
  p <- A <- rep(NA,nsurvey)
  
  for (k in 1:nsurvey){
    
    tau_k <- tau[k]
    phi_k <- phi[k]
    
    # Calculate CDF and p
    f_d = function(dmax){
      # print(dmax)
      integrand = substitute(2*pi*dmax *(1-exp(-phi*tmax*exp(-dmax^2/tau^2))),
                             list(phi = phi_k,tau = tau_k,tmax = tmax))
      eval(integrand)
    }
    
    # Calculate CDF
    Y <- Yarray[k,1:nrint[k],1:ntint[k]]
    CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
    for (j in 1:ntint[k]){
      
      tmax = max(tarray[k,j])
      
      for (i in 1:nrint[k]){
        upper_r = rarray[k,i]
        if (upper_r == Inf) upper_r = max_r[k]
        CDF_binned[i,j] = integrate(f_d,lower=0.01,
                                    upper = upper_r,
                                    subdivisions = 500)$value
      }
    }
    
    # Difference across distance bins
    tmp1 = CDF_binned
    if (nrow(tmp1)>1){
      for (i in 2:nrint[k]){
        tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
      }
    }
    
    # Difference across time bins
    p_matrix = tmp1
    if (ncol(p_matrix)>1){
      for (j in 2:ntint[k]){
        p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
      }
    }
    
    p[k] <- sum(p_matrix)
    A[k] <- pi*max_r[k]^2
  }
  
  log_offset <- log(p)
  log_offset
}