
#' cmulti_fit_joint_withR
#'
#' @param Yarray 
#' @param rarray 
#' @param tarray 
#' @param X1 
#' @param X2 
#' @param maxdistint 
#' @param tau_inits 
#' @param phi_inits 
#' @param method 
#' @param ... 
#'
#' @return
cmulti_fit_joint_withR <- function (Yarray, # Array with dimensions (nsurvey x nrint x ntint)
                              rarray, # distance intervals for each point count
                              tarray, # time intervals for each point count
                              X1 = NULL,     # Design matrix for distance modeling
                              X2 = NULL,     # Design matrix for removal modeling
                              maxdistint = 10, # Max distance for numerical integration
                              tau_inits = NULL, 
                              phi_inits = NULL,
                              method = "Nelder-Mead", ...) {
  
  logdmultinom <- function (x, size, prob) {lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))}
  
  ## robust matrix inversion (from detect package)
  .solvenear <- function(x) {
    xinv <- try(solve(x), silent = TRUE)
    if (inherits(xinv, "try-error"))
      xinv <- as.matrix(solve(Matrix::nearPD(x)$mat))
    xinv
  }
  
  # Keep track of the data we put in initially
  # (not sure this part of the script is absolutely necessary... just for helping me sanity check)
  input_data <- list(Yarray = Yarray,
                     rarray = rarray,
                     tarray = tarray,
                     X1 = X1,
                     X2 = X2,
                     maxdistint = maxdistint)
  
  # ----------------------------
  # Only conduct analysis on point counts with non-zero total counts
  # ----------------------------
  
  Ysum <- apply(Yarray,1,sum,na.rm = TRUE)
  Ykeep <- which(Ysum > 0)
  if (length(Ykeep) != length(Ysum)){
    Yarray <- Yarray[Ykeep, , ]
    rarray<- rarray[Ykeep, ]
    tarray<- tarray[Ykeep, ]
    Ysum <- Ysum[Ykeep]
  }
  
  nsurvey <- dim(Yarray)[1] # Number of surveys
  nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
  ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  
  # Format parameters and check they are named via column names in design matrix
  if (!is.null(X1)){
    X1 <- X1[Ykeep, ]
    tau_params <- colnames(X1)
  } else {
    X1 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_params <- colnames(X1)[1] <- "log_tau"
  }
  
  if (!is.null(X2)){
    X2 <- X2[Ykeep, ]
    phi_params <- colnames(X2)
  } else {
    X2 <- matrix(1,nrow = nsurvey,ncol=1)
    phi_params <- colnames(X2)[1] <- "log_phi"
  }
  
  # Calculate maximum distance for integration at each point count
  max_r <- apply(rarray,1,max,na.rm = TRUE)
  max_r[max_r == Inf] <- maxdistint
  
  # Initial values
  if (length(tau_inits) != ncol(X1)) tau_inits <- NULL
  if (is.null(tau_inits)) {
    tau_inits <- rep(0, ncol(X1))
    names(tau_inits) <- tau_params
  }
  
  if (length(phi_inits) != ncol(X2)) phi_inits <- NULL
  if (is.null(phi_inits)) {
    phi_inits <- rep(0, ncol(X2))
    names(phi_inits) <- phi_params
  }
  
  inits <- c(tau_inits,phi_inits)
  # browser()
  # Function to calculate multinomial cell probabilities for each point count
  nll.fun <- function(params) {
    # browser()
    tau <- poisson("log")$linkinv(drop(X1 %*% params[1:length(tau_params)]))
    phi <- poisson("log")$linkinv(drop(X2 %*% params[(length(tau_params)+1):length(params)]))
    
    nll <- rep(0,nsurvey)
    
    for (k in 1:nsurvey){
      
      tau_k <- tau[k]
      phi_k <- phi[k]
      
      # Calculate CDF and p
      f_d = function(dmax){
        # browser()
        integrand = substitute(2*pi*dmax *(1-exp(-phi*tmax*exp(-dmax^2/tau^2))),
                               list(phi = phi_k,tau = tau_k,tmax = tmax))
        eval(integrand)
      }
      
      # Calculate CDF
      Y <- Yarray[k,1:nrint[k],1:ntint[k]] # Data for this survey
      
      CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
      
      for (j in 1:ntint[k]){
        
        tmax = tarray[k,j] # How many minutes have elapsed so far?
        
        for (i in 1:nrint[k]){
          upper_r = rarray[k,i] # what is maximum distance so far
          if (upper_r == Inf) upper_r = max_r[k] # could be simplified earlier in script
          
          # Integrate from 1 m from the observer (seems like integration sometimes crashes if set to 0?)
          CDF_binned[i,j] = integrate(f_d,lower=0.01,
                                      upper = upper_r, 
                                      subdivisions = 500)$value
          # print(CDF_binned[i,j])
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
      # This p_matrix gives us the expected total number of birds detected during the point count
      # if Density = 1, given particular values of phi and tau
      # p_matrix
      
      # Normalize the p_matrix to yield the multinomial cell probabilities
      p_matrix = p_matrix/sum(p_matrix)
      
      # Calculate the multinomial log likelihood for this point count
      nll[k] <- logdmultinom(Y, Ysum[k], p_matrix)
      
    } # close loop on k
    
    # browser()
    nll <- -sum(nll)
    
    if (nll %in% c(NA, NaN, Inf, -Inf)) nlimit[2] else nll
    
  }
  
  nlimit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
  
  res <- optim(inits, nll.fun, method = method, hessian = TRUE)
  
  rval <- list(input_data = input_data,
               convergence = res$convergence,
               coefficients = res$par, 
               vcov = try(.solvenear(res$hessian)), 
               loglik = -res$value)
  
  if (inherits(rval$vcov, "try-error")) rval$vcov <- matrix(NA, length(rval$coefficients), length(rval$coefficients))
  rval$results <- res
  rval
}

