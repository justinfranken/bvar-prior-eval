# functions to specify prior beliefs in a Minnesota fashion.


minnesota_prior <- function(Y, p, intercept = FALSE,
                          use_flat = FALSE,
                          dummy_pars = list(),
                          lambda = 0.2, alpha = 1.0,
                          sigma_vec = NULL) {
  #' Supports either the Minnesota prior or dummy observation-based priors, depending on the input parameters.
  #' 
  #' Parameters:
  #' - Y (matrix): The dataset (T x k), where T is the number of time points and k is the number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - intercept (logical): Whether to include an intercept term in the model (default = FALSE).
  #' - use_flat (logical): Whether to use uninformative M0 and V0 priors instead of the Minnesota prior (default = FALSE).
  #' - dummy_pars (list): Parameters for the dummy observation prior.
  #' - lambda (numeric): Overall tightness parameter for the Minnesota prior (default = 0.2).
  #' - alpha (numeric): Lag decay parameter for the Minnesota prior (default = 1.0).
  #' - sigma_vec (numeric vector): Prior scale for each variable (length k). If NULL, all scales are set to 1 (default = NULL).
  #' 
  #' Returns:
  #' - A list containing:
  #'   - M0 (matrix): The prior mean matrix.
  #'   - V0 (matrix): The prior variance matrix.
  
  k <- ncol(Y)
  
  if(!use_flat) {
    # use standard Minnesota prior
    mn_pr <- build_minnesota_M0_V0(k, p, intercept=intercept,
                                   lambda=lambda, alpha=alpha,
                                   sigma_vec=sigma_vec)
    M0 <- mn_pr$M0
    V0 <- mn_pr$V0
    return(list(M0 = M0, V0 = V0))
  } else {
    # use flat uninformative prior
    m <- if(intercept) 1 + k*p else k*p
    M0_flat <- matrix(0, nrow = m, ncol = k)
    V0_flat <- diag(m) * 1e8
    return(list(M0 = M0_flat, V0 = V0_flat))
  }
}


build_minnesota_M0_V0 <- function(k, p, intercept=FALSE,
                                  lambda=0.2, alpha=1.0,
                                  sigma_vec=NULL) {
  #' Generates the prior mean and variance matrices for a Bayesian VAR model using the Minnesota prior.
  #' 
  #' Parameters:
  #' - k (integer): Number of variables in the VAR model.
  #' - p (integer): Number of lags in the VAR model.
  #' - intercept (logical): Whether to include an intercept in the model (default = FALSE).
  #' - lambda (numeric): Overall tightness of the prior (default = 0.2).
  #' - alpha (numeric): Lag decay exponent, controlling how quickly the prior tightens for higher lags (default = 1.0).
  #' - sigma_vec (numeric vector): Prior scale for each variable (length k). If NULL, all scales are set to 1 (default = NULL).
  #' 
  #' Returns:
  #' - A list containing:
  #'   - M0 (matrix): The prior mean matrix (m x k), where (m = k*p) (or (1 + k*p) if intercept is included).
  #'   - V0 (matrix): The prior variance matrix (m x m), where (m = k*p) (or (1 + k*p) if intercept is included).
  
  if(is.null(sigma_vec)) {
    sigma_vec <- rep(1, k)
  }
  
  m <- if(intercept) 1 + k*p else k*p
  
  # M0
  M0 <- matrix(0, nrow=m, ncol=k)
  if(!intercept){
    M0[1:k, 1:k] = diag(k)
  } else{
    M0[2:(k+1), 1:k] = diag(k)
  }
  
  # V0
  V0 <- matrix(0, m, m)
  
  # intercept row => large variance
  if(intercept){
    V0[1,1] <- 10^6
  }
  
  for(ell in 1:p){
    for(j in 1:k){
      col_index <- (if(intercept) 1 else 0) + (ell-1)*k + j
      for(i in 1:k){
        if(i == j){
          # own-lag
          prior_var_ij <- (lambda^2) / (ell^alpha)
        } else {
          # cross-lag
          prior_var_ij <- (lambda^2) / (ell^alpha) * (sigma_vec[i]^2 / sigma_vec[j]^2)
        }
        V0[col_index, col_index] <- prior_var_ij
      }
    }
  }
  
  return(list(M0 = M0, V0 = V0))
}


create_dummy_observations <- function(Y, p, intercept = FALSE,
                                          delta = 1.0,
                                          gamma = 1.0,
                                          prior_mean = NULL) {
  #' Generates dummy observations for a VAR model, incorporating priors for the sum of coefficients and initial observations.
  #' 
  #' Parameters:
  #' - Y (matrix): The actual data matrix (T x k), where T is the number of observations, and k is the number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - intercept (logical): Whether to include an intercept in the model (default = FALSE).
  #' - delta (numeric): Tightness parameter for the initial observational prior (default = 1.0).
  #' - gamma (numeric): Tightness parameter for the sum-of-coefficients prior (default = 1.0).
  #' - prior_mean (numeric vector): Mean vector for the prior (length k). If NULL, the mean of the first p rows of Y is used (default = NULL).
  #' 
  #' Returns:
  #' - A list containing:
  #'   - Y_dum (matrix): Dummy response observations.
  #'   - X_dum (matrix): Dummy design matrix, including lagged values and optionally an intercept term.
  
  T_real <- nrow(Y)
  k <- ncol(Y)
  
  if (is.null(prior_mean)) {
    prior_mean <- colMeans(Y[1:p, , drop = FALSE], na.rm = TRUE)
  }
  
  # Sum-of-Coefficients Dummy
  Y_sum <- diag((gamma)^-1 * prior_mean)
  X_sum <- matrix(0, nrow = k, ncol = p * k)
  for (i in 1:p) {
    X_sum[, ((i - 1) * k + 1):(i * k)] <- Y_sum
  }
  if (intercept) X_sum <- cbind(rep(0, nrow(X_sum)), X_sum)
  
  # Dummy Initial Observational Prior
  Y_init <- matrix((delta)^-1 * prior_mean, nrow = 1)
  X_init <- matrix(0, nrow = 1, ncol = p * k + as.integer(intercept))
  if (intercept) X_init[1, 1] <- 1 / (delta)^-1
  X_init[1, (1+as.integer(intercept)):(as.integer(intercept) + p * k)] <- (delta)^-1 * rep(prior_mean, times = p)
  
  # Combine Dummies
  Y_dum <- rbind(Y_sum, Y_init)
  X_dum <- rbind(X_sum, X_init)
  
  return(list(Y_dum = Y_dum, X_dum = X_dum))
}
