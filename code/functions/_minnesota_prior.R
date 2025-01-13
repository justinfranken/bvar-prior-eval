# functions to specify prior beliefs in a Minnesota fashion.


minnesota_prior <- function(Y, p, intercept = FALSE,
                          use_flat = FALSE,
                          lag_mean = 1,
                          s2_diag = NULL,
                          pi1 = 0.2, pi3 = 1, pi4 = 1000) {
  #' Supports either the Minnesota prior or dummy observation-based priors, depending on the input parameters.
  #' 
  #' Parameters:
  #' - Y (matrix): The dataset (T x k), where T is the number of time points and k is the number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - intercept (logical): Whether to include an intercept term in the model (default = FALSE).
  #' - use_flat (logical): Whether to use uninformative M0 and V0 priors instead of the Minnesota prior (default = FALSE).
  #' - lag_mean (numeric): Mean for the diagonal elements of the lag coefficients (default = 1).
  #' - s2_diag (vector): Prior variances for residuals (default = computed from the data).
  #' - pi1 (numeric): Overall tightness parameter for the prior (default = 0.2).
  #' - pi3 (numeric): Lag decay parameter for the prior (default = 1).
  #' - pi4 (numeric): Variance scaling for the intercept (default = 1000).
  #' 
  #' Returns:
  #' - A list containing:
  #'   - M0 (matrix): The prior mean matrix.
  #'   - V0 (matrix): The prior covariance matrix for VAR coefficients (m x m), where (m = k*p) or (1 + k*p) if intercept is included).   #'   - v0 (numeric): The degrees of freedom for the inverse-Wishart prior.
  #'   - S0 (matrix): The scale matrix for the inverse-Wishart prior (k x k).
  #'   - v0 (numeric): The degrees of freedom for the inverse-Wishart prior.
  
  k <- ncol(Y)
  v0 <- k + 2
  
  if(!use_flat) {
    # -- use standard Minnesota prior
    M0 <- helper_build_M0(k, p, intercept, lag_mean)
    temp_cov <- helper_build_cov_matrix(Y = Y, p = p, v = v0,
                                        s2_diag = s2_diag,
                                        intercept = intercept,
                                        pi1 = pi1,
                                        pi3 = pi3,
                                        pi4 = pi4)
    V0 <- temp_cov$V0
    S0 <- temp_cov$S0
    return(list(M0 = M0, V0 = V0, S0 = S0, v0 = v0))
  } else {
    # -- use flat uninformative prior
    m <- if(intercept) 1 + k*p else k*p
    M0_flat <- matrix(0, nrow = m, ncol = k)
    V0_flat <- diag(m) * 1e8
    S0_flat <- diag(k) * 1e-3
    return(list(M0 = M0_flat, V0 = V0_flat, S0 = S0_flat, v0 = v0))
  }
}


helper_build_M0 <- function(k, p, intercept = FALSE, lag_mean = 1){
  #' Builds the prior mean matrix M0 in a Minnesota fashion.
  #' 
  #' Parameters:
  #' - k (integer): Number of variables in the BVAR model.
  #' - p (integer): Number of lags in the BVAR model.
  #' - intercept (logical): Whether to include an intercept term in the model (default = FALSE).
  #' - lag_mean (numeric): Mean for the diagonal elements of the lag coefficients (default = 1). 
  #' Set to < 1 if data is potentially stationary.
  #' 
  #' Returns:
  #' - A matrix M0 with dimensions (m x k), where (m = k*p) (or (1 + k*p) if intercept is included).
  
  m <- if(intercept) 1 + k*p else k*p
  
  M0 <- matrix(0, nrow=m, ncol=k)
  
  if(!intercept){
    M0[1:k, 1:k] = diag(k) * lag_mean
  } else{
    M0[2:(k+1), 1:k] = diag(k) * lag_mean
  }
  
  return(M0)
}


helper_build_cov_matrix <- function(Y, p, v, s2_diag = NULL, intercept = FALSE, pi1 = 0.2, pi3 = 1, pi4 = 1000) {
  #' Builds the covariance matrix and parameters for the normal-Wishart prior.
  #' 
  #' Parameters:
  #' - Y (matrix): The dataset (T x k), where T is the number of time points and k is the number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - v (numeric): The degrees of freedom for the inverse-Wishart prior.
  #' - s2_diag (vector): Prior variances for residuals (default = computed from the data).
  #' - intercept (logical): Whether to include an intercept term in the model (default = FALSE).
  #' - pi1 (numeric): Overall tightness parameter for the prior (default = 0.2).
  #' - pi3 (numeric): Lag decay parameter for the prior (default = 1).
  #' - pi4 (numeric): Variance scaling for the intercept (default = 1000).
  #' 
  #' Returns:
  #' - A list containing:
  #'   - S0 (matrix): The scale matrix for the inverse-Wishart prior (k x k).
  #'   - v0 (numeric): The degrees of freedom for the inverse-Wishart prior.
  #'   - V0 (matrix): The prior covariance matrix for VAR coefficients (m x m), where (m = k*p) or (1 + k*p) if intercept is included).
  
  k <- ncol(Y)
  m <- if (intercept) 1 + k * p else k * p

  if(is.null(s2_diag)) {
    s2_diag <- hhelper_compute_sigma_vec(Y, p)
  }
  
  # S0
  S <- diag(0, k)
  for (i in seq_len(k)) {
    S[i, i] <- (v - k - 1) * s2_diag[i]
  }
  
  # Omega
  diag_elements <- numeric(m)
  for (lag in seq_len(p)) {
    diag_elements[((lag - 1) * k + 1):(lag * k)] <- pi1^2 / ((lag ^ pi3)^2 * s2_diag)
  }
  if (intercept) {
    diag_elements <- c(pi1^2 * pi4^2, diag_elements)
  }
  Omega <- diag(diag_elements, nrow = m)
  
  return(list(S0 = S, V0 = Omega))
}


hhelper_compute_sigma_vec <- function(Y, p) {
  #' Computes residual variance for each variable in an AR model.
  #' 
  #' Parameters:
  #' - Y (matrix): A dataset (T x k), where T is the number of time points and k is the number of variables.
  #' - p (integer): The lag order for the AR(p) model fitted to each variable.
  #' 
  #' Returns:
  #' - A numeric vector of length k, containing the residual variance for each variable.
  
  apply(Y, 2, function(y) {
    fit <- ar(y, order.max = p, aic = FALSE, method = "ols")
    mean(fit$resid^2, na.rm = TRUE) # Residual variance
  })
}

