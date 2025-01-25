# functions to specify prior beliefs for BVAR estimation.


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
  #' - s2_diag (vector): Prior standard deviations for residuals (default = computed from the data).
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
    temp_cov <- helper_build_cov_matrix_minnesota(Y = Y, p = p, v0 = v0,
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


hierarchical_prior <- function(Y, p, intercept = FALSE,
                            lag_mean = 1,
                            s2_diag = NULL,
                            pi1 = 0.2, pi4 = 1000) {
  #' Sets prior beliefs for hierarchical BVAR estimation.
  #' 
  #' Parameters:
  #' - Y (matrix): The dataset (T x k), where T is the number of time points and k is the number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - intercept (logical): Whether to include an intercept term in the model (default = FALSE).
  #' - lag_mean (numeric): Mean for the diagonal elements of the lag coefficients (default = 1).
  #' - s2_diag (vector): Prior standard deviations for residuals (default = computed from the data).
  #' - pi1 (numeric): Overall tightness parameter for the prior (default = 0.2).
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
  
  # -- use Minnesota prior
  M0 <- helper_build_M0(k, p, intercept, lag_mean)
  temp_cov <- helper_build_cov_matrix_hierachical(Y = Y, p = p, v0 = v0,
                                                s2_diag = s2_diag,
                                                intercept = intercept,
                                                pi1 = pi1,
                                                pi4 = pi4)
  V0 <- temp_cov$V0
  S0 <- temp_cov$S0
  return(list(M0 = M0, V0 = V0, S0 = S0, v0 = v0))
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
    diag(M0[1:k, ]) <- lag_mean
  } else {
    diag(M0[2:(k+1), ]) <- lag_mean
  }
  
  return(M0)
}


helper_build_cov_matrix_minnesota <- function(Y, p, v0, s2_diag = NULL, intercept = FALSE, pi1 = 0.2, pi3 = 1, pi4 = 1000) {
  #' Builds the covariance matrix and parameters for the normal-Wishart prior.
  #' 
  #' Parameters:
  #' - Y (matrix): The dataset (T x k), where T is the number of time points and k is the number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - s2_diag (vector): Prior standard deviations for residuals (default = computed from the data).
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
  S0 <- diag((v0 - k - 1) * s2_diag^2, nrow = k)

  # V0
  repeated_s2 <- rep(s2_diag, times = p)
  repeated_lags2 <- rep((1:p)^pi3, each = k)
  diag_elements  <- pi1^2 / (repeated_lags2 * repeated_s2)^2
  
  if (intercept) {
    diag_elements <- c(pi1^2 * pi4^2, diag_elements)
  }
  
  # construct the diagonal matrix V0
  V0 <- diag(diag_elements, nrow = m, ncol = m)
  
  return(list(S0 = S0, V0 = V0))
}


helper_build_cov_matrix_hierachical <- function(Y, p, v0, s2_diag = NULL, intercept = FALSE, pi1 = 0.2, pi4 = 1000) {
  #' Builds the covariance matrix and parameters for hierarchical estimation.
  #' 
  #' Parameters:
  #' - Y (matrix): The dataset (T x k), where T is the number of time points and k is the number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - s2_diag (vector): Prior standard deviations for residuals (default = computed from the data).
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
  
  if (is.null(s2_diag)) {
    s2_diag <- hhelper_compute_sigma_vec(Y, p)
  }
  
  # S0
  S0 <- diag((v0 - k - 1) * s2_diag, nrow = k)
  
  # V0
  repeated_s2 <- rep(s2_diag, times = p)
  repeated_lags2 <- rep((1:p)^2, each = k)
  diag_elements  <- (pi1^2 * (v0 - k - 1)) / (repeated_lags2 * repeated_s2)

  if (intercept) {
    diag_elements <- c(pi4^2 * (v0 - k - 1), diag_elements)
  }
  
  # construct the diagonal matrix V0
  V0 <- diag(diag_elements, nrow = m, ncol = m)
  
  return(list(S0 = S0, V0 = V0))
}


create_dummy_observations <- function(Y, p, intercept = FALSE,
                                      mu = 1.0,
                                      gamma = 1.0,
                                      prior_mean = NULL) {
  #' Generates dummy observations for a BVAR model, incorporating priors for the sum of coefficients and initial observations.
  #' 
  #' Parameters:
  #' - Y (matrix): The actual data matrix (T x k), where T is the number of observations, and k is the number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - intercept (logical): Whether to include an intercept in the model (default = FALSE).
  #' - mu (numeric): Tightness parameter for the initial observational prior (default = 1.0).
  #' - gamma (numeric): Tightness parameter for the sum-of-coefficients prior (default = 1.0).
  #' - prior_mean (numeric vector): Mean vector for the prior (length k). If NULL, the mean of the first p rows of Y is used (default = NULL).
  #' 
  #' Returns:
  #' - A list containing:
  #'   - Y_dum (matrix): Dummy response observations.
  #'   - X_dum (matrix): Dummy design matrix, including lagged values and optionally an intercept term.
  
  k <- ncol(Y)
  
  if (is.null(prior_mean)) {
    prior_mean <- colMeans(Y[1:p, , drop = FALSE], na.rm = TRUE)
  }

  # sum of coefficients dummy
  Y_sum <- diag((1/gamma) * prior_mean, nrow = k)
  X_sum <- do.call(cbind, replicate(p, Y_sum, simplify = FALSE))

  if (intercept) {
    X_sum <- cbind(rep(0, k), X_sum)
  }
  
  # initial observations dummy
  Y_init <- matrix((1/mu) * prior_mean, nrow = 1)
  
  if (intercept) {
    X_init <- matrix(c(1 / mu, (1/mu) * rep(prior_mean, p)), nrow = 1)
  } else {
    X_init <- matrix((1/mu) * rep(prior_mean, p), nrow = 1)
  }
  
  # Combine into final dummy observations
  Y_dum <- rbind(Y_sum, Y_init)
  X_dum <- rbind(X_sum, X_init)
  
  return(list(Y_dum = Y_dum, X_dum = X_dum))
}


hhelper_compute_sigma_vec <- function(Y, p) {
  #' Computes residual standard deviation for each variable in an AR(p) model.
  #' 
  #' Parameters:
  #' - Y (matrix): A dataset (T x k), where T is the number of time points and k is the number of variables.
  #' - p (integer): Number of lags for the AR model used.
  #' 
  #' Returns:
  #' - A numeric vector of length k, containing the residual standard deviation for each variable.
  
  apply(Y, 2, function(y) {
    fit <- ar(y, order.max = p, aic = FALSE, method = "ols")
    # residual standard deviation
    sqrt(mean(fit$resid^2, na.rm = TRUE))
  })
}

