# markov chain monte carlo simulation functions to sample from posterior distributions


run_gs_minnesota <- function(Y, X, M0, V0, S0, nu0,
                                  n_draws=5000, burnin=1000) {
  #' Performs Gibbs Sampling to estimate a BVAR model using a Minnesota prior for the VAR coefficients and an inverse-Wishart prior for the error covariance matrix. 
  #' The posterior samples for the VAR coefficients Phi and the error covariance matrix Sigma are stored and returned after discarding the burn-in period.
  #' 
  #' Parameters:
  #' - Y (matrix): Response matrix with dimensions (T_eff x k), where T_eff is the effective number of observations and k is the number of variables.
  #' - X (matrix): Design matrix with dimensions (T_eff x m), where m is the number of predictors (lags and intercept, if included).
  #' - M0 (matrix): Prior mean matrix for the VAR coefficients (m x k).
  #' - V0 (matrix): Prior covariance matrix for the VAR coefficients (m x m).
  #' - S0 (matrix): Prior scale matrix for the error covariance (k x k).
  #' - nu0 (numeric): Prior degrees of freedom for the error covariance matrix.
  #' - n_draws (integer): Total number of MCMC draws (default = 5000).
  #' - burnin (integer): Number of initial MCMC draws to discard as burn-in (default = 1000).
  #' 
  #' Returns:
  #' - A list containing:
  #'   - Phi (array): Posterior draws of the VAR coefficient matrix with dimensions (m x k x (n_draws - burnin)).
  #'   - Sigma (array): Posterior draws of the error covariance matrix with dimensions (k x k x (n_draws - burnin)).
  
  T_eff <- nrow(Y)
  k     <- ncol(Y)
  m     <- ncol(X)
  
  # storage
  Phi_store   <- array(NA, dim=c(m, k, n_draws))
  Sigma_store <- array(NA, dim=c(k, k, n_draws))
  
  # initialize
  Sigma_current <- diag(apply(Y, 2, var))
  
  # -- precompute terms that do not change across iterations
  XX <- crossprod(X)
  V0_inv <- chol2inv(chol(V0))
  V1_inv <- V0_inv + XX
  V1 <- chol2inv(chol(V1_inv))
  
  # posterior mean of Phi, M1
  XTY <- crossprod(X, Y)
  M1  <- V1 %*% (V0_inv %*% M0 + XTY)
  
  # for Sigma draws
  XX_inv <- chol2inv(chol(XX))
  temp <- V0 + XX_inv
  inv_middle_term <- chol2inv(chol(temp))
  
  
  ## -- draw from posterior
  for(iter in 1:n_draws) {
    # draw Phi given Sigma
    Phi_current <- posterior_draw_Phi(M1, Sigma_current, V1)
    
    # draw Sigma given Phi
    Sigma_current <- posterior_draw_Sigma(Y, X, Phi_current, M0, S0, nu0, inv_middle_term)
    
    # store
    Phi_store[,,iter]   <- Phi_current
    Sigma_store[,,iter] <- Sigma_current
  }
  
  # drop burn-in
  keep_idx <- seq.int(from = burnin + 1, to = n_draws)
  Phi_out   <- Phi_store[,, keep_idx, drop=FALSE]
  Sigma_out <- Sigma_store[,, keep_idx, drop=FALSE]
  
  return(list(
    Phi   = Phi_out,
    Sigma = Sigma_out
  ))
}
