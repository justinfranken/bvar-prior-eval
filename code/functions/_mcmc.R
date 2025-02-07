# markov chain monte carlo simulation functions to sample from posterior distributions


run_gs_minnesota <- function(Y, X, M0, V0, S0, nu0,
                                  n_draws = 5000, burnin = 1000) {
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
  k <- ncol(Y)
  m <- ncol(X)
  
  # storage
  Phi_store <- array(NA, dim=c(m, k, n_draws))
  Sigma_store <- array(NA, dim=c(k, k, n_draws))

  # -- precompute terms that do not change across iterations
  XX <- crossprod(X)
  V0_inv <- chol2inv(chol(V0))
  V1_inv <- V0_inv + XX
  V1 <- chol2inv(chol(V1_inv))
  
  # posterior mean of Phi, M1
  XTY <- crossprod(X, Y)
  M1 <- V1 %*% (V0_inv %*% M0 + XTY)
  
  # - for Sigma draws
  XX_inv <- chol2inv(chol(XX))
  Phi_ols <- crossprod(XX_inv,XTY)
  temp <- V0 + XX_inv
  diffPhi <- M0 - Phi_ols
  inv_middle_term <- chol2inv(chol(temp))
  E <- Y - X %*% Phi_ols
  SSE <- crossprod(E)
  second_term <- crossprod(diffPhi, inv_middle_term %*% diffPhi)
  
  # posterior scale
  S1  <- S0 + SSE + second_term
  
  # posterior dof
  nu1 <- nu0 + nrow(Y)


  ## -- draw from posterior
  for(iter in 1:n_draws) {
    # draw Sigma
    Sigma_current <- posterior_draw_Sigma(S1, nu1)
    
    # draw Phi given Sigma
    Phi_current <- posterior_draw_Phi(M1, Sigma_current, V1)
    
    # store
    Phi_store[,,iter] <- Phi_current
    Sigma_store[,,iter] <- Sigma_current
  }
  
  # drop burn-in
  keep_idx <- seq.int(from = burnin + 1, to = n_draws)
  Phi_out <- Phi_store[,, keep_idx, drop = FALSE]
  Sigma_out <- Sigma_store[,, keep_idx, drop = FALSE]
  
  return(list(
    Phi = Phi_out,
    Sigma = Sigma_out
  ))
}


run_mh_hierarch <- function(Y, X, p, intercept, hyper_params, 
                            n_thin, scale_hess, adjust_burn, acc_upper,
                            acc_lower, acc_change,
                            pi4 = 1000, lag_mean = 1,
                            n_draws = 5000, burnin = 1000){
  #' This function performs Metropolis-Hastings sampling to estimate the posterior
  #' distribution of hyper parameters in a hierarchical Bayesian framework. It initializes
  #' the MCMC loop by optimizing the posterior density, adjusts the proposal distribution
  #' based on acceptance rates during the burn-in phase, and stores the sampled draws
  #' after thinning. The function handles parameter constraints and ensures numerical
  #' stability throughout the sampling process.
  #'
  #' Parameters:
  #' - Y (matrix): The observed data matrix with dimensions T x k, where T is the number of observations and k is the number of variables.
  #' - X (matrix): The predictor matrix with dimensions T x n, where n is the number of predictors.
  #' - p (integer): The lag order or another relevant integer parameter used in the hierarchical_prior function.
  #' - intercept (logical): Indicates whether an intercept term is included in the model.
  #' - hyper_params (list): A list containing hyperparameters for the prior distributions of hyper parameters. 
  #'   It should include elements like pi1_param, mu_param, gamma_param, and s2_diag_param, each containing mode, a, and b values.
  #' - n_thin (integer): The thinning interval for the MCMC samples to reduce autocorrelation.
  #' - scale_hess (numeric): A scaling factor applied to the Hessian matrix to adjust the proposal covariance.
  #' - adjust_burn (numeric): A factor determining how much of the burn-in period is used for adjusting the proposal distribution.
  #' - acc_upper (numeric): The upper threshold for the acceptance rate to adjust the proposal covariance.
  #' - acc_lower (numeric): The lower threshold for the acceptance rate to adjust the proposal covariance.
  #' - acc_change (numeric): The amount by which to adjust the proposal covariance based on the acceptance rate.
  #' - pi4 (numeric, default = 1000): Variance scaling for the intercept.
  #' - lag_mean (integer, default = 1): Mean for the diagonal elements of the lag coefficients.
  #' - n_draws (integer, default = 5000): The total number of MCMC draws to generate.
  #' - burnin (integer, default = 1000): The number of initial MCMC iterations to discard as burn-in.
  #'
  #' Returns:
  #' - A list containing:
  #'   - Phi: An array of sampled Phi coefficient matrices with dimensions (k*p x k x (n_draws - burnin)/n_thin).
  #'   - Sigma: An array of sampled Sigma covariance matrices with dimensions (k x k x (n_draws - burnin)/n_thin).
  #'   - hyper_parameters: A matrix of sampled hyper parameters with dimensions ((n_draws - burnin)/n_thin x (3 + k)).
  #'   - optim_values: The optimized parameter values from the initial optimization step.
  #'   - log_post_vec: A vector of log-posterior values corresponding to the sampled hyper parameters.
  #'   - acceptance_rate: The overall acceptance rate of the MCMC sampler after the burn-in.
  
  # ----------------------------------------------------------------------------
  # initialize MCMC loop
  # ----------------------------------------------------------------------------
  
  # --- run optim
  optim_par <- c(hyper_params$pi1_param$mode, hyper_params$mu_param$mode, 
                 hyper_params$gamma_param$mode, hyper_params$s2_diag_param$mode)
  lower_bound <- c(1e-04, 1e-04, 1e-04, hyper_params$s2_diag_param$mode / 100)
  upper_bound <- c(5, 50, 50, hyper_params$s2_diag_param$mode * 100)
  
  res <- optim(
    par = optim_par,
    fn = optim_log_post_delta,
    method = "L-BFGS-B",
    lower = lower_bound,
    upper = upper_bound,
    control = list("fnscale" = -1),
    Y = Y, X = X, p = p, intercept = intercept, pi4 = pi4,
    hyper_params = hyper_params, lower_bound = lower_bound,
    upper_bound = upper_bound
  )
  
  H <- diag(res$par) * scale_hess
  J <- (upper_bound - lower_bound) * exp(res$par) / (1+exp(res$par))^2
  J <- diag(J)
  HH <- J %*% H %*% t(J)
  
  # make sure HH is positive definite
  HH_eig <- eigen(HH)
  HH_eig[["values"]] <- abs(HH_eig[["values"]])
  HH <- HH_eig
  
  # --- draw first proposal
  while (TRUE) {
    delta_draw <- as.vector(rmvn_proposal(n = 1, mean = res$par, sigma = HH))
    prior_draw <- hierarchical_prior(Y, p, intercept, 
                                     lag_mean = 1, 
                                     s2_diag = delta_draw[4:(k+3)], 
                                     pi1 = delta_draw[1], 
                                     pi4 = pi4)
    dmy_draw <- create_dummy_observations(Y, p, intercept = intercept,
                                          mu = delta_draw[2], 
                                          gamma = delta_draw[3])
    Y_aug_draw <- rbind(dmy_draw$Y_dum, Y)
    X_aug_draw <- rbind(dmy_draw$X_dum, X)
    log_post_draw <- log_posterior_delta(Y = Y_aug_draw, 
                                         X = X_aug_draw, 
                                         M0 = prior_draw$M0, 
                                         V0 = prior_draw$V0, 
                                         S0 = prior_draw$S0, 
                                         v0 = prior_draw$v0, 
                                         pi1_val = delta_draw[1], 
                                         mu_val = delta_draw[2], 
                                         gamma_val = delta_draw[3], 
                                         s2_diag_val = delta_draw[4:(k+3)], 
                                         hyper_params = hyper_params,
                                         lower_bound = lower_bound,
                                         upper_bound = upper_bound)
    if(log_post_draw > -1e16) {break}
  }
  
  # --- storage
  priors_list <- list()
  hyper_vals <- matrix(NA, nrow = (n_draws - burnin)/n_thin, ncol = 3 + k)
  log_post_vec <- rep(0, (n_draws - burnin)/n_thin)
  Phi_store <- array(NA, dim=c(ncol(X), k, (n_draws - burnin)/n_thin))
  Sigma_store <- array(NA, dim=c(k, k, (n_draws - burnin)/n_thin))
  accepted_adj <- 0 -> accepted
  
  # ----------------------------------------------------------------------------
  # run MCMC loop
  # ----------------------------------------------------------------------------
  
  for (i in seq.int(1 - burnin, n_draws - burnin)) {
    
    # --- proposal draw
    delta_temp <- as.vector(rmvn_proposal(n = 1, mean = res$par, sigma = HH))
    prior_temp <- hierarchical_prior(Y, p, intercept, 
                                     lag_mean = 1, 
                                     s2_diag = delta_temp[4:(k+3)], 
                                     pi1 = delta_temp[1], 
                                     pi4 = pi4)
    dmy_temp <- create_dummy_observations(Y, p, intercept = intercept,
                                          mu = delta_temp[2], 
                                          gamma = delta_temp[3])
    Y_aug_temp <- rbind(dmy_temp$Y_dum, Y)
    X_aug_temp <- rbind(dmy_temp$X_dum, X)
    log_post_temp <- log_posterior_delta(Y = Y_aug_temp, 
                                         X = X_aug_temp, 
                                         M0 = prior_temp$M0, 
                                         V0 = prior_temp$V0, 
                                         S0 = prior_temp$S0, 
                                         v0 = prior_temp$v0, 
                                         pi1_val = delta_temp[1], 
                                         mu_val = delta_temp[2], 
                                         gamma_val = delta_temp[3], 
                                         s2_diag_val = delta_temp[4:(k+3)], 
                                         hyper_params = hyper_params,
                                         lower_bound = lower_bound,
                                         upper_bound = upper_bound)
    
    # --- accept or reject
    if(runif(1) < exp(log_post_temp - log_post_draw)){
      log_post_draw <- log_post_temp
      delta_draw <- delta_temp
      prior_draw <- prior_temp
      Y_aug_draw <- Y_aug_temp
      X_aug_draw <- X_aug_temp
      accepted_adj <- accepted_adj + 1
      if(i > 0) {accepted <- accepted + 1}
    }
    
    # adjust hessian in burn-in phase every 10th draw
    if(i <= - as.integer(adjust_burn * burnin) && (i + burnin) %% 10 == 0){
      acc_rate <- accepted_adj / (i + burnin)
      if(acc_rate < acc_lower){
        HH$values <- HH$values * (1 - acc_change)
      } else if(acc_rate > acc_upper){
        HH$values <- HH$values * (1 + acc_change)
      }
    }
    
    # --- store draws
    if(i > 0 && i %% n_thin == 0){
      hyper_vals[i / n_thin,] <- delta_draw
      log_post_vec[i / n_thin] <- log_post_draw
      
      # -- draw posterior Sigma and Phi
      XX <- crossprod(X_aug_draw)
      V0_inv <- chol2inv(chol(prior_draw$V0))
      V1_inv <- V0_inv + XX
      V1 <- chol2inv(chol(V1_inv))
      
      # posterior mean of Phi, M1
      XTY <- crossprod(X_aug_draw, Y_aug_draw)
      M1 <- V1 %*% (V0_inv %*% prior_draw$M0 + XTY)
      
      # for Sigma draws
      XX_inv <- chol2inv(chol(XX))
      Phi_ols <- crossprod(XX_inv,XTY)
      temp <- prior_draw$V0 + XX_inv
      inv_middle_term <- chol2inv(chol(temp))
      E <- Y_aug_draw - X_aug_draw %*% Phi_ols
      SSE <- crossprod(E)
      diffPhi <- prior_draw$M0 - Phi_ols
      second_term <- crossprod(diffPhi, inv_middle_term %*% diffPhi)
      
      # posterior scale
      S1  <- prior_draw$S0 + SSE + second_term
      
      # posterior dof
      nu1 <- prior_draw$v0 + nrow(Y_aug_draw)
      
      # draw Sigma
      Sigma_current <- posterior_draw_Sigma(S1, nu1)
      
      # draw Phi given Sigma
      Phi_current <- posterior_draw_Phi(M1, Sigma_current, V1)
      
      # store coefficients
      Phi_store[,,i / n_thin] <- Phi_current
      Sigma_store[,,i / n_thin] <- Sigma_current
    }
  }
  
  acceptance_rate <- accepted / (n_draws - burnin)
  
  return(list(
    Phi = Phi_store,
    Sigma = Sigma_store,
    hyper_parameters = hyper_vals,
    optim_values = res$par,
    log_post_vec = log_post_vec,
    acceptance_rate = acc_rate
  ))
}


run_gs_ssvs <- function(Y, X, intercept = TRUE, lag_mean = 1, 
                        tau0 = 0.01, tau1 = 10, delta_prob = 0.8,
                        n_draws = 5000, burnin = 1000) {
  #' Runs a Gibbs sampler for stochastic search variable selection (SSVS) in a Bayesian VAR.
  #'
  #' Parameters:
  #' - Y (matrix): The observed data matrix with dimensions T_eff x k (T_eff observations and k variables).
  #' - X (matrix): The design matrix (including lags) with dimensions T_eff x m.
  #' - intercept (logical): Indicates whether an intercept term is included in the model.
  #' - lag_mean (numeric): The lag mean parameter used in forming the Minnesota prior.
  #' - tau0 (numeric): Scaling factor used when the corresponding delta is 0.
  #' - tau1 (numeric): Scaling factor used when the corresponding delta is 1.
  #' - delta_prob (numeric): The prior probability for inclusion (used in the log-prior for delta).
  #' - n_draws (integer): Total number of Gibbs sampling iterations.
  #' - burnin (integer): Number of initial draws to discard as burn-in.
  #'
  #' Returns:
  #' - A list containing:
  #'   - Phi: An array of sampled VAR coefficient matrices with dimensions (m x k x (n_draws - burnin)).
  #'   - Sigma: An array of sampled error covariance matrices with dimensions (k x k x (n_draws - burnin)).
  #'   - delta_draws: A matrix of sampled delta vectors (n_draws x (k*p) or if intercept is TRUE n_draws x (k*p +1)).
  
  T_eff <- nrow(Y)
  k <- ncol(Y)
  m <- ncol(X)
  p <- (m - as.integer(intercept)) / k
  
  XX <- crossprod(X)
  XX_inv <- solve(XX)
  Xy <- crossprod(X, Y)
  
  se_ols <- compute_sigma_ols_VAR(Y, X, XX_inv)
  
  # priors
  temp <- minnesota_prior(Y, p , intercept = intercept, lag_mean = lag_mean)
  S_lowerbar <- temp$S0
  v_lowerbar <- temp$v0
  M0 <- temp$M0
  v_upperbar <- v_lowerbar + T_eff
  Phi_ols <- crossprod(XX_inv, Xy)
  residuals <- Y - X %*% Phi_ols
  diff_Phi <- M0 - Phi_ols
  S <- crossprod(residuals)
  
  # initialize delta
  delta_init <- rep(1, ncol(X))  # or some random choice
  
  # precompute everything for the initial delta
  precomp <- precompute_from_delta(delta = delta_init,
                                   T = T_eff, k = k,
                                   S_lowerbar = S_lowerbar, v_lowerbar = v_lowerbar,
                                   tau0 = tau0, tau1 = tau1, se_ols = se_ols,
                                   XX = XX, XX_inv = XX_inv,
                                   diff_Phi = diff_Phi, S = S,
                                   delta_prob = delta_prob)
  
  # storage
  store_delta <- matrix(NA, nrow = n_draws, ncol = length(delta_init))
  Phi_store <- array(NA, dim=c(m, k, n_draws))
  Sigma_store <- array(NA, dim=c(k, k, n_draws))
  
  # gibbs sampling
  for (iter in seq_len(n_draws)) {
    # draw delta and new priors / posteriors
    precomp <- update_delta_rank1(precomp,
                                  T = T_eff, k = k,
                                  S_lowerbar = S_lowerbar, v_lowerbar = v_lowerbar,
                                  tau0 = tau0, tau1 = tau1, se_ols = se_ols,
                                  XX = XX, XX_inv = XX_inv,
                                  diff_Phi = diff_Phi, S = S,
                                  delta_prob = delta_prob)
    # store delta
    store_delta[iter, ] <- precomp$delta
    
    # draw coefficients Phi and Sigma
    if(iter > burnin){
      # draw Sigma
      Sigma_current <- posterior_draw_Sigma(precomp$S_upperbar, v_upperbar)
      
      # draw Phi given Sigma
      M1 <- precomp$P_inv %*% (precomp$inv_Omega_lowerbar %*% M0 + Xy)
      Phi_current <- posterior_draw_Phi(M1, Sigma_current, precomp$P_inv)
      
      # store coefficients
      Phi_store[,,iter] <- Phi_current
      Sigma_store[,,iter] <- Sigma_current
    }
  }
  
  # drop burn-in
  keep_idx <- seq.int(from = burnin + 1, to = n_draws)
  Phi_out <- Phi_store[,, keep_idx, drop = FALSE]
  Sigma_out <- Sigma_store[,, keep_idx, drop = FALSE]
  
  return(list(Phi = Phi_out,
              Sigma = Sigma_out,
              delta_draws = store_delta
  ))
}
