# function which calls all functions to run a BVAR analysis


run_bvar_minnesota <- function(Yraw, p, intercept = FALSE,
                     use_dummies = FALSE,
                     dummy_pars = list(mu = 1.0, gamma = 1.0, prior_mean = NULL),
                     use_flat = FALSE,
                     lag_mean = 1, pi1 = 0.2, pi3 = 1, pi4 = 1000,
                     sigma_vec = NULL,
                     n_draws = 5000, burnin = 1000) {
  #' Runs a Bayesian VAR model with Minnesota style priors and inverse Wishart distributed covariance.
  #' 
  #' Parameters:
  #' - Yraw (matrix): The dataset (T x k), where T is the number of time points and k is the number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - intercept (logical): Whether to include an intercept term in the model (default = FALSE).
  #' - use_dummies (logical): Whether to augment the dataset with dummy observations (default = FALSE).
  #' - dummy_pars (list): Parameters for the dummy observations, including:
  #'   - mu (numeric): Tightness for the initial observational prior (default = 1.0).
  #'   - gamma (numeric): Tightness for the sum-of-coefficients prior (default = 1.0).
  #'   - prior_mean (numeric vector): Mean vector for the prior. If NULL, it is derived from Yraw (default = NULL).
  #' - use_flat (logical): Whether to use flat priors for the coefficients (default = FALSE).
  #' - lag_mean (numeric): Mean for the diagonal elements of the lag coefficients (default = 1).
  #' - pi1 (numeric): Overall tightness parameter for the prior (default = 0.2).
  #' - pi3 (numeric): Lag decay parameter for the prior (default = 1).
  #' - pi4 (numeric): Variance scaling for the intercept (default = 1000).
  #' - sigma_vec (vector): Prior variances for residuals (default = computed from the data).
  #' - n_draws (integer): Total number of MCMC draws (default = 5000).
  #' - burnin (integer): Number of initial MCMC draws to discard as burn-in (default = 1000).
  #' 
  #' Returns:
  #' - A list containing posterior draws:
  #'   - Phi (array): Posterior draws of the VAR coefficient matrix with dimensions (m x k x (n_draws - burnin)).
  #'   - Sigma (array): Posterior draws of the error covariance matrix with dimensions (k x k x (n_draws - burnin)).
  
  
  # build (Y, X), possibly with dummies
  ld_aug <- create_lagged_data(Yraw, p, intercept,
                               use_dummies,
                               dummy_pars)
  Y <- ld_aug$Y
  X <- ld_aug$X
  
  # prepare prior
  prior <- minnesota_prior(Yraw, p, intercept,
                           use_flat,
                           lag_mean, 
                           sigma_vec, pi1, pi3, pi4)
  M0 <- prior$M0
  V0 <- prior$V0
  v0 <- prior$v0
  S0 <- prior$S0
  
  # run the sampler
  res_gibbs <- run_gs_minnesota(
    Y, X,
    M0, V0,
    S0, v0,
    n_draws = n_draws,
    burnin = burnin
  )
  
  return(res_gibbs)
}


run_bvar_hierarch <- function(Yraw, p, intercept, mh_params, hyper_params,
                              pi4 = 1000, lag_mean = 1,
                              n_draws = 5000, burnin = 1000){
  #' Runs a hierarchical Bayesian Vector Autoregression (BVAR) model using Metropolis-Hastings MCMC.
  #'
  #' Parameters:
  #' - Yraw (matrix): The raw observed data matrix with dimensions corresponding to observations and variables.
  #' - p (integer): The lag order for the Vector Autoregression model.
  #' - intercept (logical): Indicates whether an intercept term should be included in the model.
  #' - mh_params (list): A list containing Metropolis-Hastings parameters in the following order:
  #'   - n_thin (integer): The thinning interval for the MCMC samples to reduce autocorrelation.
  #'   - scale_hess (numeric): A scaling factor applied to the Hessian matrix to adjust the proposal covariance.
  #'   - adjust_burn (numeric): A factor determining how much of the burn-in period is used for adjusting the proposal distribution.
  #'   - acc_upper (numeric): The upper threshold for the acceptance rate to adjust the proposal covariance.
  #'   - acc_lower (numeric): The lower threshold for the acceptance rate to adjust the proposal covariance.
  #'   - acc_change (numeric): The amount by which to adjust the proposal covariance based on the acceptance rate.
  #' - hyper_params (list): A list containing hyperparameters for the prior distributions of hyper parameters. 
  #' It should include elements like pi1_param, mu_param, gamma_param, and s2_diag_param, each containing mode, a, and b values.
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
  #'   - acceptance_rate: The overall acceptance rate of the MCMC sampler.
  
  ld <- create_lagged_data(Yraw, p, intercept,
                               use_dummies = FALSE)
  Y <- ld$Y
  X <- ld$X
  
  n_thin <- mh_params[[1]]
  scale_hess <- mh_params[[2]]
  adjust_burn <- mh_params[[3]]
  acc_upper <- mh_params[[4]]
  acc_lower <- mh_params[[5]]
  acc_change <- mh_params[[6]]
  
  res_mh <- run_mh_hierarch(Y, X, p, intercept, hyper_params, n_thin, scale_hess, 
                            adjust_burn, acc_upper, acc_lower, acc_change,
                            pi4 = pi4, lag_mean = lag_mean,
                            n_draws = n_draws, burnin = burnin)
  return(res_mh)
}


run_bvar_ssvs <- function(Yraw, p, 
                          intercept = TRUE, 
                          use_dummies = FALSE,
                          dummy_pars = list(mu = 1.0, gamma = 1.0, prior_mean = NULL),
                          lag_mean = 1,
                          tau0 = 1/10,
                          tau1 = 10,
                          delta_prob = 0.8,
                          n_draws = 5000,
                          burnin = 1000){
  #' Runs a Bayesian VAR model with stochastic search variable selection (SSVS).
  #'
  #' Parameters:
  #' - Yraw (matrix): The raw observed data matrix.
  #' - p (integer): The lag order for the VAR model.
  #' - intercept (logical): Indicates whether an intercept term should be included.
  #' - use_dummies (logical): Indicates whether to include dummy observations.
  #' - dummy_pars (list): A list of parameters for dummy observations, including:
  #'      - mu (numeric): Tightness for the initial observational prior (default = 1.0).
  #'      - gamma (numeric): Tightness for the sum-of-coefficients prior (default = 1.0).
  #'      - prior_mean: An optional prior mean for the coefficients.
  #' - lag_mean (numeric): The lag mean parameter used in the Minnesota prior.
  #' - tau0 (numeric): The scaling factor when a delta element equals 0.
  #' - tau1 (numeric): The scaling factor when a delta element is 1.
  #' - delta_prob (numeric): The prior probability for variable inclusion (used in the log-prior for delta).
  #' - n_draws (integer): Total number of Gibbs sampling iterations.
  #' - burnin (integer): Number of initial iterations to discard as burn-in.
  #'
  #' Returns:
  #' - A list containing:
  #'     delta_draws: A matrix of sampled delta vectors (n_draws x length(delta)).
  #'     Phi: An array of sampled VAR coefficient matrices with dimensions (ncol(X) x k x (n_draws - burnin)).
  #'     Sigma: An array of sampled error covariance matrices with dimensions (k x k x (n_draws - burnin)).
  ld_aug <- create_lagged_data(Yraw, p, intercept,
                               use_dummies,
                               dummy_pars)
  Y <- ld_aug$Y
  X <- ld_aug$X
  
  res_ssvs <- run_gs_ssvs(Y = Y, 
                          X = X, 
                          intercept = intercept, 
                          lag_mean = lag_mean, 
                          tau0 = tau0, 
                          tau1 = tau1, 
                          delta_prob = delta_prob, 
                          n_draws = n_draws, 
                          burnin = burnin)
  
  return(res_ssvs)
}
