# function which calls all functions to run a BVAR analysis


run_bvar <- function(Yraw, p, intercept = FALSE,
                     use_dummies = FALSE,
                     dummy_pars = list(delta = 1.0, gamma = 1.0, prior_mean = NULL),
                     use_flat = FALSE,
                     lag_mean = 1, pi1 = 0.2, pi3 = 1, pi4 = 1000,
                     sigma_vec = NULL,
                     n_draws = 5000, burnin = 1000) {
  #' Runs a Bayesian VAR model with flexible priors and options.
  #' 
  #' Parameters:
  #' - Yraw (matrix): The dataset (T x k), where T is the number of time points and k is the number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - intercept (logical): Whether to include an intercept term in the model (default = FALSE).
  #' Note: If working with an intercept, scaling up % values to larger numbers (0.05 -> 5) can resolve non-positive matrix issues.
  #' - use_dummies (logical): Whether to augment the dataset with dummy observations (default = FALSE).
  #' - dummy_pars (list): Parameters for the dummy observations, including:
  #'   - delta (numeric): Tightness for the initial observational prior (default = 1.0).
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
