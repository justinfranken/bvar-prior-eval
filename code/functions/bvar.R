# function which calls all functions to run a BVAR analysis


run_bvar <- function(Yraw, p, intercept = FALSE,
                     use_dummies = FALSE,
                     dummy_pars = list(delta = 1.0, gamma = 1.0, prior_mean = NULL),
                     use_flat = FALSE,
                     lambda = 0.2, alpha = 1.0,
                     sigma_vec = NULL,
                     S0 = diag(ncol(Yraw)), nu0 = ncol(Yraw) + 2,
                     n_draws = 5000, burnin = 1000) {
  #' Runs a Bayesian VAR model with flexible priors and options.
  #' 
  #' Parameters:
  #' - Yraw (matrix): The dataset (T x k), where T is the number of time points and k is the number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - intercept (logical): Whether to include an intercept term in the model (default = FALSE).
  #' - use_dummies (logical): Whether to augment the dataset with dummy observations (default = FALSE).
  #' - dummy_pars (list): Parameters for the dummy observations, including:
  #'   - delta (numeric): Tightness for the initial observational prior (default = 1.0).
  #'   - gamma (numeric): Tightness for the sum-of-coefficients prior (default = 1.0).
  #'   - prior_mean (numeric vector): Mean vector for the prior. If NULL, it is derived from Yraw (default = NULL).
  #' - use_flat (logical): Whether to use flat priors for the coefficients (default = FALSE).
  #' - lambda (numeric): Overall tightness of the Minnesota prior (default = 0.2).
  #' - alpha (numeric): Lag decay exponent for the Minnesota prior (default = 1.0).
  #' - sigma_vec (numeric vector): Prior scale for each variable (length k). If NULL, all scales are set to 1 (default = NULL).
  #' - S0 (matrix): Prior scale matrix for the error covariance (k x k, default = identity matrix of size k).
  #' - nu0 (numeric): Prior degrees of freedom for the error covariance (default = k + 2, where k is the number of variables).
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
                           dummy_pars,
                           lambda, alpha,
                           sigma_vec)
  M0 <- prior$M0
  V0 <- prior$V0
  
  if(use_flat){
    S0 <- diag(ncol(Yraw)) * 1e-3
  }
  
  # run the sampler
  res_gibbs <- run_gs_minnesota(
    Y, X,
    M0, V0,
    S0, nu0,
    n_draws = n_draws,
    burnin = burnin
  )
  
  return(res_gibbs)
}
