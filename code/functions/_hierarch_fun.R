# functions needed for hierarchical bvar analysis


calc_gamma_pdf_params <- function(m, s) {
  #' This function determines the shape (a) and rate (b) parameters of a Gamma
  #' distribution by solving a quadratic equation derived from the provided mode (m)
  #' and scale (s). 
  #'
  #' Parameters:
  #' - m (numeric): The mode of the Gamma distribution.
  #' - s (numeric): The scale parameter, related to the standard deviation of the distribution.
  #'
  #' Returns:
  #' - A list containing:
  #'   - a (numeric): The shape parameter of the Gamma distribution.
  #'   - b (numeric): The rate parameter of the Gamma distribution.
  #'   - mode (numeric): The mode m as provided in the input.
  
  A <- 1
  B <- m
  C <- -s^2
  
  disc <- B^2 - 4*A*C

  b_val <- (-B + sqrt(disc)) / (2 * A)
  
  a_val <- m / b_val + 1
  
  return(list(a = a_val, b = b_val, mode = m))
}


log_delta_pdf <- function(delta, a, b, inverse = FALSE) {
  #' This function calculates the log-PDF of Delta based on the Gamma distribution.
  #' It allows for an inverse.
  #'
  #' Parameters:
  #' - delta (numeric): The gamma distributed value at which to evaluate the PDF.
  #' - a (numeric): The shape parameter of the Gamma distribution.
  #' - b (numeric): The scale parameter of the Gamma distribution.
  #' - inverse (logical, default = FALSE): If TRUE, computes the log-PDF for the inverse of Delta.
  #'
  #' Returns:
  #' - val (numeric): The logarithm of the PDF evaluated at the specified delta. If inverse is TRUE,
  #'   it returns the log-PDF of the inverse of delta.
  
  if(inverse){
    val <- a * log(b) - (a + 1) * log(delta) - b / delta - lgamma(a)
  } else{
    val <- dgamma(delta, shape = a, scale = b, log = TRUE)
  }
  return(val)
}


log_ml_Y_density <- function(X, M, P, Q, v) {
  #' Computes the logarithm of the marginal likelihood density for Y.
  #'
  #' Parameters:
  #' - X (matrix): The observed data matrix of dimensions p x q.
  #' - M (matrix): The mean matrix of dimensions p x q, representing the expected values.
  #' - P (matrix): The precision matrix of dimensions p x p, associated with the distribution of X.
  #' - Q (matrix): The precision matrix of dimensions q x q, associated with the distribution of X.
  #' - v (numeric): The degrees of freedom parameter, influencing the shape of the distribution.
  #'
  #' Returns:
  #' - log_pdf (numeric): The logarithm of the marginal likelihood density evaluated at the provided X.

  p <- nrow(X) # T  
  q <- ncol(X) # k
  
  XM <- X - M

  PXM <- P %*% XM
  S_mat <- Q + t(XM) %*% PXM
  
  # log-determinants
  chol_P <- chol(P)
  logdetP <- 2 * sum(log(diag(chol_P)))
  
  chol_Q <- chol(Q)
  logdetQ <- 2 * sum(log(diag(chol_Q)))
  
  chol_S <- chol(S_mat)
  logdetS <- 2 * sum(log(diag(chol_S)))
  
  # sum of gamma terms
  i <- 1:q
  gamma_terms <- lgamma((v+1-i)/2) - lgamma((v+p+1-i)/2)
  gamma_sum <- sum(gamma_terms)
  
  # put it all together
  log_pdf <- - (p*q/2)*log(pi) + (q/2)*logdetP + (v/2)*logdetQ - gamma_sum - ((v + p)/2)*logdetS
  
  return(log_pdf)
}


log_posterior_delta <- function(Y, X, M0, V0, S0, v0, pi1_val, mu_val, gamma_val, 
                                s2_diag_val, hyper_params, upper_bound, lower_bound){
  #' This function calculates the log-posterior probability of the hyper parameters
  #' (pi1_val, mu_val, gamma_val, s2_diag_val) given the observed data Y,
  #' predictor matrix X, and prior parameters. It evaluates the marginal likelihood
  #' of Y and incorporates the prior distributions of the hyper parameters to obtain
  #' the posterior density. The function also enforces parameter bounds, returning a
  #' large negative value if any parameter falls outside the specified limits.
  #'
  #' Parameters:
  #' - Y (matrix): The observed data matrix of dimensions T_eff x k.
  #' - X (matrix): The predictor matrix of dimensions T_eff x m.
  #' - M0 (matrix): The prior mean matrix of dimensions compatible with X %*% M0.
  #' - V0 (matrix): The prior covariance matrix used in calculating P_term, of dimensions compatible with X.
  #' - S0 (matrix): The prior precision matrix of dimensions q x q, used in the marginal likelihood calculation.
  #' - v0 (numeric): The prior degrees of freedom parameter for the marginal likelihood.
  #' - pi1_val (numeric): The value of the hyper parameter π₁ for which the posterior density is evaluated.
  #' - mu_val (numeric): The value of the hyper parameter μ for which the posterior density is evaluated.
  #' - gamma_val (numeric): The value of the hyper parameter γ for which the posterior density is evaluated.
  #' - s2_diag_val (numeric vector): The values of the hyper parameters σ²₁, σ²₂, ..., σ²ₖ for which the posterior density is evaluated.
  #' - hyper_params (list): A list containing hyperparameters for the prior distributions of hyper parameters, structured as:
  #'   - pi1_param (list): Contains a and b for π₁'s Gamma prior.
  #'   - mu_param (list): Contains a and b for μ's Gamma prior.
  #'   - gamma_param (list): Contains a and b for γ's Gamma prior.
  #'   - s2_diag_param (list): Contains a and b for each σ²ᵢ's Gamma prior.
  #' - upper_bound (numeric vector): The upper bounds for the hyper parameters. Must be the same length as hyper.
  #' - lower_bound (numeric vector): The lower bounds for the hyper parameters. Must be the same length as hyper.
  #'
  #' Returns:
  #' - log_posterior (numeric): The logarithm of the posterior density evaluated at the specified hyper parameters.
  #'   Returns -1e18 if any hyper parameter is outside the specified bounds, effectively assigning zero probability.
  
  hyper <- c(pi1_val, mu_val, gamma_val, s2_diag_val)
  
  if(any(lower_bound > hyper | hyper > upper_bound)){
    return(-1e18)
  }
  
  T_eff <- nrow(Y)
  k <- ncol(Y)
  
  # --- log(m(Y|delta))
  M_term <- X %*% M0
  
  # P_term
  temp <- diag(T_eff) + X %*% V0 %*% t(X)
  P_term <- chol2inv(chol(temp))
  
  # log marginal likelihood
  log_ml <- log_ml_Y_density(Y, M_term, P_term, S0, v0)
  
  # --- log(p(delta))
  # extract hyperparameters
  pi1_a <- hyper_params[["pi1_param"]]$a
  pi1_b <- hyper_params[["pi1_param"]]$b
  
  mu_a <- hyper_params[["mu_param"]]$a
  mu_b <- hyper_params[["mu_param"]]$b
  
  gamma_a <- hyper_params[["gamma_param"]]$a
  gamma_b <- hyper_params[["gamma_param"]]$b
  
  s2_diag_a <- hyper_params[["s2_diag_param"]]$a
  s2_diag_b <- hyper_params[["s2_diag_param"]]$b
  
  # compute log PDFs
  log_pi1_pdf <- log_delta_pdf(pi1_val, pi1_a, pi1_b)
  log_mu_pdf <- log_delta_pdf(mu_val, mu_a, mu_b)
  log_gamma_pdf <- log_delta_pdf(gamma_val, gamma_a, gamma_b)
  log_s2_diag_pdf <- log_delta_pdf(s2_diag_val, s2_diag_a, s2_diag_b, inverse = TRUE)

  return(log_ml + log_pi1_pdf + log_mu_pdf + log_gamma_pdf + sum(log_s2_diag_pdf))
}


optim_log_post_delta <- function(par, Y, X, p, intercept, pi4, hyper_params, lower_bound, upper_bound) {
  #' Computes the logarithm of the posterior density for hyper parameters used 
  #' only in the optimization step (optim()).
  #'
  #' Parameters:
  #' - par (numeric vector): A vector of hyper parameters.
  #' - Y (matrix): The observed data matrix of dimensions T_eff x k, where T_eff is the number of observations and k is the number of variables.
  #' - X (matrix): The predictor matrix of dimensions T_eff x n, where n is the number of predictors.
  #' - p (integer): The lag order or another relevant integer parameter used in the hierarchical_prior function.
  #' - intercept (logical): Indicates whether an intercept term is included in the model.
  #' - pi4 (numeric): A specific hyper parameter value used in the hierarchical_prior function.
  #' - hyper_params (list): A list containing hyperparameters for the prior distributions of hyper parameters.
  #' - lower_bound (numeric vector): The lower bounds for the hyper parameters (pi1_val, mu_val, gamma_val, and each s2_diag_val).
  #' - upper_bound (numeric vector): The upper bounds for the hyper parameters (pi1_val, mu_val, gamma_val, and each s2_diag_val).
  #'
  #' Returns:
  #' - lp (numeric): The logarithm of the posterior density evaluated at the specified hyper parameters. 
  #'   Returns -1e18 if any parameter in par is outside the specified bounds.
  
  k <- ncol(Y)
  
  pi1_val <- par[1]
  mu_val <- par[2]
  gamma_val <- par[3]
  s2_diag_val <- par[4:(3 + k)]
  
  # compute priors with hyper parameters
  prior_new <- hierarchical_prior(Y, p, intercept, 
                                  lag_mean = 1, 
                                  s2_diag = s2_diag_val, 
                                  pi1 = pi1_val, 
                                  pi4 = pi4)
  M0_new <- prior_new$M0
  V0_new <- prior_new$V0
  S0_new <- prior_new$S0
  v0_new <- prior_new$v0
  dmy_new <- create_dummy_observations(Y, p, intercept, mu = mu_val, gamma = gamma_val)
  Y_aug_new <- rbind(dmy_new$Y_dum, Y)
  X_aug_new <- rbind(dmy_new$X_dum, X)
  
  # call original log_posterior_delta
  lp <- log_posterior_delta(
    Y = Y_aug_new, X = X_aug_new, M0 = M0_new, V0 = V0_new, S0 = S0_new, v0 = v0_new,
    pi1_val = pi1_val, mu_val = mu_val, gamma_val = gamma_val,
    s2_diag_val = s2_diag_val,
    hyper_params = hyper_params, lower_bound = lower_bound, upper_bound = upper_bound
  )
  
  return(lp)
}
