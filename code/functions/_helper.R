# data adjusting functions for bvar analysis


create_lagged_data <- function(Y, p, intercept = FALSE,
                                            use_dummies = FALSE,
                                            dummy_pars = list(delta = 1.0, gamma = 1.0, prior_mean = NULL)) {
  #' Creates the response matrix Y and design matrix X for a VAR model, incorporating lagged observations from the real data. 
  #' If specified, it also augments the matrices with dummy observations.
  #' 
  #' Parameters:
  #' - Y (matrix): The dataset (T x k), where T is the number of time points and k is the number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - intercept (logical): Whether to include an intercept term in the design matrix (default = FALSE).
  #' - use_dummies (logical): Whether to include dummy observations in the augmented matrices (default = FALSE).
  #' - dummy_pars (list): Parameters for the dummy observations, including:
  #'   - delta (numeric): Tightness for the initial observational prior (default = 1.0).
  #'   - gamma (numeric): Tightness for the sum-of-coefficients prior (default = 1.0).
  #'   - prior_mean (numeric vector): Mean vector for the prior. If NULL, it is derived from Yraw (default = NULL).
  #' 
  #' Returns:
  #' - A list containing:
  #'   - Y (matrix): The augmented response matrix, including dummy observations if use_dummies is TRUE.
  #'   - X (matrix): The augmented design matrix, including dummy observations if use_dummies is TRUE.
  
  # real-data-lagged
  ld <- helper_lagged_data(Y, p, intercept = intercept)
  Y_real <- ld$y
  X_real <- ld$X
  
  if(!use_dummies) {
    return(list(Y = Y_real, X = X_real))
  } else {
    # build dummy blocks
    dums <- do.call(create_dummy_observations,
                    args = c(list(Y=Y, p=p, intercept=intercept),
                             dummy_pars))
    Y_dum <- dums$Y_dum
    X_dum <- dums$X_dum
    
    # append
    Y_aug <- rbind(Y_dum, Y_real)
    X_aug <- rbind(X_dum, X_real)
    
    return(list(Y = Y_aug, X = X_aug))
  }
}



helper_lagged_data <- function(data, p, intercept = FALSE) {
  #' Generates the response matrix y and the design matrix X from a given dataset, 
  #' preparing the data for VAR model estimation. The design matrix can optionally include an intercept term.
  #' 
  #' Parameters:
  #' - data (matrix): A T x k matrix of observations, where T is the number of time points and k is the number of variables.
  #' - p (integer): Number of lags to include in the design matrix.
  #' - intercept (logical): Whether to include a column of ones as an intercept in the design matrix (default = FALSE).
  #' 
  #' Returns:
  #' - A list with two elements:
  #'   - y: A (T-p) x k matrix, where rows correspond to observations and columns to variables.
  #'   - X: A (T-p) x (k*p) matrix (or (T-p) x (1 + k*p) if intercept is TRUE), representing the design 
  #'   matrix with lagged values and optionally an intercept term.
  
  T <- nrow(data)
  k <- ncol(data)
  
  lagged_data <- embed(data, p + 1)
  
  # current Y_t is always the first k columns
  y <- lagged_data[, 1:k, drop=FALSE]
  lag_cols <- lagged_data[, (k + 1):ncol(lagged_data), drop=FALSE]
  
  if (intercept) {
    X <- cbind(1, lag_cols)
  } else {
    X <- lag_cols
  }
  
  return(list(y = y, X = X))
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
