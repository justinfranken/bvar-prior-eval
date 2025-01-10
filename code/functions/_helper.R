# helper functions in bvar analysis


create_lagged_data <- function(Y, p, intercept = FALSE,
                                            use_dummies = FALSE,
                                            dummy_pars = list()) {
  #' Creates the response matrix Y and design matrix X for a VAR model, incorporating lagged observations from the real data. 
  #' If specified, it also augments the matrices with dummy observations.
  #' 
  #' Parameters:
  #' - Y (matrix): The dataset (T x k), where T is the number of time points and k is the number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - intercept (logical): Whether to include an intercept term in the design matrix (default = FALSE).
  #' - use_dummies (logical): Whether to include dummy observations in the augmented matrices (default = FALSE).
  #' - dummy_pars (list): Parameters for the dummy observations, passed to create_dummy_observations (default = empty list).
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
