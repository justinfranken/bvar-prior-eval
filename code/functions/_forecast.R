# functions to predict from coefficient draws and measure errors.


bvar_forecast <- function(
    Y, 
    phi_draws, 
    Sigma_draws, 
    p, 
    h, 
    alpha = 0.05, 
    intercept = TRUE,
    draw_shocks = TRUE
) {
  #' Generates forecasts from a Bayesian VAR model based on posterior draws.
  #'
  #' Parameters:
  #' - Y (matrix): A T x k data matrix. Only the last p rows are used as initial conditions.
  #' - phi_draws (array): An array of posterior draws of the VAR coefficients, with dimensions 
  #'   ((k*p)+1) x k x M, where the first row represents the intercept.
  #' - Sigma_draws (array): An array of posterior draws of the error covariance matrix, with 
  #'   dimensions k x k x M.
  #' - p (integer): The number of lags in the VAR model.
  #' - h (integer): The forecast horizon, i.e., the number of periods ahead to forecast.
  #' - alpha (numeric, default = 0.05): The significance level used to compute forecast 
  #'   prediction intervals (e.g., a 95% interval when alpha = 0.05).
  #' - intercept (logical, default = TRUE): If TRUE, include an intercept in the 
  #'   forecast. If FALSE, the intercept is set to zero and the first row of phi_draws is 
  #'   ignored in the calculation.
  #' - draw_shocks: If FALSE, no shocks will be drawn. If TRUE, shocks will be drawn.
  #'
  #' Returns:
  #' - A list containing:
  #'   - median: An h x k matrix of median forecasts for each horizon and variable.
  #'   - lower: An h x k matrix of lower quantile forecasts (alpha/2 quantile).
  #'   - upper: An h x k matrix of upper quantile forecasts (1 - alpha/2 quantile).
  #'   - draws: An array of all forecast draws with dimensions h x k x M.
  
  Tfull <- nrow(Y)
  k     <- ncol(Y)
  M     <- dim(phi_draws)[3]
  
  # Storage for forecast draws: h steps, k variables, M draws
  forecast_draws <- array(NA, dim = c(h, k, M))
  
  # Extract the last p observations from the data
  Y_init <- Y[(Tfull - p + 1):Tfull, , drop = FALSE]
  
  # Define the offset and how we handle the intercept
  offset <- if (intercept) 1 else 0
  
  for (j in 1:M) {
    
    # current draws
    Phi_j   <- phi_draws[ , , j]  # ( (k*p) + 1 ) x k
    Sigma_j <- Sigma_draws[ , , j]
    chol_Sig_j <- chol(Sigma_j)
    
    # Prepare matrix to hold forecasts for draw j
    Y_hat_j <- matrix(0, nrow = h, ncol = k)
    current_stack <- rbind(Y_init, matrix(NA, nrow = h, ncol = k))
    
    # forecast recursively
    for (step in 1:h) {
      
      # Intercept: either take the first row of Phi_j or set to zero
      if (intercept) {
        inter <- Phi_j[1, ]
      } else {
        inter <- rep(0, k)
      }
      
      # Summation of p lags
      lag_contrib <- rep(0, k)
      for (lag_i in 1:p) {
        # Identify the rows in Phi_j that correspond to the lag_i-th lag
        # If we have an intercept row, we shift by 1
        row_start <- offset + (lag_i - 1) * k + 1
        row_end   <- row_start + k - 1
        
        Phi_lag_i <- Phi_j[row_start:row_end, ]
        
        # multiply by the appropriate row in current_stack
        lag_contrib <- lag_contrib + Phi_lag_i %*% current_stack[p + step - lag_i, ]
      }
      
      # random errors
      eps_j <- crossprod(chol_Sig_j, rnorm(k))
      if(!draw_shocks) eps_j = rep(0, k)
      
      # new forecast
      Y_new <- inter + lag_contrib + eps_j
      current_stack[p + step, ] <- Y_new
      
      # save it
      Y_hat_j[step, ] <- Y_new
    }
    
    # store forecasts for this draw
    forecast_draws[, , j] <- Y_hat_j
  }
  
  # compute median and prediction quantiles
  fc_median <- apply(forecast_draws, c(1, 2), median)
  fc_lower  <- apply(forecast_draws, c(1, 2), quantile, probs = alpha/2)
  fc_upper  <- apply(forecast_draws, c(1, 2), quantile, probs = 1 - alpha/2)
  
  return(list(
    median  = fc_median,
    lower   = fc_lower,
    upper   = fc_upper,
    draws   = forecast_draws
  ))
}


fcst_rmse <- function(fcst_obj, Ypred, h, k){
  #' Computes the root mean squared error (RMSE) for forecasts.
  #'
  #' Parameters:
  #' - fcst_obj (list): A forecast object from bvar_forecast().
  #' - Ypred (matrix): A matrix of true values for the forecast period, with dimensions h x k.
  #' - h (integer): The forecast horizon (number of steps ahead).
  #' - k (integer): The number of variables being forecast.
  #'
  #' Returns:
  #' - A list containing:
  #'   - row_rmse: An h x 1 matrix of RMSE values computed for each forecast horizon (across variables).
  #'   - col_rmse: A 1 x k matrix of RMSE values computed for each variable (across forecast horizons).
  #'   - all_rmse: A single numeric value representing the overall RMSE across all forecasts.
  
  # h-step ahead forecast
  row_rmse_obj <- matrix(rep(NA, h), ncol = 1, dimnames = list(paste0("T+", 1:h), "Y1 ... Yk"))
  for (i in 1:h) {
    row_rmse_obj[i,] <- sqrt(mean((fcst_obj$median[i,] - Ypred[i,])^2))
  }
  
  # variable forecast
  col_rmse_obj <- matrix(rep(NA, k), nrow = 1, dimnames = list("T1 ... Th", paste0("Y", 1:k)))
  for (i in 1:k) {
    col_rmse_obj[,i] <- sqrt(mean((fcst_obj$median[,i] - Ypred[,i])^2))
  }
  
  # total rmse
  all_rmse_obj <- sqrt(mean((fcst_obj$median - Ypred)^2))
  
  return(list(
    row_rmse = row_rmse_obj,
    col_rmse = col_rmse_obj,
    all_rmse = all_rmse_obj
  ))
}


fcst_pred_int_acc <- function(fcst_obj, Ypred, h){
  #' Evaluates the forecast prediction interval accuracy.
  #'
  #' Parameters:
  #' - fcst_obj (list): A forecast object from bvar_forecast().
  #' - Ypred (matrix): A matrix of true values for the forecast period, with dimensions h x k.
  #' - h (integer): The forecast horizon (number of steps ahead).
  #'
  #' Returns:
  #' - A matrix of dimensions h x 1 where each entry is 1 if the true value for that forecast horizon 
  #'   falls within the prediction interval for at least one variable, and 0 otherwise.
  
  result <- matrix(0, nrow = h, ncol = 1, dimnames = list(paste0("T+", 1:h), "Y1 ... Yk"))
  
  for (i in 1:h) {
    temp <- (fcst_obj$lower[i,] <= Ypred[i,]) & (Ypred[i,] <= fcst_obj$upper[i,])
    
    # if in pred interval 1, else 0
    result[i,] <- ifelse(any(temp), 1, 0)
  }
  
  return(result)
}
