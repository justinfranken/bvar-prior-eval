# functions to predict coefficient draws from raw data


predict_bvar <- function(Phi_store, Sigma_store, 
                         Yraw, p, H = 4, 
                         draw_shocks = TRUE,
                         intercept = FALSE,
                         n_cores = 1) {
  #' This function produces H-step ahead forecasts for a Bayesian VAR model using posterior draws of 
  #' the coefficient matrix (Phi_store) and covariance matrix (Sigma_store). Forecasts are computed 
  #' for each posterior draw, given the last p observations from the dataset.
  #' 
  #' Parameters:
  #' - Phi_store (array): A 3D array of VAR coefficients with dimensions [p*k, k, n_post_draws], 
  #' where n_post_draws is the number of posterior draws.
  #' - Sigma_store (array): A 3D array of covariance matrices with dimensions [k, k, n_post_draws].
  #' - Yraw (matrix): The full dataset (T x k), where T is the number of observations, and k is the 
  #' number of variables.
  #' - p (integer): The lag order of the VAR model.
  #' - H (integer): Forecast horizon, i.e., how many steps ahead to predict (default = 4).
  #' - draw_shocks (logical): Whether to sample random shocks from the covariance matrix for each 
  #' forecast step (default = TRUE).
  #' - intercept (logical): Whether the Phi_store includes an intercept term (default = FALSE).
  #' - n_cores (integer): Number of parallel R sessions (processes) that will be used to execute prediction.
  #' 
  #' Returns:
  #' - A 3D array of forecasts with dimensions [H, k, n_post_draws], where H is the forecast horizon, 
  #' k is the number of variables, and n_post_draws is the number of posterior draws.

  k <- ncol(Yraw)
  T_data <- nrow(Yraw)
  Y_last_p <- Yraw[(T_data - p + 1):T_data, , drop = FALSE]
  n_post_draws <- dim(Phi_store)[3]
  
  # initiate future.apply 
  plan(multisession, workers = n_cores)
  
  forecast_list <- future_lapply(seq_len(n_post_draws), function(i) {
    Phi_draw_i   <- Phi_store[ , , i]
    Sigma_draw_i <- Sigma_store[ , , i]
    
    helper_forecast_from_draw(
      Phi_draw_init = Phi_draw_i,
      Sigma_draw    = Sigma_draw_i,
      Y_last_p      = Y_last_p,
      H             = H,
      draw_shocks   = draw_shocks,
      intercept     = intercept
    )
  })
  # reset to standard environment
  plan(sequential)  
  
  Y_forecasts <- array(
    unlist(forecast_list),
    dim = c(H, k, n_post_draws),
    dimnames = list(paste0("T+", 1:H), colnames(Yraw), NULL)
  )
  return(Y_forecasts)
}



helper_forecast_from_draw <- function(Phi_draw_init, Sigma_draw,
                               Y_last_p,
                               H = 4,
                               draw_shocks = TRUE,
                               intercept = FALSE) {
  #' This function produces H-step ahead forecasts based on a VAR(p) model, given the last p 
  #' observations, a draw of the coefficient matrix (Phi_draw), and the covariance matrix (Sigma_draw). 
  #' Optionally, it can include an intercept and draw random shocks to incorporate predictive uncertainty.
  #' 
  #' Parameters:
  #' - Phi_draw_init (matrix): A (p*k) x k matrix of VAR coefficients, optionally including an intercept as the first row.
  #' - Sigma_draw (matrix): A k x k covariance matrix for the shocks.
  #' - Y_last_p (matrix): A p x k matrix of the last p observations, ordered in ascending time (row 1 is the oldest, row p is the most recent).
  #' - H (integer): Forecast horizon, i.e., how many steps ahead to predict (default = 4).
  #' - draw_shocks (logical): Whether to include random shocks in the forecast (default = TRUE).
  #' - intercept (logical): Whether Phi_draw_init includes an intercept in the first row (default = FALSE).
  #' 
  #' Returns:
  #' - A H x k matrix of forecasts, with rows representing time steps and columns representing variables.
  
  k <- ncol(Y_last_p)
  p <- nrow(Phi_draw_init) / k
  
  Phi_intercept <- rep(0, k)
  if (intercept) {
    Phi_intercept <- Phi_draw_init[1, ]
    Phi_draw      <- Phi_draw_init[-1, ]
  } else {
    Phi_draw      <- Phi_draw_init
  }
  
  Phi_list <- hhelper_reshape_phi(Phi_draw, k, p)
  
  Y_future <- matrix(NA, nrow = H, ncol = k,
                     dimnames = list(paste0("T+", 1:H), colnames(Y_last_p)))
  
  # Start buffer from last p observed data
  Y_buffer <- Y_last_p
  
  for (h in seq_len(H)) {
    # sum_{ell=1 to p} (Phi_list[[ell]] %*% Y_buffer[p-(ell-1), ])
    Y_next_mean <- Reduce(
      f = `+`,
      x = lapply(seq_len(p), function(ell) {
        Phi_list[[ell]] %*% Y_buffer[p - (ell - 1), ]
      })
    )
    
    if (draw_shocks) {
      e_t <- mvrnorm(1, mu = rep(0, k), Sigma = Sigma_draw)
    } else {
      e_t <- rep(0, k)
    }
    
    Y_next <- as.numeric(Phi_intercept + Y_next_mean + e_t)
    Y_future[h, ] <- Y_next
    
    # Update buffer
    Y_buffer <- rbind(Y_buffer[-1, ], Y_next)
  }
  return(Y_future)
}


hhelper_reshape_phi <- function(Phi_draw, k, p) {
  #' This function takes a matrix with dimensions (p*k) x k and reshapes it into a list of p 
  #' matrices, each with dimensions k x k. Each matrix corresponds to a specific lag in a VAR model.
  #' 
  #' Parameters:
  #' - Phi_draw (matrix): A (p*k) x k matrix to be reshaped.
  #' - k (integer): Number of variables in the VAR model.
  #' - p (integer): Number of lags in the VAR model.
  #' 
  #' Returns:
  #' - A list of p matrices, where each matrix has dimensions k x k.
  
  Phi_list <- vector("list", length = p)
  for (ell in seq_len(p)) {
    row_start <- (ell - 1) * k + 1
    row_end   <- ell * k
    Phi_list[[ell]] <- Phi_draw[row_start:row_end, , drop = FALSE]
  }
  return(Phi_list)
}
