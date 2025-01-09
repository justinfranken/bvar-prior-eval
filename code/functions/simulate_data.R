# simulating financial growth rates for Monte Carlo Simulations


sim_y <- function(K,
                  p, 
                  T,
                  # --- PARAMETERS FOR TIME SERIES SHOCKS ----------------------
                  init_sd = 0.015,
                  # --- individual shocks -----------------------
                  min_indiv_shocks = 0,
                  max_indiv_shocks = 0,
                  indiv_shock_length_min = 1,
                  indiv_shock_length_max = 3,
                  # --- high volatility period ------------------
                  min_high_vol_periods = 0,
                  max_high_vol_periods = 0,
                  high_vol_period_length_min = 4,
                  high_vol_period_length_max = 8,
                  # --- exogenous shocks ------------------------
                  min_exog_shocks = 0,
                  max_exog_shocks = 0,
                  exog_shock_length_min = 1,
                  exog_shock_length_max = 2){
  #' Creates T simulations of a K*1 vector y of a VAR(p) process.
  #' 
  #' Parameters:
  #' - K (integer): How many variables are included.
  #' - p (integer): How many lags are included.
  #' - T (integer): How many observations are going to be simulated.
  #' - init_sd (scalar): Initial volatility for first p lags.
  #' - min_..._shocks, max_..._shocks (integer): Range of how many individual/ exogenous
  #'                                    shocks can appear in each variable's time series.
  #' - ..._length_min, ..._length_max (scalar): Range of durations for the shocks.
  #' 
  #' Returns:
  #' - K*1 matrix of simulated y's.

  # compute coefficient matrix A
  A_list <- helper_coefficient_generator(p, K)

  # compute shocks
  eps <- helper_eps_generator(K = K, 
                              T = T+p,
                              min_indiv_shocks = min_indiv_shocks,
                              max_indiv_shocks = max_indiv_shocks,
                              indiv_shock_length_min = indiv_shock_length_min,
                              indiv_shock_length_max = indiv_shock_length_max,
                              min_high_vol_periods = min_high_vol_periods,
                              max_high_vol_periods = max_high_vol_periods,
                              high_vol_period_length_min = high_vol_period_length_min,
                              high_vol_period_length_max = high_vol_period_length_max,
                              min_exog_shocks = min_exog_shocks,
                              max_exog_shocks = max_exog_shocks,
                              exog_shock_length_min = exog_shock_length_min,
                              exog_shock_length_max = exog_shock_length_max)
  
  # initialize series
  Y <- matrix(0, nrow = T + p, ncol = K)
  for (i in 1:p) {
    Y[i, ] <- rnorm(K, mean = 0, sd = init_sd)
  }
  
  # simulate series
  for (t in (p+1):(T + p)) {
    Y_t <- rep(0, K)
    for (lag in 1:p) {
      Y_t <- Y_t + A_list[[lag]] %*% Y[t-lag, ]
    }
    Y[t, ] <- Y_t + eps[t, ]
  }
  
  return(Y[(p+1):(T + p), ])
}


helper_coefficient_generator <- function(p, 
                                         K, 
                                         d_min = 0.1,
                                         d_max = 0.3, 
                                         off_d_mean = 0, 
                                         off_d_sd = 0.05,
                                         persistence_factor = 0.4) {
  #' Helper function which creates a list of p many K*K VAR coefficient 
  #' matrices, which decay lineally with increasing lags p.
  #' 
  #' Parameters:
  #' - p (integer): How many lags are included.
  #' - K (integer): How many variables are included.
  #' - d_min_coef, d_max_coef (scalar): Range of the diagonal elements in the coefficient matrix.
  #' - off_d_mean_coef, off_d_sd_coef (scalar): Mean and standard deviation of off-diagonal
  #'                                            elements in the coefficient matrix.
  #' - persistence_factor (scalar): Determines how strongly coefficients of lag p+1
  #'                                 depend on those of lag p (default = 0.5).
  #' 
  #' Returns:
  #' - A list of p-many K*K VAR coefficient matrices.
  
  A_list <- vector("list", p)
  prev_A <- NULL
  
  for (i in 1:p) {
    
    diag_min <- d_min / i
    diag_max <- d_max / i
    off_diag_mean <- off_d_mean / i
    off_diag_sd   <- off_d_sd / i
    
    A_raw <- matrix(0, nrow=K, ncol=K)

    # Generate elements with persistence
    for (r in 1:K) {
      for (col in 1:K) {
        if (r != col) {
          if (!is.null(prev_A)) {
            prev_sign <- sign(prev_A[r, col])
            base_value <- rnorm(1, mean = persistence_factor * prev_A[r, col], sd = off_diag_sd)
            A_raw[r, col] <- base_value + prev_sign * abs(rnorm(1, mean = off_diag_mean, sd = off_diag_sd))
          } else {
            # Initial lag coefficients sampled independently
            A_raw[r, col] <- rnorm(1, mean = off_diag_mean, sd = off_diag_sd)
          }
        }
        else {
          if(!is.null(prev_A)){
            prev_sign <- sign(prev_A[r, col])
            base_value <- runif(1, min = persistence_factor * diag_min, max = persistence_factor * diag_max)
            A_raw[r, col] <- base_value + prev_sign * abs(runif(1, min = diag_min, max = diag_max))
          }
          else{
            # Initial lag coefficients sampled independently
            diag_sign <- sample(c(-1,1), 1, prob = c(0.90, 0.10))
            A_raw[r, col] <- runif(1, min = diag_min, max = diag_max) * diag_sign
          }
        }
      }
    }
    
    # Update the list and store the current matrix as previous
    A_list[[i]] <- A_raw
    prev_A <- A_raw
  }
  
  # Check if stable
  A_list <- hhelper_make_stationary_A(A_list, K)
  
  return(A_list)
}



helper_eps_generator <- function(K, 
                                 T, 
                                 shock_diag_min = 0.05,
                                 shock_diag_max = 0.2, 
                                 mean_vola = 0.025, 
                                 sd_vola = 0.005,
                                 # --- individual shocks -----------------------
                                 min_indiv_shocks,
                                 max_indiv_shocks,
                                 indiv_shock_length_min,
                                 indiv_shock_length_max,
                                 indiv_shock_ampl_min = 0.075,
                                 indiv_shock_ampl_max = 0.2,
                                 indiv_shock_decay = 0.5,
                                 # --- high volatility period ------------------
                                 min_high_vol_periods,
                                 max_high_vol_periods,
                                 high_vol_period_length_min,
                                 high_vol_period_length_max,
                                 high_vol_strength_min = 2,
                                 high_vol_strength_max = 3.5,
                                 # --- exogenous shocks ------------------------
                                 min_exog_shocks,
                                 max_exog_shocks,
                                 exog_shock_length_min,
                                 exog_shock_length_max,
                                 exog_shock_ampl_min = 0.075,
                                 exog_shock_ampl_max = 0.15,
                                 exog_shock_decay = 0.5,
                                 exog_shock_shift_range = 3) {
  #' Helper function to creates a T*K matrix of shocks. First, draws from a 
  #' multivariate normal distribution with a correlation structure Rho. Can add
  #' high vola periods and exogenous and individual shocks, where exogenous shocks
  #' are correlated with high vola periods if they exist.
  #'
  #' Parameters:
  #' - K (integer): Number of variables.
  #' - T (integer): Number of observations.
  #' - shock_diag_min, shock_diag_max (scalar): Range for the off-diagonal correlations.
  #' - mean_vola, sd_vola (scalar): Mean and SD for drawing each firm's volatility.
  #' --- shock parameters ------------------------------------------------------
  #' - min_indiv_shocks, max_indiv_shocks, in_indiv_shocks, 
  #'   max_indiv_shocks, in_indiv_shocks, max_indiv_shocks (integer): 
  #'    - Range of how many "shocks" can appear in each variable's time series.
  #' - indiv_shock_length_min, indiv_shock_length_max, exog_shock_length_min, 
  #'   exog_shock_length_max, high_vol_period_length_min, 
  #'   high_vol_period_length_max (integer): 
  #'    - Range of durations for these shocks.
  #' - indiv_shock_ampl_min, indiv_shock_ampl_max, exog_shock_ampl_min, 
  #'   exog_shock_ampl_max (scalar): 
  #'    - Range for the amplitude of each shock.
  #' - high_vol_strength_min, high_vol_strength_max (scalar):
  #'    - Range for the intensity of increased volatility.
  #' - indiv_shock_decay, exog_shock_decay (scalar): Decay factor for the shock.
  #' - exog_shock_shift_range (integer): how many periods around the high-vol 
  #'   segment we shift exogenous shocks.
  #'
  #' Returns:
  #' - A T*K matrix containing the simulation shocks.
  
  # create baseline correlation matrix
  Rho <- matrix(0, nrow = K, ncol = K)
  diag(Rho) <- 1
  for (r in 1:(K-1)) {
    for (col in (r+1):K) {
      x <- runif(1, min = shock_diag_min, max = shock_diag_max)
      Rho[r, col] <- x
      Rho[col, r] <- x
    }
  }

  # construct the covariance matrix Sigma
  sd_vec <- rnorm(K, mean = mean_vola, sd = sd_vola)
  Sigma  <- diag(sd_vec) %*% Rho %*% diag(sd_vec)
  
  # generate the base normal random vectors with covariance Sigma
  cholSigma <- chol(Sigma)
  Z         <- matrix(rnorm(T * K), nrow = T, ncol = K)
  eps       <- Z %*% cholSigma
  
  
  # --- PERIODS OF HIGHER VOLATILITY -------------------------------------------
  number_high_vol_periods <- sample(min_high_vol_periods:max_high_vol_periods, 1)
  
  hv_periods <- data.frame(
    hv_start  = integer(number_high_vol_periods),
    hv_length = integer(number_high_vol_periods)
  )
  
  if (number_high_vol_periods > 0) {
    for (i in 1:number_high_vol_periods) {
      hv_length <- sample(high_vol_period_length_min:high_vol_period_length_max, 1)
      hv_start  <- sample(seq_len(T - hv_length), 1)
      hv_multiplier <- runif(1, min = high_vol_strength_min, max = high_vol_strength_max)
      
      range_idx <- hv_start:(hv_start + hv_length - 1)
      eps[range_idx, ] <- eps[range_idx, ] * hv_multiplier
      
      # store these in hv_periods so we can reuse them
      hv_periods$hv_start[i]  <- hv_start
      hv_periods$hv_length[i] <- hv_length
    }
  }
  
  
  # --- INDIVIDUAL SHOCKS ------------------------------------------------------
  num_indiv_shocks_vec <- sample(min_indiv_shocks:max_indiv_shocks, K, replace = TRUE)
  
  for (k in 1:K) {
    n_indiv_shocks <- num_indiv_shocks_vec[k]
    
    if (n_indiv_shocks > 0) {
      for (s in 1:n_indiv_shocks) {
        indiv_shock_length  <- sample(indiv_shock_length_min:indiv_shock_length_max, 1)
        indiv_shock_start   <- sample(seq_len(T - indiv_shock_length), 1)
        
        indiv_shock_amplitude <- runif(1, indiv_shock_ampl_min, indiv_shock_ampl_max)
        indiv_shock_sign      <- sample(c(-1, 1), 1)
        indiv_shock_amplitude <- indiv_shock_amplitude * indiv_shock_sign
        
        for (i in 0:(indiv_shock_length - 1)) {
          t_index <- indiv_shock_start + i
          eps[t_index, k] <- eps[t_index, k] + indiv_shock_amplitude * (indiv_shock_decay^i)
        }
      }
    }
  }
  
  
  # --- EXOGENOUS SHOCKS -------------------------------------------------------
  n_exog_shocks <- sample(min_exog_shocks:max_exog_shocks, 1)
  
  if (n_exog_shocks > 0) {
    for (s in 1:n_exog_shocks) {
      exog_shock_length    <- sample(exog_shock_length_min:exog_shock_length_max, 1)
      exog_shock_amplitude <- runif(K, exog_shock_ampl_min, exog_shock_ampl_max)
      
      # pick a random vola period to center around
      if (number_high_vol_periods > 0) {
        idx <- sample(seq_len(nrow(hv_periods)), 1)
        center <- hv_periods$hv_start[idx]
        
        # shift the exog shock start
        offset <- sample(seq(-exog_shock_shift_range, exog_shock_shift_range), 1)
        exog_start <- center + offset
      } else {
        # if no vola periods, just pick random start
        exog_start <- sample(seq_len(T - exog_shock_length), 1)
      }
      
      # safety check
      if (exog_start < 1) exog_start <- 1
      if (exog_start + exog_shock_length - 1 > T) {
        exog_shock_length <- T - exog_start + 1
      }
      
      # sign can differ for each variable
      sign_vec <- sample(c(-1, 1), 1)
      
      # apply to all K variables
      for (i in 0:(exog_shock_length - 1)) {
        t_index <- exog_start + i
        for (k in 1:K) {
          eps[t_index, k] <- eps[t_index, k] + 
            exog_shock_amplitude[k] * (exog_shock_decay^i) * sign_vec
        }
      }
    }
  }
  
  return(eps)
}


hhelper_make_stationary_A <- function(A_list, K) {
  #' Helper helper function which adjusts the coefficient matrix A such that it 
  #' is stationary if needed.
  #' 
  #' Parameters:
  #' - A_list (list): List of length p with K x K coefficient matrices from 
  #'                  helper_coefficient_generator.
  #' 
  #' Returns: 
  #' - Adjusted A_list, if it was unstationary.

  p <- length(A_list)
  A_sum <- Reduce(`+`, A_list)
  
  repeat {
    # also check spectral radius
    eigvals <- eigen(A_sum)$values
    rho     <- max(Mod(eigvals))
    
    # scale down if necessary
    if (rho >= 0.99) {
      # scale factor = 1.01 * rho or something slightly bigger than rho
      scale_factor <- 1.01 * rho
      if (scale_factor < 1.01) scale_factor <- 1.01
      
      for (i in seq_len(p)) {
        A_list[[i]] <- A_list[[i]] / scale_factor
      }
      A_sum <- A_sum / scale_factor
    } else {
      break
    }
  }

  return(A_list)
}




K = 10
p = 3
T = 40

# no shocks
ts.plot(sim_y(K, p, T, 
              min_indiv_shocks = 0, 
              max_indiv_shocks = 0, 
              min_high_vol_periods = 0, 
              max_high_vol_periods = 0,
              min_exog_shocks = 0,
              max_exog_shocks = 0), 
        col=1:K, main= "No shocks", ylim = c(-0.25, 0.25))
abline(h = 0, col = "black", lty = 2)

# individual shocks
ts.plot(sim_y(K, p, T, 
              min_indiv_shocks = 0, 
              max_indiv_shocks = 2, 
              min_high_vol_periods = 0, 
              max_high_vol_periods = 0,
              min_exog_shocks = 0,
              max_exog_shocks = 0), 
        col=1:K, main= "Individual shocks", ylim = c(-0.25, 0.25))
abline(h = 0, col = "black", lty = 2)

# exogenous shocks
ts.plot(sim_y(K, p, T, 
              min_indiv_shocks = 0, 
              max_indiv_shocks = 0, 
              min_high_vol_periods = 1, 
              max_high_vol_periods = 2,
              min_exog_shocks = 1,
              max_exog_shocks = 2), 
        col=1:K, main= "Exogenous shocks", ylim = c(-0.25, 0.25))
abline(h = 0, col = "black", lty = 2)

# individual and exogenous shocks
ts.plot(sim_y(K, p, T,
              min_indiv_shocks = 0, 
              max_indiv_shocks = 2, 
              min_high_vol_periods = 1, 
              max_high_vol_periods = 2,
              min_exog_shocks = 1,
              max_exog_shocks = 2), 
        col=1:K, main= "Individual and exogenous shocks", ylim = c(-0.25, 0.25))
abline(h = 0, col = "black", lty = 2)
