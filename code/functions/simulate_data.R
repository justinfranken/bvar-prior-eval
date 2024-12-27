# simulating financial growth rates for Monte Carlo Simulations


sim_y <- function(K,
                  p, 
                  T,
                  target_mu = 0.015,
                  sd_of_mu = 0.01,
                  # --- parameters for coefficient matrix A ---------------------
                  d_min_coef = 0.2, 
                  d_max_coef = 0.4, 
                  off_d_mean_coef = 0.05, 
                  off_d_sd_coef = 0.01,
                  # --- parameters for time series shocks ----------------------
                  shock_diag_min = 0.3, 
                  shock_diag_max = 0.6, 
                  mean_vola = 0.015, 
                  sd_vola = 0.005,
                  min_indiv_shocks = 0,
                  max_indiv_shocks = 2,
                  shock_length_min = 1,
                  shock_length_max = 5,
                  shock_ampl_min   = 0.05,
                  shock_ampl_max   = 0.15,
                  shock_decay      = 0.5){
  #' Creates T simulations of a K*1 vector y of a VAR(p) process.
  #' 
  #' Parameters:
  #' - K (integer): How many variables are included.
  #' - p (integer): How many lags are included.
  #' - T (integer): How many observations are going to be simulated.
  #' - target_mu, sd_of_mu (scalar): Average and standard deviation of growth rate from y.
  #' - d_min_coef, d_max_coef (scalar): Range of the diagonal elements in the coefficient matrix.
  #' - off_d_mean_coef, off_d_sd_coef (scalar): Mean and standard deviation of off-diagonal
  #'                                            elements in the coefficient matrix.
  #' - shock_diag_min, shock_diag_max (scalar): Range for the off-diagonal correlations.
  #' - mean_vola, sd_vola (scalar): Mean and SD for drawing each firm's volatility.
  #' - min_indiv_shocks, max_indiv_shocks (integer): Range of how many big shocks can appear
  #'                                                 in each variable's time series.
  #' - shock_length_min, shock_length_max (scalar): Range of durations for the big shocks.
  #' - shock_ampl_min, shock_ampl_max (scalar): Range for the initial amplitude of each shock.
  #' - shock_decay (scalar): Decay factor for the shock. 
  #' 
  #' Returns:
  #' - K*1 matrix of simulated y's.

  # compute coefficient matrix A
  A_list <- helper_coefficient_generator(p, 
                                         K, 
                                         d_min_coef, 
                                         d_max_coef, 
                                         off_d_mean_coef, 
                                         off_d_sd_coef)
  
  # compute average growth rates
  c <- helper_average_growth(target_mu, 
                             sd_of_mu, 
                             A_list, 
                             K)
  
  # compute shocks
  eps <- helper_eps_generator(K, 
                              T, 
                              shock_diag_min, 
                              shock_diag_max, 
                              mean_vola, 
                              sd_vola,
                              min_indiv_shocks,
                              max_indiv_shocks,
                              shock_length_min,
                              shock_length_max,
                              shock_ampl_min,
                              shock_ampl_max,
                              shock_decay)
  
  # initial p lags
  Y <- matrix(0, nrow = T, ncol = K)
  for (i in 1:p) {
    Y[i,] <- c + rnorm(K, mean = 0, sd = sd_vola)
  }
  
  # generate the series
  for (t in (p+1):T) {
    # compute the VAR component
    Y_t <- c
    for (lag in 1:p) {
      Y_t <- Y_t + A_list[[lag]] %*% Y[t-lag,]
    }
    # add the shock
    Y[t,] <- Y_t + eps[t,]
  }
  
  return(Y)
}


helper_coefficient_generator <- function(p, 
                                  K, 
                                  d_min = 0.2, 
                                  d_max = 0.4, 
                                  off_d_mean = 0.02, 
                                  off_d_sd=0.01){
  #' Helper function which creates a list of p many K*K VAR coefficient 
  #' matrices, which decay lineally with increasing lags p.
  #' 
  #' Parameters:
  #' - p (integer): How many lags are included.
  #' - K (integer): How many variables are included.
  #' - d_min_coef, d_max_coef (scalar): Range of the diagonal elements in the coefficient matrix.
  #' - off_d_mean_coef, off_d_sd_coef (scalar): Mean and standard deviation of off-diagonal
  #'                                            elements in the coefficient matrix.
  #' 
  #' Returns:
  #' - A list of p-many K*K VAR coefficient matrices.

  A_list <- vector("list", p)
  for (i in 1:p) {
    
    diag_min <- d_min / i
    diag_max <- d_max / i
    off_diag_mean <- off_d_mean / i
    off_diag_sd   <- off_d_sd / i
    
    A_raw <- matrix(0, nrow=K, ncol=K)
    
    # persisting diagonal elements
    diag(A_raw) <- runif(K, min=diag_min, max=diag_max)
    
    # off-diagonal elements
    for (r in 1:K) {
      for (col in 1:K) {
        if (r != col) {
          A_raw[r, col] <- rnorm(1, mean=off_diag_mean, sd=off_diag_sd) 
        }
      }
    }
    
    A_list[[i]] <- A_raw
  }
  
  # check if stable
  A_list <- hhelper_make_stationary_A(A_list, K)
  
  return(A_list)
}


helper_average_growth <- function(target_mu, sd_of_mu, A_list, K){
  #' Helper function to create average growth vector for the VAR(p) process.
  #' 
  #' Parameters:
  #' - target_mu (scalar): Average growth rate of y.
  #' - sd_of_mu (scalar): Standard deviation of average growth rates of each 
  #'                      variable K in y.
  #' - A_list (list): List output from helper_coefficient_generator.
  #' - K (integer): Number of variables in y.
  #' 
  #' Returns:
  #' - Vector of average growth for the VAR(p) process.
  
  A_sum <- Reduce(`+`, A_list)
  
  mu <- rnorm(K, mean = target_mu, sd = sd_of_mu)
  
  tmp <- diag(K) - A_sum
  c_vec <- tmp %*% mu
  
  return(c_vec)
}


helper_eps_generator <- function(K, 
                                 T, 
                                 shock_diag_min = 0.5, 
                                 shock_diag_max = 0.85, 
                                 mean_vola = 0.2, 
                                 sd_vola = 0.005,
                                 min_indiv_shocks = 0,
                                 max_indiv_shocks = 3,
                                 shock_length_min = 1,
                                 shock_length_max = 5,
                                 shock_ampl_min   = 0.05,
                                 shock_ampl_max   = 0.15,
                                 shock_decay      = 0.5) {
  #' Helper function to creates a T*K matrix of shocks. First, draws from a 
  #' multivariate normal distribution with a correlation structure Rho. Then adds 
  #' short-lived "individual shocks" for each variable K, which revert quickly.
  #'
  #' Parameters:
  #' - K (integer): Number of variables.
  #' - T (integer): Number of observations.
  #' - shock_diag_min, shock_diag_max (scalar): Range for the off-diagonal correlations.
  #' - mean_vola, sd_vola (scalar): Mean and SD for drawing each firm's volatility.
  #' - min_indiv_shocks, max_indiv_shocks (integer): Range of how many "big shocks" 
  #'        can appear in each variable's time series.
  #' - shock_length_min, shock_length_max (scalar): Range of durations for these big shocks.
  #' - shock_ampl_min, shock_ampl_max (scalar): Range for the initial amplitude of each shock.
  #' - shock_decay (scalar): Decay factor for the shock. 
  #'
  #' Returns:
  #' - A T*K matrix containing the simulation shocks.
  
  # create baseline correlation matrix
  Rho <- matrix(0, nrow = K, ncol = K)
  diag(Rho) <- 1
  
  # off-diagonal: random draws (for upper triangle), then mirror
  for (r in 1:(K-1)) {
    for (col in (r+1):K) {
      x <- runif(1, min = shock_diag_min, max = shock_diag_max)
      Rho[r, col] <- x
      Rho[col, r] <- x
    }
  }
  
  # adjust Rho if it is not positive (semi) definite
  Rho <- as.matrix(nearPD(Rho, corr = TRUE)$mat)
  
  # standard deviations: each firm has random volatility
  sd_vec <- rnorm(K, mean = mean_vola, sd = sd_vola)
  Sigma  <- diag(sd_vec) %*% Rho %*% diag(sd_vec)
  
  # generate the base normal random vectors with covariance Sigma
  cholSigma <- chol(Sigma)
  Z         <- matrix(rnorm(T * K), nrow = T, ncol = K)
  eps       <- Z %*% cholSigma
  
  # add individual larger shocks 
  num_shocks_vec <- sample(min_indiv_shocks:max_indiv_shocks, K, replace = TRUE)
  
  for (k in seq_len(K)) {
    n_shocks <- num_shocks_vec[k]
    
    if (n_shocks > 0) {
      for (s in 1:n_shocks) {
        shock_length  <- sample(shock_length_min:shock_length_max, 1)
        shock_start   <- sample(seq_len(T - shock_length), 1)
        
        shock_amplitude <- runif(1, shock_ampl_min, shock_ampl_max)
        shock_sign      <- sample(c(-1, 1), 1)  
        shock_amplitude <- shock_amplitude * shock_sign
        
        # reverting back to the mean
        for (i in 0:(shock_length - 1)) {
          t_index <- shock_start + i
          eps[t_index, k] <- eps[t_index, k] + shock_amplitude * (shock_decay^i)
        }
      }
    }
  }
  
  return(eps)
}


hhelper_make_stationary_A <- function(A_list, K, tol = 1e-10) {
  #' Helper helper function which adjusts the coefficient matrix A such that it 
  #' is stationary if needed.
  #' 
  #' Parameters:
  #' - A_list (list): List of length p with K x K coefficient matrices from 
  #'                  helper_coefficient_generator.
  #' - tol (scalar): small tolerance for numerical checks.
  #' 
  #' Returns: 
  #' - Adjusted A_list, if it was unstationary.

  p <- length(A_list)
  A_sum <- Reduce(`+`, A_list)
  
  repeat {
    # check if (I - A_sum) is invertible
    tmp <- diag(K) - A_sum
    d <- det(tmp)

    # also check spectral radius
    eigvals <- eigen(A_sum)$values
    rho     <- max(Mod(eigvals))
    
    # scale down if necessary
    if ((abs(d) < tol) || (rho >= 1)) {
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




K = 5
p = 5
T = 50
ts.plot(sim_y(K, p, T), col=1:K, main="Simulated Growth Rates of Comparable Firms")
mean(colMeans(sim_y(K, p, T)))
colMeans(sim_y(K, p, T))
