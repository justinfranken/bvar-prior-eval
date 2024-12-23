# simulating financial growth rates for Monte Carlo Simulations


coefficient_generator <- function(p, 
                                  K, 
                                  d_min = 0.2, 
                                  d_max = 0.4, 
                                  off_d_mean = 0.02, 
                                  off_d_sd=0.01){
  #' Creates a list of p many K x K VAR coefficient matrices, which decay with increasing lags p.
  #' 
  #' Parameters:
  #' - p (scalar): How many lags are included.
  #' - K (scalar): How many variables are included.
  #' - d_min (scalar): Minimum of the diagonal elements in the coefficient matrix.
  #' - d_max (scalar): Maximum of the diagonal elements in the coefficient matrix.
  #' - off_d_mean (scalar): Mean of off-diagonal elements in the coefficient matrix.
  #' - off_d_sd (scalar): Standard deviation of off-diagonal elements in the coefficient matrix.
  #' 
  #' Returns:
  #' - A list of VAR coefficient matrices.

  A_list <- vector("list", p)
  for (i in 1:p) {
    
    diag_min <- d_min / i
    diag_max <- d_max / i
    off_diag_mean <- off_d_mean / i
    off_diag_sd   <- off_d_sd / i
    
    A_raw <- matrix(0, nrow=K, ncol=K)
    
    # persisting diagonal elements:
    diag(A_raw) <- runif(K, min=diag_min, max=diag_max)
    
    # off-diagonal elements:
    for (r in 1:K) {
      for (col in 1:K) {
        if (r != col) {
          A_raw[r, col] <- rnorm(1, mean=off_diag_mean, sd=off_diag_sd) 
        }
      }
    }
    # check eigenvalues for stability: If max eigenvalue >= 1, scale down
    eigenvals <- eigen(A_raw)$values
    if (max(Mod(eigenvals)) >= 1) {
      scale_factor <- 1.1 * max(Mod(eigenvals))
      A_raw <- A_raw / scale_factor
    }
    
    A_list[[i]] <- A_raw
  }
  return(A_list)
}


eps_generator <- function(K, 
                          T, 
                          shock_diag_min = 0.5, 
                          shock_diag_max = 0.85, 
                          mean_vola = 0.2, 
                          sd_vola = 0.005){
  #' Simulates the shocks of the VAR model, epsilon.
  #' 
  #' Parameters:
  #' - K (scalar): How many variables are included.
  #' - T (scalar): How many observations are going to be simulated.
  #' - shock_diag_min (scalar): Minimum of the off diagonal element in Rho.
  #' - shock_diag_max (scalar): Maximum of the off diagonal element in Rho.
  #' - mean_vola (scalar): Mean of the firms volatility.
  #' - sd_vola (scalar): Standard deviation of the firms volatility.
  #' 
  #' Returns:
  #' - A T*K matrix containing our simulation epsilons.
  
  Rho <- matrix(0, nrow = K, ncol = K)
  
  diag(Rho) <- 1
  
  # Off-diagonal: random draws only for upper triangle, then mirror
  for (r in 1:(K-1)) {
    for (col in (r+1):K) {
      x <- runif(1, min = shock_diag_min, max = shock_diag_max)
      Rho[r, col] <- x
      Rho[col, r] <- x
    }
  }
  
  # Set standard deviations: assume the leading firm has slightly lower volatility
  sd_vec <-  rnorm(K, mean = mean_vola, sd = sd_vola) 
  Sigma <- diag(sd_vec) %*% Rho %*% diag(sd_vec)
  
  # We must generate normal random vectors with covariance Sigma.
  # Without external libraries, we can use a Cholesky decomposition:
  cholSigma <- chol(Sigma)
  Z <- matrix(rnorm(T*K), T, K)
  eps <- Z %*% cholSigma  # eps now has covariance Sigma
  
  return(eps)
}


sim_y <- function(K, 
                  p, 
                  T,
                  avg_growth_mean = 0.015,
                  avg_growth_sd = 0.01,
                  d_min_coef = 0.2, 
                  d_max_coef = 0.4, 
                  off_d_mean_coef = 0.02, 
                  off_d_sd_coef=0.01,
                  diag_min_shock = 0.001, 
                  diag_max_shock = 0.4, 
                  mean_vola_shock = 0.02, 
                  sd_vola_shock = 0.005){
  #' Creates T simulations of a K*1 vector y of a VAR(p) process.
  #' 
  #' Parameters:
  #' - K (scalar): How many variables are included.
  #' - p (scalar): How many lags are included.
  #' - T (scalar): How many observations are going to be simulated.
  #' - avg_growth_mean (scalar): Average growth rate of y.
  #' - avg_growth_sd (scalar): Standard deviation of average growth rate of y.
  #' - d_min_coef (scalar): Minimum of the diagonal elements in the coefficient matrix.
  #' - d_max_coef (scalar): Maximum of the diagonal elements in the coefficient matrix.
  #' - off_d_mean_coef (scalar): Mean of off-diagonal elements in the coefficient matrix.
  #' - off_d_sd_coef (scalar): Standard deviation of off-diagonal elements in the coefficient matrix.
  #' - diag_min_shock (scalar): Minimum of the off diagonal element in Rho.
  #' - diag_max_shock (scalar): Maximum of the off diagonal element in Rho.
  #' - mean_vola_shock (scalar): Mean of the firms volatility.
  #' - sd_vola_shock (scalar): Standard deviation of the firms volatility.
  #' 
  #' Returns:
  #' - K*1 matrix of simulated y's.
  
  # y average growth rate
  c <- rnorm(K, mean = avg_growth_mean, sd = avg_growth_sd)
  
  # compute coefficient matrix A
  A_list <- coefficient_generator(p, K, d_min_coef, d_max_coef, off_d_mean_coef, off_d_sd_coef)
  
  # compute shocks eps
  eps <- eps_generator(K, T, diag_min_shock, diag_max_shock, mean_vola_shock, sd_vola_shock)
  
  # initial p lags
  Y <- matrix(0, nrow = T, ncol = K)
  for (i in 1:p) {
    Y[i,] <- c + rnorm(K, mean = 0, sd = sd_vola_shock)
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

K = 5
p = 5
T = 50
ts.plot(sim_y(K, p, T), col=1:K, main="Simulated Growth Rates of Comparable Firms")
legend("topright", legend=paste("Firm", 1:K), col=1:K, lty=1)
