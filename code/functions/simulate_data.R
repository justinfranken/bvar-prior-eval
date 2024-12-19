# simulating financial growth rates for Monte Carlo Simulations


coefficient_generator <- function(p, K, d_min = 0.2, d_max = 0.4, off_d_mean = 0.02, off_d_sd=0.01){
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
  #' - a list of VAR coefficient matrices.

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

