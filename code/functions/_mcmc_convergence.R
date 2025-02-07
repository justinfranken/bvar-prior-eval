# functions to assess convergence of mcmc chains


geweke_coef <- function(out_obj, coef, intercept = TRUE, alpha = 0.05){
  #' Performs a Geweke convergence diagnostic on MCMC draws of BVAR coefficients, according
  #' to Geweke (1992).
  #'
  #' Parameters:
  #' - out_obj (list): An object containing MCMC draws.
  #' - coef (integer): Either 1 for Phi or 2 for Sigma. 
  #' - intercept (logical, default = TRUE): Indicates whether an intercept is included in the BVAR model.
  #' - alpha (numeric, default = 0.05): The significance level used to evaluate the z-scores.
  #'
  #' Returns:
  #' - A list containing:
  #'    - logicals: A matrix (m x k) of logical values, where TRUE indicates that the absolute Geweke
  #'                z-score for that coefficient is less than the critical value from the standard normal
  #'                distribution (i.e., convergence is acceptable at the given alpha level).
  #'    - z_values: A matrix (m x k) of the rounded Geweke z-scores for each coefficient.
  
  k <- dim(out_obj[[coef]][,,1])[2]
  m <- dim(out_obj[[coef]][,,1])[1]
  p <- (m - as.integer(intercept)) / k
  
  # matrix names
  varNames <- paste0("Y", 1:k)
  rowNames <- c("Intercept")
  for (i in 1:p) {
    rowNames <- c(rowNames, paste0("Lag", rep(i, k)))
  }
  if(!intercept) rowNames <- rowNames[-1]
  
  # storage
  out_logical <- matrix(NA, nrow = m, ncol = k,
                        dimnames = list(rowNames, varNames))
  out_values <- matrix(NA, nrow = m, ncol = k,
                       dimnames = list(rowNames, varNames))
  
  # get geweke scores
  for (row in seq_len(m)) {
    for (col in seq_len(k)) {
      temp <- geweke.diag(as.mcmc(out_obj[[coef]][row,col,]))$z
      out_logical[row,col] = abs(temp) < qnorm(1 - alpha / 2)
      out_values[row,col] = round(temp, digits = 3)
    }
  }
  
  # count how many TRUEs exist
  prct_TRUE <- sum(out_logical) / (((p * k) + as.integer(intercept)) * k)
  
  return(list(
    logicals = out_logical,
    z_values = out_values,
    prct_TRUE = prct_TRUE
  ))
}


geweke_hyp <- function(out_obj, alpha = 0.05){
  #' Performs the Geweke convergence diagnostic on hyperparameter draws, according
  #' to Geweke (1992).
  #'
  #' Parameters:
  #' - out_obj (list): An object containing MCMC output hyperparameters.
  #' - alpha (numeric, default = 0.05): The significance level for the convergence test.
  #'
  #' Returns:
  #' - A list with the following components:
  #'   - logicals: A logical vector indicating for each hyperparameter whether the absolute Geweke z-score
  #'               is less than the critical value.
  #'   - z_values: A numeric vector of the Geweke z-scores for each hyperparameter.
  #'   - prct_TRUE: The proportion of hyperparameters that passed the convergence diagnostic.
  
  len_hyp <- dim(out_obj$hyper_parameters)[2]
  
  # storage
  out_logical <- rep(NA, len_hyp)
  out_values <- rep(NA, len_hyp)
  
  # get geweke scores
  for (col in seq_len(len_hyp)) {
    temp <- geweke.diag(out_obj$hyper_parameters[,col])$z
    out_logical[col] <- abs(temp) < qnorm(1 - alpha / 2)
    out_values[col] <- temp
  }
  
  # count how many TRUEs exist
  prct_TRUE <- sum(out_logical) / len_hyp
  
  return(list(
    logicals = out_logical,
    z_values = out_values,
    prct_TRUE = prct_TRUE
  ))
}


traceplot_coef <- function(out_obj, coef, coef_row, coef_col){
  #' Plots the trace of the BVAR coefficients over the MCMC chain.
  #'
  #' Parameters:
  #' - out_obj (list): A list containing the posterior draws of the BVAR coefficients.
  #' - coef (integer): Either 1 for Phi or 2 for Sigma. 
  #' - coef_row (integer): The row index in the Phi matrix corresponding to the coefficient to plot.
  #' - coef_col (integer): The column index in the Phi matrix corresponding to the coefficient to plot.
  #'
  #' Returns:
  #' - A traceplot of Phi is produced. The plot shows:
  #'    - The chain draws on the x-axis.
  #'    - The coefficient values on the y-axis.
  #'    - A red dashed line at zero.
  #'    - A green dashed line at the median value of the coefficient draws.

  plot(
    out_obj[[coef]][coef_row,coef_col,], 
    type = "l", 
    main = paste0("Traceplot of ", ifelse(coef == 1, "Phi", "Sigma"),"[", coef_row, ",", coef_col, "]",
    " with median ", 
    round(median(out_obj[[coef]][coef_row,coef_col,]), digits = 2),
    " (green)"),
  xlab = "Itteration j",
  ylab = "Coefficient value")
  abline(h = 0, col = "red", lty = 2)
  abline(h = median(out_obj[[coef]][coef_row,coef_col,]), col = "green", lty = 2)
}


traceplot_hyper <- function(out_obj, hyp_para){
  #' Plots the trace of a hyperparameter over the MCMC chain.
  #'
  #' Parameters:
  #' - out_obj (list): A list containing MCMC output of the run_bvar_hierarch() function.
  #' - hyp_para (integer): The column index for the hyperparameter to be plotted.
  #'
  #' Returns:
  #' - A traceplot is produced. The plot shows:
  #'    - Chain draws on the x-axis.
  #'    - Hyperparameter values on the y-axis.
  #'    - A red dashed line at 0.
  #'    - A green dashed line at the median of the hyperparameter draws.
  
  plot(
    out_obj$hyper_parameters[,hyp_para], 
    type = "l",
    main = paste0("Traceplot of hyperparameter ", 
                  hyp_para, " with median ", 
                  round(median(out_obj$hyper_parameters[,hyp_para]), digits = 2),
                  " (green) and acceptance rate ",
                  round(
                    length(
                      unique(out_obj$hyper_parameters[,hyp_para]
                      )
                    )/length(
                      out_obj$hyper_parameters[,hyp_para]), 
                    digits = 2)
    ),
    xlab = "Itteration j",
    ylab = "Hyperparameter value"
  )
  abline(h = 0, col = "red", lty = 2)
  abline(h = median(out_obj$hyper_parameters[,hyp_para]), col = "green", lty = 2)
}

traceplot_hyper <- function(out_obj, hyp_para){
  #' Plots the trace of a hyperparameter over the MCMC chain.
  #'
  #' Parameters:
  #' - out_obj (list): A list containing MCMC output of the run_ssvs() function.
  #' - hyp_para (integer): The column index in for the hyperparameter to be plotted.
  #'
  #' Returns:
  #' - A traceplot is produced. The plot shows:
  #'    - Chain draws on the x-axis.
  #'    - Hyperparameter values on the y-axis.
  #'    - A red dashed line at 0.
  #'    - A green dashed line at the median of the hyperparameter draws.
  plot(out_obj$delta_draws[,hyp_para], 
       type = "l",
       main = paste0("Traceplot of hyperparameter ", hyp_para, 
                     " with mean ", mean(out_obj$delta_draws[,hyp_para]), " (green)"),
       xlab = "Itteration j",
       ylab = "Hyperparameter value",
       ylim = c(0, 1)
  )
  abline(h = mean(out_obj$delta_draws[,hyp_para]), col = "green", lty = 2)
  abline(h = 0, col = "red", lty = 2)
}


brooks_plot <- function(out_obj, B, coef, coef_row, coef_col){
  #' Generates a Brooks-style CUSUM hairiness plot for a selected BVAR coefficient.
  #'
  #' Parameters:
  #' - out_obj (list): An MCMC output object that contains the posterior draws of BVAR coefficients.
  #' - B (integer): The burn-in index, indicating the number of initial draws to discard.
  #' - coef (integer): Either 1 for Phi or 2 for Sigma. 
  #' - coef_row (integer): The row index of the coefficient to be plotted.
  #' - coef_col (integer): The column index of the coefficient to be plotted.
  #'
  #' Returns:
  #' - A Brooks-style CUSUM hairiness plot for the selected coefficient is produced.
  
  chain <- out_obj[[coef]][coef_row,coef_col,]
  
  # get smoothness running means
  res <- helper_cusum_hairiness(chain, B)
  D_vals <- res$D
  lower <- res$L
  upper <- res$U
  
  # plot only after burn-in
  t_seq <- (B+2):length(D_vals)
  
  # create plot
  plot(
    x    = t_seq,
    y    = D_vals[t_seq],
    type = "l",
    xlab = "Iteration t",
    ylab = expression(D[t]),
    main = paste0("CUSUM Hairiness Measure of ", ifelse(coef == 1, "Phi", "Sigma"), "[", coef_row, ",", coef_col, "]")
  )
  lines(t_seq, lower[t_seq], col = "red", lty = 2)
  lines(t_seq, upper[t_seq], col = "red", lty = 2)
  abline(h = 0.5, col = "gray")
}


helper_cusum_hairiness <- function(chain, B, alpha = 0.05) {
  #' Computes CUSUM-based hairiness indicators and confidence bounds for an MCMC chain,
  #' according to Brooks (1998).
  #'
  #' Parameters:
  #' - chain (numeric vector): A vector of MCMC draws of length R.
  #' - B (integer): Burn-in index (number of initial draws to discard).
  #' - alpha (numeric, default = 0.05): Significance level for the confidence bounds.
  #'
  #' Returns:
  #' - A list containing:
  #'   - D: A numeric vector of running means (D_t) computed for t = B+2, ..., R.
  #'   - L: A numeric vector of lower confidence bounds for D_t.
  #'   - U: A numeric vector of upper confidence bounds for D_t.
  
  R <- length(chain)
  
  # posterior mean estimate after burn-in
  muhat <- mean(chain[(B+1):R])
  
  # compute the CUSUM
  S <- numeric(R + 1)
  S[1:B] <- NA
  S[B] <- 0
  
  for (t in (B+1):R) {
    S[t] <- S[t-1] + (chain[t] - muhat)
  }
  
  # hairiness indicator d_j
  
  d <- numeric(R + 1)
  for (j in (B+2):(R-1)) {
    if ((S[j-1] > S[j] && S[j] < S[j+1]) ||
        (S[j-1] < S[j] && S[j] > S[j+1])) {
      d[j] <- 1
    } else {
      d[j] <- 0
    }
  }
  
  # running means of d_j
  D <- rep(NA, R + 1)
  for (t in (B+2):R) {
    D[t] <- mean(d[(B+1):(t-1)])
  }
  
  # confidence bounds
  zval <- qnorm(1 - alpha / 2)
  L <- rep(NA, R + 1)
  U <- rep(NA, R + 1)
  for (t in (B+2):R) {
    width <- zval * sqrt(1 / (4 * (t - B - 1)))
    L[t] <- -width + 0.5
    U[t] <-  width + 0.5
  }
  
  return(list(
    D     = D,
    L     = L,
    U     = U
  ))
}
