# functions needed for ssvs bvar analysis


sherman_morrison_update <- function(A_inv, A_logdet, u, v) {
  #' Computes a rank-1 update of the inverse matrix and its log-determinant using the Sherman-Morrison formula.
  #'
  #' Parameters:
  #' - A_inv (matrix): The inverse of matrix A.
  #' - A_logdet (numeric): The logarithm of the determinant of matrix A.
  #' - u (vector or matrix): The vector u in the rank-1 update.
  #' - v (vector or matrix): The vector v in the rank-1 update.
  #'
  #' Returns:
  #' - A list containing:
  #'   - A_inv: The updated inverse of matrix A_new = A + u v^T.
  #'   - logdet: The updated log-determinant.

  denom <- as.numeric(1 + (t(v) %*% A_inv %*% u))
  # if (abs(denom) < 1e-15) {
  #   stop("Rank-1 update denominator numerically zero; check stability.")
  # }
  
  A_inv_u <- A_inv %*% u
  vT_A_inv <- t(v) %*% A_inv
  
  # new inverse
  A_new_inv <- A_inv - (A_inv_u %*% vT_A_inv) / denom
  
  # log-det update
  A_new_logdet <- A_logdet + log(abs(denom))
  
  return(list(A_inv = A_new_inv, logdet = A_new_logdet))
}


build_Omega_lowerbar <- function(delta, tau0, tau1, se_ols) {
  #' Constructs a diagonal matrix Omega_lowerbar based on delta and scaling parameters.
  #'
  #' Parameters:
  #' - delta (numeric vector): A vector indicating which scaling factor to use for each diagonal element.
  #' - tau0 (numeric): The scaling factor used when the corresponding delta is zero.
  #' - tau1 (numeric): The scaling factor used when the corresponding delta is non-zero.
  #' - se_ols (numeric vector): The standard error from an OLS regression, used as a baseline scaling factor.
  #'
  #' Returns:
  #' - A diagonal matrix with diagonal entries computed as:
  #'   - (tau0 * se_ols)^2 if the corresponding delta equals zero;
  #'   - (tau1 * se_ols)^2 otherwise.
  
  diag_entries <- ifelse(delta == 0,
                         (tau0 * se_ols)^2,
                         (tau1 * se_ols)^2)
  return(diag(diag_entries))
}


build_inv_Omega_lowerbar <- function(delta, tau0, tau1, se_ols) {
  #' Constructs the inverse of the diagonal matrix Omega_lowerbar.
  #'
  #' Parameters:
  #' - delta (numeric vector): A vector indicating which scaling factor to use for each diagonal element.
  #' - tau0 (numeric): The scaling factor used when the corresponding delta is zero.
  #' - tau1 (numeric): The scaling factor used when the corresponding delta is non-zero.
  #' - se_ols (numeric vector): The standard error from an OLS regression, used as a baseline scaling factor.
  #'
  #' Returns:
  #' - A diagonal matrix representing the inverse of Omega_lowerbar, with diagonal entries:
  #'   - 1/((tau0 * se_ols)^2) if the corresponding delta equals zero;
  #'   - 1/((tau1 * se_ols)^2) otherwise.
  
  diag_entries <- ifelse(delta == 0,
                         1/((tau0 * se_ols)^2),
                         1/((tau1 * se_ols)^2))
  return(diag(diag_entries))
}


pi_delta <- function(delta, p = 0.8) {
  #' Computes the log-prior probability of delta under a Bernoulli model in log-scale.
  #'
  #' Parameters:
  #' - delta (numeric vector): A binary vector (elements are 0 or 1) representing the state indicators.
  #' - p (numeric, default = 0.8): The success probability for the Bernoulli distribution.
  #'
  #' Returns:
  #' - A numeric value representing the log-prior probability of delta, computed as:
  
  m <- length(delta)
  s <- sum(delta)
  return(s*log(p) + (m - s)*log(1 - p))
}


compute_sigma_ols_VAR <- function(Y, X, XX_inv) {
  #' This function computes the maximum standard errors for each equation in a Vector Autoregression (VAR)
  #' model estimated by Ordinary Least Squares (OLS).
  #' 
  #' Parameters:
  #' - Y (matrix): A T x k matrix of dependent variables.
  #' - X (matrix): A T x m matrix of regressors.
  #' - XX_inv (matrix): The inverse of X'X, used in OLS estimation.
  #' 
  #' Returns:
  #' - A vector of length m, containing the maximum standard errors for each equation.
  
  T <- nrow(Y)
  k <- ncol(Y)
  m <- ncol(X)
  
  all_se <- vector("list", k)
  
  for (i in seq_len(k)) {
    y_i <- Y[, i]
    
    b_i <- XX_inv %*% (t(X) %*% y_i)
    
    resid_i <- y_i - X %*% b_i
    sigma2_i <- c(crossprod(resid_i)) / (T - m)
    
    var_b_i <- sigma2_i * XX_inv
    
    se_b_i <- sqrt(diag(var_b_i))
    
    all_se[[i]] <- se_b_i
  }
  
  temp <- matrix(unlist(all_se), ncol = m, byrow = TRUE)
  out <- do.call("pmax", as.data.frame(t(temp)))
  
  return(out)
}


precompute_from_delta <- function(delta, 
                                  T, k,
                                  S_lowerbar, v_lowerbar,
                                  tau0, tau1, se_ols,
                                  XX, XX_inv, 
                                  diff_Phi, S,
                                  delta_prob) {
  #' Computes several initial matrices and log-determinants based on delta and other model parameters.
  #'
  #' Parameters:
  #' - delta (numeric vector): A binary vector (or a vector of indicators) used to select scaling factors.
  #' - T (integer): The effective number of observations.
  #' - k (integer): The number of variables.
  #' - S_lowerbar (matrix): The prior scale matrix for S.
  #' - v_lowerbar (numeric): The degrees of freedom associated with the prior scale matrix.
  #' - tau0 (numeric): Scaling factor applied when a delta element equals 0.
  #' - tau1 (numeric): Scaling factor applied when a delta element is 1.
  #' - se_ols (numeric vector): The standard error from an OLS regression used as a baseline scaling factor.
  #' - XX (matrix): The product X'X from the design matrix.
  #' - XX_inv (matrix): The inverse of XX.
  #' - diff_Phi (matrix or vector): The difference between the current Phi and its prior mean.
  #' - S (matrix): Squared error matrix t(epsilon) %*% epsilon.
  #' - delta_prob (numeric): The prior probability parameter for delta (used in the log-prior calculation).
  #'
  #' Returns:
  #' - A list containing:
  #'   - delta: The input delta.
  #'   - Omega_lowerbar: The diagonal matrix computed using delta, tau0, tau1, and se_ols.
  #'   - inv_Omega_lowerbar: The inverse of Omega_lowerbar.
  #'   - M: The matrix M = Omega_lowerbar + XX_inv.
  #'   - M_inv: The inverse of M.
  #'   - logdet_M: The log-determinant of M.
  #'   - P: The matrix P = inv_Omega_lowerbar + XX.
  #'   - P_inv: The inverse of P.
  #'   - logdet_P: The log-determinant of P.
  #'   - logdet_Omega_lowerbar: The log-determinant of Omega_lowerbar.
  #'   - S_upperbar: The posterior updated scale matrix.
  #'   - logdet_S_upperbar: The log-determinant of S_upperbar.
  #'   - val: The log posterior kernel value g(delta).
  
  Omega_lowerbar <- build_Omega_lowerbar(delta, tau0, tau1, se_ols)
  inv_Omega_lowerbar <- build_inv_Omega_lowerbar(delta, tau0, tau1, se_ols)
  
  # Omega
  logdet_Omega_lowerbar <- sum(log(diag(Omega_lowerbar)))
  
  # M
  M <- Omega_lowerbar + XX_inv
  M_chol <- chol(M)
  M_inv <- chol2inv(M_chol)
  logdet_M <- 2 * sum(log(diag(M_chol)))
  
  # P
  P <- inv_Omega_lowerbar + XX
  P_chol <- chol(P)
  P_inv <- chol2inv(P_chol)
  logdet_P <- 2 * sum(log(diag(P_chol)))
  
  # S
  S_upperbar <- S_lowerbar + S
  S_upperbar <- S_upperbar + t(diff_Phi) %*% M_inv %*% diff_Phi
  S_upperbar_chol <- chol(S_upperbar)
  logdet_S_upperbar <- 2 * sum(log(diag(S_upperbar_chol)))
  
  # log posterior kernel g(delta)
  val <- ( (logdet_Omega_lowerbar + logdet_P) * (-k/2) ) +
    ( logdet_S_upperbar * ( -(v_lowerbar + T)/2 ) ) +
    pi_delta(delta, delta_prob)
  
  # return list with everything
  return(list(delta = delta,
              Omega_lowerbar = Omega_lowerbar,
              inv_Omega_lowerbar = inv_Omega_lowerbar,
              M = M, M_inv = M_inv, logdet_M = logdet_M,
              P = P, P_inv = P_inv, logdet_P = logdet_P,
              logdet_Omega_lowerbar = logdet_Omega_lowerbar,
              S_upperbar = S_upperbar,
              logdet_S_upperbar = logdet_S_upperbar,
              val = val))
}


toggle_delta_rank1 <- function(precomp_list,
                               i,
                               T, k,
                               S_lowerbar, v_lowerbar,
                               tau0, tau1, se_ols,
                               XX, XX_inv,
                               diff_Phi, S,
                               delta_prob) {
  #' Toggles the i-th element of delta and updates precomputed matrices using a rank-1 update.
  #'
  #' Parameters:
  #' - precomp_list (list): A list of precomputed matrices and quantities for the current delta,
  #'   including delta, Omega_lowerbar, inv_Omega_lowerbar, M_inv, logdet_M, P_inv, logdet_P,
  #'   logdet_Omega_lowerbar, S_upperbar, and its log-determinant (val).
  #' - i (integer): The index of the delta element to toggle.
  #' - T (integer): The effective number of observations.
  #' - k (integer): The number of variables.
  #' - S_lowerbar (matrix): The prior scale matrix used in the computation of S_upperbar.
  #' - v_lowerbar (numeric): The prior degrees of freedom.
  #' - tau0 (numeric): The scaling factor used when delta equals 0.
  #' - tau1 (numeric): The scaling factor used when delta is 1.
  #' - se_ols (numeric vector): A vector of OLS standard errors, used for scaling.
  #' - XX (matrix): The matrix product X'X.
  #' - XX_inv (matrix): The inverse of XX.
  #' - diff_Phi (matrix or vector): The difference between the current Phi and its prior mean.
  #' - S (matrix): Squared error matrix t(epsilon) %*% epsilon.
  #' - delta_prob (numeric): The prior probability parameter for delta, used in the log-prior computation.
  #'
  #' Returns:
  #' - A list with the updated precomputed quantities after toggling the i-th element of delta:
  #'   - delta: The updated delta vector.
  #'   - Omega_lowerbar: The updated diagonal matrix Omega_lowerbar.
  #'   - inv_Omega_lowerbar: The updated inverse of Omega_lowerbar.
  #'   - M_inv: The updated inverse of matrix M.
  #'   - logdet_M: The updated log-determinant of M.
  #'   - P_inv: The updated inverse of matrix P.
  #'   - logdet_P: The updated log-determinant of P.
  #'   - logdet_Omega_lowerbar: The updated log-determinant of Omega_lowerbar.
  #'   - S_upperbar: The updated scale matrix S_upperbar.
  #'   - logdet_S_upperbar: The updated log-determinant of S_upperbar.
  #'   - val: The updated log posterior kernel value g(delta).

  delta_old <- precomp_list$delta
  old_val   <- precomp_list$val
  
  old_di <- delta_old[i]
  if (old_di == 0) {
    # flipping 0->1
    d   <- (tau1 * se_ols[i])^2 - (tau0 * se_ols[i])^2
    d_i <- (1/((tau1 * se_ols[i])^2)) - (1/((tau0 * se_ols[i])^2))
    delta_new_i <- 1
  } else {
    # flipping 1->0
    d   <- (tau0 * se_ols[i])^2 - (tau1 * se_ols[i])^2
    d_i <- (1/((tau0 * se_ols[i])^2)) - (1/((tau1 * se_ols[i])^2))
    delta_new_i <- 0
  }
  
  # new delta
  delta_new <- delta_old
  delta_new[i] <- delta_new_i
  

  e_i <- numeric(length(delta_old))
  e_i[i] <- 1
  
  # Sherman–Morrison on M_inv
  sm_M    <- sherman_morrison_update(precomp_list$M_inv,
                                     precomp_list$logdet_M,
                                     u = e_i * d, 
                                     v = e_i) 
  M_inv_new <- sm_M$A_inv
  logdet_M_new <- sm_M$logdet
  
  # Sherman–Morrison on P_inv
  sm_P    <- sherman_morrison_update(precomp_list$P_inv,
                                     precomp_list$logdet_P,
                                     u = e_i * d_i,
                                     v = e_i)
  P_inv_new <- sm_P$A_inv
  logdet_P_new <- sm_P$logdet
  
  # Omega lowerbar
  old_ii <- precomp_list$Omega_lowerbar[i,i]
  new_ii <- old_ii + d
  logdet_Omega_lowerbar_new <- precomp_list$logdet_Omega_lowerbar - log(old_ii) + log(new_ii)
  new_Omega_lowerbar <- precomp_list$Omega_lowerbar
  new_Omega_lowerbar[i,i] <- new_ii
  new_inv_Omega_lowerbar <- precomp_list$inv_Omega_lowerbar
  new_inv_Omega_lowerbar[i,i] <- 1/new_ii
  
  # build S_upperbar_new
  S_up_new <- S_lowerbar + S
  S_up_new <- S_up_new + t(diff_Phi) %*% M_inv_new %*% diff_Phi
  S_up_new_chol <- chol(S_up_new)
  logdet_S_up_new <- 2 * sum(log(diag(S_up_new_chol)))
  
  # compute g(delta_new)
  val_new <- ( (logdet_Omega_lowerbar_new + logdet_P_new) * (-k/2) ) +
    ( logdet_S_up_new * (-(v_lowerbar + T)/2) ) +
    pi_delta(delta_new, delta_prob)
  
  # return a new precomp list
  return(list(delta = delta_new,
              Omega_lowerbar = new_Omega_lowerbar,
              inv_Omega_lowerbar = new_inv_Omega_lowerbar,
              M_inv = M_inv_new, logdet_M = logdet_M_new, 
              P_inv = P_inv_new, logdet_P = logdet_P_new,
              logdet_Omega_lowerbar = logdet_Omega_lowerbar_new,
              S_upperbar = S_up_new,
              logdet_S_upperbar = logdet_S_up_new,
              val = val_new))
}


update_delta_rank1 <- function(precomp,
                               T, k,
                               S_lowerbar, v_lowerbar,
                               tau0, tau1, se_ols,
                               XX, XX_inv,
                               diff_Phi, S,
                               delta_prob) {
  #' Updates the delta vector via sequential rank-1 toggles using a Bernoulli acceptance step.
  #'
  #' Parameters:
  #' - precomp (list): A list of precomputed matrices and quantities for the current delta, including:
  #'     delta, Omega_lowerbar, inv_Omega_lowerbar, M_inv, logdet_M, P_inv, logdet_P,
  #'     logdet_Omega_lowerbar, S_upperbar, and val (the log-posterior kernel value).
  #' - T (integer): The effective number of observations.
  #' - k (integer): The number of variables.
  #' - S_lowerbar (matrix): The prior scale matrix for inverse Wishart.
  #' - v_lowerbar (numeric): Prior degrees of freedom associated with S_lowerbar.
  #' - tau0 (numeric): Scaling factor used when a delta element equals 0.
  #' - tau1 (numeric): Scaling factor used when a delta element is 1.
  #' - se_ols (numeric vector): A vector of OLS standard errors used for scaling.
  #' - XX (matrix): The product X'X from the design matrix.
  #' - XX_inv (matrix): The inverse of XX.
  #' - diff_Phi (matrix or vector): The difference between the current Phi and its prior mean.
  #' - S (matrix): Squared error matrix t(epsilon) %*% epsilon.
  #' - delta_prob (numeric): The prior probability parameter for delta (used in the log-prior computation).
  #'
  #' Returns:
  #' - precomp (list): The updated precomputed list after sequentially considering toggles of each element in delta.
  #'   The list includes the updated delta vector, matrices (Omega_lowerbar, inv_Omega_lowerbar, M_inv, P_inv),
  #'   their log-determinants, S_upperbar and its log-determinant, and the updated log-posterior kernel value.
  
  m <- length(precomp$delta)
  
  for (i in seq_len(m)) {
    # current log-posterior
    u_current <- precomp$val
    
    # toggle i => new precomputation
    toggled_precomp <- toggle_delta_rank1(precomp_list = precomp,
                                          i = i,
                                          T = T, k = k,
                                          S_lowerbar = S_lowerbar, v_lowerbar = v_lowerbar,
                                          tau0 = tau0, tau1 = tau1, se_ols = se_ols,
                                          XX = XX, XX_inv = XX_inv,
                                          diff_Phi = diff_Phi, S = S,
                                          delta_prob = delta_prob)
    u_flip <- toggled_precomp$val
    
    # Bernoulli probability
    lmax <- max(u_current, u_flip)
    p_i <- exp(u_flip - lmax) / (exp(u_flip - lmax) + exp(u_current - lmax))
    
    draw_i <- rbinom(n = 1, size = 1, prob = p_i)
    if (draw_i == 1) {
      # Accept the flip => updated precomp
      precomp <- toggled_precomp
    }
  
    # else do nothing => keep old precomp
  }
  
  return(precomp)
}


# ------------------------------------------------------------------------------
# old functions before implementing efficient rank-1 updates
# ------------------------------------------------------------------------------


#' build_Omega_lowerbar <- function(delta, tau0, tau1, se_ols) {
#'   #' This function constructs a diagonal matrix  based on a vector of binary indicators (delta). 
#'   #' Each diagonal element is determined by whether delta is 0 or 1, scaling the respective variances.
#'   #' 
#'   #' Parameters:
#'   #' - delta (vector): A vector of 0/1 indicators of length k, determining which scale to apply.
#'   #' - tau0, tau1 (scalars): Hyperparameters controlling the spike and slab adjustment given delta.
#'   #' - se_ols (vectors): A vector length k, providing standard deviation errors of the OLS residuals
#'   #'  for the coefficients based on delta.
#'   #'  
#'   #'  Returns:
#'   #'  - A k x k diagonal matrix where each diagonal entry is depending on delta.
#'   
#'   diag_entries <- ifelse(delta == 0, (tau0 * se_ols)^2, (tau1 * se_ols)^2)
#'   
#'   return(diag(diag_entries))
#' }
#' 
#' 
#' compute_sigma_ols_VAR <- function(Y, X, XX_inv) {
#'   #' This function computes the maximum standard errors for each equation in a Vector Autoregression (VAR)
#'   #' model estimated by Ordinary Least Squares (OLS).
#'   #' 
#'   #' Parameters:
#'   #' - Y (matrix): A T x k matrix of dependent variables.
#'   #' - X (matrix): A T x m matrix of regressors.
#'   #' - XX_inv (matrix): The inverse of X'X, used in OLS estimation.
#'   #' 
#'   #' Returns:
#'   #' - A vector of length m, containing the maximum standard errors for each equation.
#'   
#'   T <- nrow(Y)
#'   k <- ncol(Y)
#'   m <- ncol(X)
#'   
#'   all_se <- vector("list", k)
#'   
#'   for (i in seq_len(k)) {
#'     y_i <- Y[, i]
#'     
#'     b_i <- XX_inv %*% (t(X) %*% y_i)
#'     
#'     resid_i <- y_i - X %*% b_i
#'     sigma2_i <- c(crossprod(resid_i)) / (T - m)
#'     
#'     var_b_i <- sigma2_i * XX_inv
#'     
#'     se_b_i <- sqrt(diag(var_b_i))
#'     
#'     all_se[[i]] <- se_b_i
#'   }
#'   
#'   temp <- matrix(unlist(all_se), ncol = m, byrow = TRUE)
#'   out <- do.call("pmax", as.data.frame(t(temp)))
#'   
#'   return(out)
#' }
#' 
#' 
#' calc_g <- function(delta,
#'                    T, k,
#'                    S_lowerbar,
#'                    v_lowerbar,
#'                    tau0, tau1, se_ols,
#'                    delta_prob,
#'                    XX, XX_inv, Xy, diff_Phi, S) {
#'   #' This function computes the posterior density function g(delta, Y) for a given binary 
#'   #' vector delta in a Bayesian VAR model.
#'   #' 
#'   #' Parameters:
#'   #' - delta (vector): Current 0/1 vector of length k.
#'   #' - T (integer): Number of observations.
#'   #' - k (integer): Number of equations in the VAR model.
#'   #' - S_lowerbar (matrix): Prior scale part of the inverse-Wishart prior.
#'   #' - v_lowerbar (numeric): Degrees of freedom in the inverse-Wishart prior.
#'   #' - tau0, tau1 (scalars): Hyperparameters controlling the spike and slab adjustment given delta.
#'   #' - se_ols (vectors): A vector length k, providing standard deviation errors of the OLS residuals
#'   #'  for the coefficients based on delta.
#'   #' - delta_prob (function or vector): Prior mass function for delta.
#'   #' - XX, XX_inv (matrices): Regressor matrix and its inverse.
#'   #' - Xy (matrix): Matrix product of X'Y.
#'   #' - diff_Phi (matrix): Difference between prior mean and estimated coefficients.
#'   #' - S (matrix): Residual variance matrix.
#'   #' 
#'   #' Returns: 
#'   #' - A scalar value representing g(delta, Y).
#'   
#'   # compute updatet prior matrices based on new given delta
#'   Omega_lowerbar_delta <- build_Omega_lowerbar(delta, tau0, tau1, se_ols)
#'   inv_Omega_lowerbar_delta <- chol2inv(chol(Omega_lowerbar_delta))
#'   Omega_upperbar_delta <- solve(inv_Omega_lowerbar_delta + XX)
#'   
#'   M <- Omega_lowerbar_delta + XX_inv
#'   M_inv <- chol2inv(chol(M))
#'   
#'   S_upperbar_delta <- S_lowerbar + S + t(diff_Phi) %*% M_inv %*% diff_Phi
#'   
#'   # compute determinant factors 
#'   Omega_lowerbar_chol <- chol(Omega_lowerbar_delta)
#'   log_det_Omega_lowerbar <- 2 * sum(log(diag(Omega_lowerbar_chol)))
#'   
#'   Omega_upperbar_chol <- chol(Omega_upperbar_delta)
#'   log_det_Omega_upperbar <- 2 * sum(log(diag(Omega_upperbar_chol)))
#'   
#'   S_upperbar_chol <- chol(S_upperbar_delta)
#'   log_det_S_upperbar <- 2 * sum(log(diag(S_upperbar_chol)))
#'   
#'   # compute likelihood in log scale
#'   val <- ((log_det_Omega_lowerbar - log_det_Omega_upperbar) * (-k/2)) +
#'     (log_det_S_upperbar * (-(v_lowerbar + T)/2)) +
#'     pi_delta(delta, delta_prob)
#'   
#'   return(val)
#' }
#' 
#' 
#' update_delta <- function(delta_current,
#'                          T, k,
#'                          S_lowerbar,
#'                          v_lowerbar,
#'                          tau0, tau1, se_ols,
#'                          delta_prob,
#'                          XX, XX_inv, Xy, diff_Phi, S) {
#'   #' This function sequentially updates the delta vector using a Gibbs sampling procedure in a Bayesian VAR model.
#'   #' 
#'   #' Parameters:
#'   #' - delta (vector): Current 0/1 vector of length k.
#'   #' - T (integer): Number of observations.
#'   #' - k (integer): Number of equations in the VAR model.
#'   #' - S_lowerbar (matrix): Prior scale part of the inverse-Wishart prior.
#'   #' - v_lowerbar (numeric): Degrees of freedom in the inverse-Wishart prior.
#'   #' - tau0, tau1 (scalars): Hyperparameters controlling the spike and slab adjustment given delta.
#'   #' - se_ols (vectors): A vector length k, providing standard deviation errors of the OLS residuals
#'   #'  for the coefficients based on delta.
#'   #' - delta_prob (function or vector): Prior mass function for delta.
#'   #' - XX, XX_inv (matrices): Regressor matrix and its inverse.
#'   #' - Xy (matrix): Matrix product of X'Y.
#'   #' - diff_Phi (matrix): Difference between prior mean and estimated coefficients.
#'   #' - S (matrix): Residual variance matrix.
#'   #' 
#'   #' Returns:
#'   #' - An updated binary vector delta_new after the Gibbs sampling step.
#'   
#'   m <- length(delta_current)
#'   
#'   # We will update in sequence i=1,...,m.
#'   delta_new <- delta_current
#'   
#'   for (i in seq_len(m)) {
#'     
#'     # 1) Temporarily force delta_i = 0
#'     delta_zero <- delta_new
#'     delta_zero[i] <- 0
#'     u_0 <- calc_g(delta_zero, T, k,
#'                   S_lowerbar,
#'                   v_lowerbar,
#'                   tau0, tau1, se_ols,
#'                   delta_prob,
#'                   XX, XX_inv, Xy, diff_Phi, S)
#'     
#'     # 2) Temporarily force delta_i = 1
#'     delta_one <- delta_new
#'     delta_one[i] <- 1
#'     u_1 <- calc_g(delta_one, T, k,
#'                   S_lowerbar,
#'                   v_lowerbar,
#'                   tau0, tau1, se_ols,
#'                   delta_prob,
#'                   XX, XX_inv, Xy, diff_Phi, S)
#'     
#'     # 3) Compute Bernoulli probability
#'     lmax <- max(u_1, u_0)
#'     u_1_log_diff <- u_1 - lmax
#'     u_0_log_diff <- u_0 - lmax
#'     u_1_exp_diff <- exp(u_1_log_diff)
#'     u_0_exp_diff <- exp(u_0_log_diff)
#'     p_i <- u_1_exp_diff / (u_1_exp_diff + u_0_exp_diff)
#'     
#'     # 4) Draw delta_i ~ Bernoulli(p_i)
#'     delta_new[i] <- rbinom(n=1, size=1, prob=p_i)
#'   }
#'   
#'   return(delta_new)
#' }
