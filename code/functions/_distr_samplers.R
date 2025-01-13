# functions which sample from different distributions


posterior_draw_Phi <- function(M1, Sigma, V1) {
  #' Draws posterior samples of the BVAR coefficient matrix (Phi).
  #' 
  #' Parameters:
  #' - M1 (matrix): Updated mean matrix of dimensions m x k.
  #' - Sigma (matrix): Covariance matrix of dimensions k x k.
  #' - V1 (matrix): Prior updated covariance matrix of dimensions m x m.
  #' 
  #' Returns:
  #' - A matrix Phi_draw sampled from the posterior distribution of the BVAR coefficients, with dimensions (m x k).

  m <- nrow(M1)
  k <- ncol(M1)
  E <- chol(V1) %*% matrix(rnorm(m*k), m, k) %*% chol(Sigma)
  return(M1 + E)
}


posterior_draw_Sigma <- function(Y, X, Phi, M0, S0, nu0, inv_middle_term) {
  #' Draws posterior samples of the error covariance matrix (Sigma).
  #' 
  #' Parameters:
  #' - Y (matrix): Response matrix with dimensions (T_eff x k), where T_eff is the effective number of observations, and k is the number of variables.
  #' - X (matrix): Design matrix with dimensions (T_eff x m), where m is the number of predictors.
  #' - Phi (matrix): Current draw of the BVAR coefficient matrix, with dimensions (m x k).
  #' - M0 (matrix): Prior mean matrix for Phi (m x k).
  #' - S0 (matrix): Prior scale matrix for Sigma (k x k).
  #' - nu0 (numeric): Prior degrees of freedom for Sigma.
  #' - inv_middle_term(matrix): Precomputed term.
  #' 
  #' Returns:
  #' - A matrix Sigma_draw sampled from the posterior inverse-Wishart distribution, with dimensions (k x k).

  E <- Y - X %*% Phi
  SSE <- crossprod(E)
  diffPhi <- M0 - Phi
  
  second_term <- crossprod(diffPhi, inv_middle_term %*% diffPhi)
  
  # posterior dof
  nu1 <- nu0 + nrow(Y)
  # posterior scale
  S1  <- S0 + SSE + second_term
  
  # draw Sigma from IW(S1, nu1)
  Sigma_draw <- riwish(nu1, S1)
  return(helper_make_positive_definite(Sigma_draw))
}


riwish <- function(df, S){
  #' This function generates a single sample from the inverse-Wishart distribution.
  #' 
  #' Parameters:
  #' - df (numeric): Degrees of freedom for the inverse-Wishart distribution.
  #' - S (matrix): Scale matrix for the inverse-Wishart distribution. 
  #' If S is not positive definite, a small adjustment is made to ensure stability.
  #' 
  #' Returns:
  #' - A matrix representing a random draw from the inverse-Wishart distribution IW(S, df).
  
  helper_make_positive_definite(S)
  Sinv <- chol2inv(chol(S))
  W <- rWishart(1, df, Sinv)[,,1]
  return(chol2inv(chol(W)))
}


helper_make_positive_definite <- function(matrix, lambda = 1e-8, max_iter = 1000) {
  #' Ensures a matrix is positive definite.
  #' 
  #' Parameters:
  #' - matrix (matrix): Sigma matrix to be checked and adjusted.
  #' - lambda (numeric): The initial increment to be added to the diagonal elements (default = 1e-8).
  #' - max_iter (integer): Maximum number of iterations to attempt adjustment (default = 1000).
  #' 
  #' Returns:
  #' - A positive definite matrix obtained by adjusting the input matrix.
  #' 
  #' Throws:
  #' - An error if the matrix cannot be adjusted to be positive definite within the specified number of iterations.
  
  
  # check if the matrix is already positive definite
  is_positive_definite <- tryCatch({
    chol(matrix)
    TRUE
  }, error = function(e) {
    FALSE
  })
  
  # return the original matrix if already positive definite
  if (is_positive_definite) {
    return(matrix)  
  }
  
  # add increment to diagonal
  for (i in 1:max_iter) {
    result <- tryCatch({
      chol(matrix + diag(lambda, nrow(matrix)))
      return(matrix + diag(lambda, nrow(matrix)))
    }, error = function(e) {
      NULL
    })
    
    # stop loop if positive definite
    if (!is.null(result)) {
      return(result)
    }
    
    lambda <- lambda * 10
  }
  stop("Sigma draw could not be adjusted to positive definite matrix.\n 
       Adjust prior variance parameters or adjust data accordingly.")
}
