# functions which sample from different distributions


posterior_draw_Phi <- function(M1, Sigma, V1) {
  #' Draws posterior samples of the BVAR coefficient matrix (Phi).
  #' 
  #' Parameters:
  #' - M1 (matrix): Posterior mean matrix of dimensions m x k.
  #' - Sigma (matrix): Covariance matrix of dimensions k x k.
  #' - V1 (matrix): Posterior covariance matrix of dimensions m x m.
  #' 
  #' Returns:
  #' - A matrix Phi_draw sampled from the posterior distribution of the BVAR coefficients, with dimensions (m x k).

  m <- nrow(M1)
  k <- ncol(M1)
  E <- chol(V1) %*% matrix(rnorm(m*k), m, k) %*% chol(Sigma)
  return(M1 + E)
}


posterior_draw_Sigma <- function(S1, nu1) {
  #' Draws posterior samples of the error covariance matrix (Sigma).
  #' 
  #' Parameters:
  #' - S1 (matrix): Posterior scale matrix for Sigma (k x k).
  #' - nu1 (numeric): Posterior degrees of freedom for Sigma.
  #' 
  #' Returns:
  #' - A matrix Sigma_draw sampled from the posterior inverse-Wishart distribution, with dimensions (k x k).

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


rmvn_proposal <- function(n, mean, sigma) {
  #' Generates random samples from a multivariate normal distribution.
  #'
  #' Parameters:
  #' - n (integer): The number of samples to generate.
  #' - mean (numeric vector): The mean vector of the multivariate normal distribution, of length `m`.
  #' - sigma (list): A list containing the eigen decomposition of the covariance matrix, with:
  #'   - values (numeric vector): The eigenvalues of length m.
  #'   - vectors (matrix): The eigenvectors matrix of dimensions m x m.
  #'
  #' Returns:
  #' - A matrix of dimension n x m where each row is a sampled vector from the specified 
  #' multivariate normal distribution.
  
  m <- length(sigma[["values"]])
  R <- t(sigma[["vectors"]] %*%
           (t(sigma[["vectors"]]) * sqrt(sigma[["values"]])))
  out <- matrix(rnorm(n * m), nrow = n, ncol = m, byrow = TRUE) %*% R
  out <- sweep(out, 2, mean, "+")
  return(out)
}


helper_make_positive_definite <- function(matrix, lambda = 1e-8, max_iter = 1000) {
  #' Ensures a matrix is positive definite.
  #' 
  #' Parameters:
  #' - matrix (matrix): Matrix to be checked and adjusted.
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
