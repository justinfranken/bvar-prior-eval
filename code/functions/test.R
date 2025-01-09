##############################################################################
# 0) Setup
##############################################################################
# Load packages
library(MASS)    # for mvrnorm
library(Matrix)  # for solve, etc.
library(parallel)

##############################################################################
# 1) Example Data
##############################################################################


k      <- 3     # number of variables
p      <- 3     # lag order
T_data <- 105 + p   # number of observations

# True coefficients (Phi_true) in companion form
# For simplicity, let's define a random stable VAR(3)
# or just define some manual matrix.
Phi_true_1 <- matrix(c(0.5, 0.1, 0.0,
                       0.0, 0.6, 0.1,
                       0.1, 0.0, 0.4), nrow=k, byrow=TRUE)
Phi_true_2 <- 0.1 * diag(k)
Phi_true_3 <- matrix(0, nrow=k, ncol=k)  # third lag is zero, say

# True Sigma
Sigma_true <- matrix(c(1.0, 0.2, 0.0,
                       0.2, 0.8, 0.1,
                       0.0, 0.1, 1.2), k, k)

# Simulate data
Yraw <- matrix(0, nrow=T_data, ncol=k)

# Generate shocks
E <- mvrnorm(T_data, mu=rep(0,k), Sigma=Sigma_true)

# Start from zeros
for(t in (p+1):T_data){
  Yraw[t, ] <- (Phi_true_1 %*% Yraw[t-1,] +
                  Phi_true_2 %*% Yraw[t-2,] +
                  Phi_true_3 %*% Yraw[t-3,]) + E[t,]
}
Yraw <- Yraw[-(1:p),]

#k <- 2
#p <- 2
#T_data <- 105
#Yraw <- sim_y(k, p, T_data, 
#              min_indiv_shocks = 0, 
#              max_indiv_shocks = 0, 
#              min_high_vol_periods = 0, 
#              max_high_vol_periods = 0,
#              min_exog_shocks = 0,
#              max_exog_shocks = 0)



Ypred <- Yraw[101:105,]
Yraw <- Yraw[1:100,]

##############################################################################
# 2) Create Lagged Data (Y, X)
##############################################################################


create_lagged_data <- function(data, p, intercept = FALSE) {
  # data:      (T x k) matrix of observations
  # p:         number of lags
  # intercept: logical; if TRUE, include a leading column of 1's in X
  #
  # Output:
  #    y: (T - p) x k     -- rows are time, columns are variables
  #    X: (T - p) x [k*p  or (1 + k*p)]   -- design matrix, depending on intercept
  
  T <- nrow(data)
  k <- ncol(data)
  
  # embed(...) produces [T - p] rows, each containing [1 + p] "blocks" of length k
  # i.e. row i = (Y_i, Y_{i+1}, ..., Y_{i+p}) if p+1 steps
  # so 'lagged_data' is (T-p) x [k * (p+1)] dimension
  lagged_data <- embed(data, p + 1)  # each row: Y_t, Y_{t-1}, ..., Y_{t-p}
  
  # The "current" Y_t is always the first k columns
  y <- lagged_data[, 1:k, drop=FALSE]  # => (T-p) x k
  
  # The lagged values are the remaining columns
  # i.e. columns (k+1) to end = (Y_{t-1}, ..., Y_{t-p})
  lag_cols <- lagged_data[, (k + 1):ncol(lagged_data), drop=FALSE]
  
  if (intercept) {
    # Build design matrix: [Intercept, Y_{t-1}, ..., Y_{t-p}]
    X <- cbind(1, lag_cols)
  } else {
    # Build design matrix: [Y_{t-1}, ..., Y_{t-p}] only
    X <- lag_cols
  }
  
  return(list(y = y, X = X))
}


create_dummy_observations_glp <- function(Y, p, intercept = FALSE,
                                          delta = 1.0,  # tightness for initial observational prior
                                          gamma = 1.0,  # tightness for sum-of-coefficients prior
                                          prior_mean = NULL) {
  T_real <- nrow(Y)
  k <- ncol(Y)
  
  if (is.null(prior_mean)) {
    prior_mean <- colMeans(Y[1:p, , drop = FALSE], na.rm = TRUE)
  }
  
  # Sum-of-Coefficients Dummy
  Y_sum <- diag((gamma)^-1 * prior_mean)
  X_sum <- matrix(0, nrow = k, ncol = p * k)
  for (i in 1:p) {
    X_sum[, ((i - 1) * k + 1):(i * k)] <- Y_sum
  }
  if (intercept) X_sum <- cbind(rep(0, nrow(X_sum)), X_sum)
  
  # Dummy Initial Observational Prior
  Y_init <- matrix((delta)^-1 * prior_mean, nrow = 1)
  X_init <- matrix(0, nrow = 1, ncol = p * k + as.integer(intercept))
  if (intercept) X_init[1, 1] <- 1 / (delta)^-1
  X_init[1, (1+as.integer(intercept)):(as.integer(intercept) + p * k)] <- (delta)^-1 * rep(prior_mean, times = p)
  
  # Combine Dummies
  Y_dum <- rbind(Y_sum, Y_init)
  X_dum <- rbind(X_sum, X_init)
  
  return(list(Y_dum = Y_dum, X_dum = X_dum))
}

create_lagged_data_with_dummies <- function(Y, p, intercept = FALSE,
                                            use_dummies = FALSE,
                                            dummy_pars = list()) {
  # 1) Real-data-lagged
  ld <- create_lagged_data(Y, p, intercept = intercept)
  Y_real <- ld$y
  X_real <- ld$X
  
  if(!use_dummies) {
    return(list(Y = Y_real, X = X_real))
  } else {
    # 2) Build dummy blocks
    dums <- do.call(create_dummy_observations_glp,
                    args = c(list(Y=Y, p=p, intercept=intercept),
                             dummy_pars))
    Y_dum <- dums$Y_dum
    X_dum <- dums$X_dum
    
    # 3) Append
    Y_aug <- rbind(Y_dum, Y_real)
    X_aug <- rbind(X_dum, X_real)
    
    return(list(Y = Y_aug, X = X_aug))
  }
}


##############################################################################
# 3) Constructing a Minnesota(-like) Prior
##############################################################################

build_minnesota_M0_V0 <- function(k, p, intercept=FALSE,
                                  lambda=0.2, alpha=1.0,
                                  sigma_vec=NULL) {
  # k: number of variables
  # p: number of lags
  # intercept: whether to include an intercept
  # lambda: overall tightness
  # alpha : lag decay exponent
  # sigma_vec: prior scale for each variable (length k).
  #            If NULL, set them all to 1 for demonstration.
  
  if(is.null(sigma_vec)) {
    sigma_vec <- rep(1, k)
  }
  
  m <- if(intercept) 1 + k*p else k*p
  
  # M0 is k x m
  M0 <- matrix(0, nrow=m, ncol=k)
  if(!intercept){
    M0[1:k, 1:k] = diag(k)
  } else{
    M0[2:(k+1), 1:k] = diag(k)
  }
  
  # V0 is m x m
  V0 <- matrix(0, m, m)
  
  # 1) Intercept row => large variance
  if(intercept){
    V0[1,1] <- 10^6
  }
  
  # 2) Fill the lag coefficients
  for(ell in 1:p){
    for(j in 1:k){
      # column index for (lag ell, var j)
      col_index <- (if(intercept) 1 else 0) + (ell-1)*k + j
      for(i in 1:k){
        if(i == j){
          # own-lag
          prior_var_ij <- (lambda^2) / (ell^alpha)
        } else {
          # cross-lag
          prior_var_ij <- (lambda^2) / (ell^alpha) * (sigma_vec[i]^2 / sigma_vec[j]^2)
        }
        V0[col_index, col_index] <- prior_var_ij
      }
    }
  }
  
  return(list(M0 = M0, V0 = V0))
}


prepare_prior <- function(Y, p, intercept = FALSE,
                          use_dummies = FALSE,
                          dummy_pars = list(),
                          lambda = 0.2, alpha = 1.0,
                          sigma_vec = NULL) {
  k <- ncol(Y)
  
  if(!use_dummies) {
    # Use standard Minnesota prior
    mn_pr <- build_minnesota_M0_V0(k, p, intercept=intercept,
                                   lambda=lambda, alpha=alpha,
                                   sigma_vec=sigma_vec)
    M0 <- mn_pr$M0
    V0 <- mn_pr$V0
    return(list(M0 = M0, V0 = V0))
  } else {
    # Use GLP Dummy Observations with flat prior
    # No need to build M0 and V0 from build_minnesota_M0_V0
    # Instead, set M0 and V0 to be non-informative
    m <- if(intercept) 1 + k*p else k*p
    M0_flat <- matrix(0, nrow = m, ncol = k)   # Flat prior mean
    V0_flat <- diag(m) * 1e8                    # Flat prior covariance
    return(list(M0 = M0_flat, V0 = V0_flat))
  }
}


##############################################################################
# 4) Gibbs Sampler Setup
##############################################################################

# (A) Sample from Inverse-Wishart
riwish <- function(df, S){
  # draw one sample from IW(S, df).
  # We can do this by drawing from Wishart(S^-1, df) and inverting.
  W <- rWishart(1, df, solve(S))[,,1]
  return(solve(W))
}

# (B) Sample from Matrix-Normal
#     Phi ~ MN(M, U, V) => vec(Phi) ~ N(vec(M), V \otimes U)
rmatrixnorm <- function(M, U, V){
  k <- nrow(M)
  m <- ncol(M)
  vecM <- as.vector(M)
  SigmaMN <- kronecker(V, U)
  vecPhi <- mvrnorm(1, mu=vecM, Sigma=SigmaMN)
  Phi_draw <- matrix(vecPhi, nrow=k, ncol=m)
  return(Phi_draw)
}


##############################################################################
# 5) Posterior Update Functions
##############################################################################

posterior_draw_Phi <- function(Y, X, Sigma, M0, V0) {
  # Y: (T_eff x k)
  # X: (T_eff x m)
  # Sigma: (k x k), current draw of Sigma
  # M0: (k x m), prior mean
  # V0: (m x m), prior covariance of coefficients (right-cov in MN)

  # Posterior for Phi: MN(M1, Sigma, V1)
  # V1 = (V0^{-1} + X'X)^{-1}
  V1_inv <- solve(V0) + t(X) %*% X
  V1     <- solve(V1_inv)
  
  M1 <- V1 %*% (solve(V0) %*% M0 + t(X) %*% Y) # => (k x m)
  
  # Sample
  Phi_draw <- rmatrixnorm(M1, Sigma, V1)
  
  return(Phi_draw)
}

posterior_draw_Sigma <- function(Y, X, Phi, M0, V0, S0, nu0) {
  T_eff <- nrow(Y)
  k     <- ncol(Y)
  
  # Residuals
  # Phi is (k x m), X is (T_eff x m), so X %*% t(Phi) is (T_eff x k)
  E <- Y - X %*% Phi
  SSE <- t(E) %*% E  # (k x k)
  
  # (Phi - M0) part
  diffPhi <- M0 - Phi  # (k x m)
  
  # We want: diffPhi * V0^{-1} * diffPhi'
  # But carefully:  V0^{-1} is (m x m), diffPhi is (k x m), so
  #   diffPhi %*% V0^{-1} => (k x m), then multiply by diffPhi' => (k x k)
  inv_middle_term <- solve(V0 + solve(t(X) %*% X))
  second_term <- t(diffPhi) %*% inv_middle_term %*% diffPhi
  
  # Posterior dof
  nu1 <- nu0 + T_eff
  # Posterior scale
  S1  <- S0 + SSE + second_term
  
  # Draw Sigma from IW(S1, nu1)
  Sigma_draw <- riwish(nu1, S1)
  return(Sigma_draw)
}


##############################################################################
# 6) Gibbs Sampler
##############################################################################

run_minnesota_IW_bvar <- function(Y, X, M0, V0, S0, nu0,
                                  n_draws=5000, burnin=1000) {
  # Y: (T_eff x k)
  # X: (T_eff x m)
  # M0, V0, S0, nu0 define the Normal-Wishart prior
  # n_draws: total MCMC draws
  # burnin : how many to discard as burn-in
  
  T_eff <- nrow(Y)
  k     <- ncol(Y)
  m     <- ncol(X)
  
  # Storage
  Phi_store   <- array(NA, dim=c(m, k, n_draws))
  Sigma_store <- array(NA, dim=c(k, k, n_draws))
  
  # Initialize
  Sigma_current <- diag(apply(Y, 2, var)) # or set to sample from IW once
  
  for(iter in 1:n_draws) {
    # Step 1: draw Phi from its conditional
    Phi_current <- posterior_draw_Phi(Y, X, Sigma_current, M0, V0)
    
    # Step 2: draw Sigma from its conditional
    Sigma_current <- posterior_draw_Sigma(Y, X, Phi_current, M0, V0, S0, nu0)
    
    # Store
    Phi_store[,,iter]   <- Phi_current
    Sigma_store[,,iter] <- Sigma_current
  }
  
  # Drop burn-in
  keep_idx <- (burnin+1):n_draws
  Phi_store   <- Phi_store[,, keep_idx]
  Sigma_store <- Sigma_store[,, keep_idx]
  
  return(list(Phi=Phi_store, Sigma=Sigma_store))
}


run_bvar <- function(Yraw, p, intercept = FALSE,
                     use_dummies = FALSE,
                     dummy_pars = list(delta = 1.0, gamma = 1.0, prior_mean = NULL),
                     lambda = 0.2, alpha = 1.0,
                     sigma_vec = NULL,
                     S0 = diag(ncol(Yraw)), nu0 = ncol(Yraw) + 2,
                     n_draws = 5000, burnin = 1000) {
  # Step 1: Build (Y, X), possibly with dummies
  ld_aug <- create_lagged_data_with_dummies(Yraw, p, intercept,
                                            use_dummies,
                                            dummy_pars)
  Y <- ld_aug$Y
  X <- ld_aug$X
  
  # Step 2: Prepare prior
  prior <- prepare_prior(Yraw, p, intercept,
                         use_dummies,
                         dummy_pars,
                         lambda, alpha,
                         sigma_vec)
  M0 <- prior$M0
  V0 <- prior$V0
  
  # Step 3: If using dummies, ensure the prior is flat
  if(use_dummies){
    # Optionally, you can adjust S0 and nu0 here if needed
    # For simplicity, we'll keep them as passed
    # You might want to set S0 to a small value if using dummies as data
    S0 <- diag(ncol(Yraw)) * 1e-3
    nu0 <- 0
  }
  
  # Step 4: Run the sampler
  res_gibbs <- run_minnesota_IW_bvar(
    Y, X,
    M0, V0,
    S0, nu0,
    n_draws = n_draws,
    burnin = burnin
  )
  
  return(res_gibbs)
}

##############################################################################
# 7) Run sampler
##############################################################################

# Use the prior from above:
# M0 (k x m), V0 (m x m)
# S0 (k x k), nu0 (scalar)
n_draws <- 5000
burnin  <- 1000


use_dummies_ex <- TRUE  # Set to FALSE to use standard Minnesota prior

# Define dummy parameters if using dummies
dummy_pars_ex <- list(delta = 1.0,
                      gamma = 1.0,
                      prior_mean = NULL)  # or your own chosen mean vector


res_gibbs <- run_bvar(
  Yraw         = Yraw,
  p            = p,
  intercept    = FALSE,
  use_dummies  = use_dummies_ex,
  dummy_pars   = dummy_pars_ex,
  lambda       = 0.2,
  alpha        = 1.0,
  sigma_vec    = NULL,
  S0           = diag(k),
  nu0          = k + 2,
  n_draws      = n_draws,
  burnin       = burnin
)



# Posterior means
Phi_post_mean   <- apply(res_gibbs$Phi, c(1,2), mean)
Sigma_post_mean <- apply(res_gibbs$Sigma, c(1,2), mean)


# Suppose:
varNames <- paste0("Y", 1:k)

# Build colnames for Phi
rowNamesPhi <- c("Intercept")
for (i in 1:k) {
  rowNamesPhi <- c(rowNamesPhi, paste0("Lag", rep(i, k)))
}

intercept <- FALSE
if(!intercept) rowNamesPhi <- rowNamesPhi[-1]

# Assign row & column names to Phi
dimnames(Phi_post_mean) <- list(rowNamesPhi, varNames)

# Assign row & column names to Sigma
dimnames(Sigma_post_mean) <- list(varNames, varNames)

# Now let's see them
Phi_post_mean
Sigma_post_mean



##############################################################################
# 8) Forecasting
##############################################################################

reshape_phi <- function(Phi_draw, k, p) {
  # Phi_draw has dimension (p*k) x k
  # We want a list of p matrices each k x k
  Phi_list <- vector("list", length = p)
  for (ell in seq_len(p)) {
    # rows for lag ell go from ( (ell-1)*k + 1 ) to ( ell*k )
    row_start <- (ell - 1) * k + 1
    row_end   <- ell * k
    Phi_list[[ell]] <- Phi_draw[row_start:row_end, , drop = FALSE]
  }
  return(Phi_list)
}

forecast_from_draw <- function(Phi_draw_init, Sigma_draw,
                               Y_last_p,  # A matrix or list of the last p observations, each is 1 x k
                               H = 5,
                               draw_shocks = TRUE,
                               intercept = FALSE) {
  # Y_last_p should be in ascending order in time,
  # e.g. Y_last_p[1, ] = Y_{T-p+1}, ..., Y_last_p[p, ] = Y_{T}.
  
  k <- ncol(Y_last_p)
  p <- nrow(Phi_draw_init) / k  # we assume dimension: (p*k) x k
  
  # Reshape Phi_draw into a list of (k x k) blocks and separate intercept
  Phi_intercept <- rep(0, k)
  if(intercept){
    Phi_intercept <- Phi_draw_init[1,]
    Phi_draw <- Phi_draw_init[-1,]
  } else{
    Phi_draw <- Phi_draw_init
  } 
  Phi_list <- reshape_phi(Phi_draw, k, p)
  
  # Storage for forecasts
  Y_future <- matrix(NA, nrow = H, ncol = k,
                     dimnames = list(paste0("T+", 1:H), colnames(Y_last_p)))
  
  # We keep a rolling buffer of the last p states for iteration
  # Start with the user-supplied last p actual data
  Y_buffer <- Y_last_p
  
  for (h in seq_len(H)) {
    # Compute mean prediction
    Y_next_mean <- numeric(k)
    for (ell in seq_len(p)) {
      # The row we need for lag ell is the row from the buffer:
      #   the buffer has last p observations with the newest at the bottom
      #   so lag ell is Y_buffer[p - (ell - 1), ]
      Y_next_mean <- Y_next_mean + Phi_list[[ell]] %*% 
        Y_buffer[p - (ell - 1),]
    }
    
    # If you want a draw (to incorporate predictive uncertainty):
    # Y_next = Y_next_mean + e_t,   e_t ~ N(0, Sigma_draw)
    # Otherwise, for the conditional mean forecast only:
    if (draw_shocks) {
      e_t <- mvrnorm(n = 1, mu = rep(0, k), Sigma = Sigma_draw)
    } else {
      e_t <- rep(0, k)
    }
    
    Y_next <- as.numeric(Phi_intercept + Y_next_mean + e_t)
    
    # Store
    Y_future[h, ] <- Y_next
    
    # Update buffer: shift up, append Y_next
    Y_buffer <- rbind(Y_buffer[-1, ], Y_next)
  }
  
  return(Y_future)
}


predict_bvar <- function(Phi_store, Sigma_store, 
                         Yraw, p, H = 5, 
                         draw_shocks = TRUE,
                         intercept = FALSE) {
  # Phi_store: array of dimension [m, k, n_post_draws]
  # Sigma_store: array of dimension [k, k, n_post_draws]
  # Yraw: full dataset (T x k)
  # p: lag order
  # H: forecast horizon
  # draw_shocks: whether to sample from Sigma for each step ahead
  
  k <- ncol(Yraw)
  T_data <- nrow(Yraw)
  
  # Last p observations from the training sample
  Y_last_p <- Yraw[(T_data - p + 1):T_data, , drop = FALSE]
  
  n_post_draws <- dim(Phi_store)[3]
  
  # We'll store forecasts for each draw in a 3D array: (H x k x n_post_draws)
  Y_forecasts <- array(NA, dim = c(H, k, n_post_draws),
                       dimnames = list(paste0("T+", 1:H),
                                       colnames(Yraw),
                                       NULL))
  
  for (i in seq_len(n_post_draws)) {
    Phi_draw_i   <- Phi_store[ , , i]
    Sigma_draw_i <- Sigma_store[ , , i]
    
    Y_forecasts[ , , i] <- forecast_from_draw(
      Phi_draw_init = Phi_draw_i,
      Sigma_draw = Sigma_draw_i,
      Y_last_p   = Y_last_p,
      H          = H,
      draw_shocks = draw_shocks,
      intercept = intercept
    )
  }
  
  return(Y_forecasts)
}



# Posterior draws:
Phi_store   <- res_gibbs$Phi   # dimension: (m, k, #draws) 
Sigma_store <- res_gibbs$Sigma # dimension: (k, k, #draws) 


Y_forecasts <- predict_bvar(
  Phi_store   = Phi_store,
  Sigma_store = Sigma_store,
  Yraw        = Yraw,
  p           = p,
  H           = 5, 
  draw_shocks = TRUE,
  intercept = FALSE
)

# Y_forecasts is (5 x k x n_post_draws)
# Summaries:
Y_forecast_mean <- apply(Y_forecasts, c(1, 2), mean)
Y_forecast_ci_lower <- apply(Y_forecasts, c(1, 2), quantile, probs = 0.025)
Y_forecast_ci_upper <- apply(Y_forecasts, c(1, 2), quantile, probs = 0.975)

Y_forecast_mean
Y_forecast_ci_lower
Y_forecast_ci_upper


