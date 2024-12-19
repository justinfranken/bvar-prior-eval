
# Set parameters
T <- 50  # number of time periods
K <- 3    # number of firms (variables)
p <- 2    # lag order

# set.seed(123) # for reproducibility

#---------------------------------------------------------
# Step 1: Construct stable VAR coefficients
# Each A_i is K x K. We choose small coefficients to ensure stability.
# Construct A_i matrices with domain structure:
# - Diagonal elements: moderate positive persistence (e.g. around 0.2)
# - Off-diagonal elements: smaller positive values (e.g. around 0.05) to represent some cross-influence
# After constructing, we will verify stability and scale if needed.

A_list <- vector("list", p)
for (i in 1:p) {
  A_raw <- matrix(0, nrow=K, ncol=K)
  
  # Diagonal: persistence
  diag(A_raw) <- runif(K, min=0.2, max=0.4)  # Own-lag coefficients in [0.1,0.3]
  
  # Off-diagonal: cross-effects
  # Small random positive values to indicate that if one firm grows, it slightly helps others in the future.
  for (r in 1:K) {
    for (col in 1:K) {
      if (r != col) {
        A_raw[r, col] <- rnorm(1, mean=0.02, sd=0.01) # Slight positive influence
      }
    }
  }
  
  # Check eigenvalues for stability: If max eigenvalue >= 1, scale down
  eigenvals <- eigen(A_raw)$values
  if (max(Mod(eigenvals)) >= 1) {
    scale_factor <- 1.1 * max(Mod(eigenvals))
    A_raw <- A_raw / scale_factor
  }
  
  A_list[[i]] <- A_raw
}

#---------------------------------------------------------
# Step 2: Set a drift (constant term)
# Assume all firms have a slightly positive average growth rate, say around 1% per period.
c <- rep(0.01, K)  # all firms grow ~1% on average

#---------------------------------------------------------
# Step 3: Generate correlated disturbance terms
# Firms are in the same sector, so their shocks are positively correlated.
Rho <- matrix(c(1,   0.7, 0.6,
                0.7, 1,   0.65,
                0.6, 0.65,1   ), nrow=K)
# Set standard deviations: assume the leading firm has slightly lower volatility
sd_vec <- c(0.02, 0.025, 0.025)
Sigma <- diag(sd_vec) %*% Rho %*% diag(sd_vec)

# We must generate normal random vectors with covariance Sigma.
# Without external libraries, we can use a Cholesky decomposition:
cholSigma <- chol(Sigma)
Z <- matrix(rnorm(T*K), T, K)
eps <- Z %*% cholSigma  # eps now has covariance Sigma

#---------------------------------------------------------
# Step 4: Simulate the VAR(p) process
Y <- matrix(0, nrow=T, ncol=K)

# Initialize first p observations:
# Start them around the drift, perhaps with small random deviations
for (i in 1:p) {
  Y[i,] <- c + rnorm(K, mean=0, sd=0.01)
}

# Generate the series
for (t in (p+1):T) {
  # Compute the VAR component
  Y_t <- c
  for (lag in 1:p) {
    Y_t <- Y_t + A_list[[lag]] %*% Y[t-lag,]
  }
  # Add the shock
  Y[t,] <- Y_t + eps[t,]
}

#---------------------------------------------------------
# Step 5: Inspect the simulated results
# Each column is a firm. They should exhibit similar growth patterns and co-movements.
ts.plot(Y, col=1:K, main="Simulated Growth Rates of Comparable Firms")
legend("topright", legend=paste("Firm", 1:K), col=1:K, lty=1)

# The time series Y now represent growth rates for three comparable firms.
# They share similar drift, are correlated, and follow a stable VAR(p) process.












