# test functions


# ---- import libraries --------------------------------------------------------

source(paste0(getwd(),"/lib.r"))


# ---- import functions --------------------------------------------------------

function_files <- c(
  "_distr_samplers.R",
  "_forecast.R",
  "_helper.R",
  "_mcmc.R",
  "_minnesota_prior.R",
  "bvar.R",
  "simulate_data.R"
)
for (i in 1:length(function_files)) {
  source(paste0(getwd(),"/functions/", function_files[i]))
}


# ---- simulate data -----------------------------------------------------------

k <- 7
p <- 3
h <- 4
T_data <- 40 + h
data <- sim_y(k, p, T_data, 
              min_indiv_shocks = 0, 
              max_indiv_shocks = 0, 
              min_high_vol_periods = 0, 
              max_high_vol_periods = 0,
              min_exog_shocks = 0,
              max_exog_shocks = 0)
data <- data*100
Ypred <- data[(T_data-h+1):T_data,]
Yraw <- data[1:(T_data-h),]

# ---- run Minnesota Gibbs Sampling --------------------------------------------

dummy_pars_ex <- list(delta = 1.0,
                      gamma = 1.0,
                      prior_mean = NULL)
intercept = FALSE

start <- Sys.time()
res_gibbs <- run_bvar(
  Yraw         = Yraw,
  p            = p,
  intercept    = intercept,
  use_dummies  = FALSE,
  dummy_pars   = dummy_pars_ex,
  use_flat     = FALSE,
  lag_mean     = 1,
  pi1          = 0.2,
  pi3          = 1,
  pi4          = 1000,
  sigma_vec    = NULL, 
  n_draws      = 5000,
  burnin       = 1000
)
mcmc <- Sys.time() - start
mcmc

# ---- inspect outputs ---------------------------------------------------------

Phi_post_median   <- apply(res_gibbs$Phi, c(1,2), median)
Sigma_post_median <- apply(res_gibbs$Sigma, c(1,2), median)

# give names to results
varNames <- paste0("Y", 1:k)
rowNamesPhi <- c("Intercept")
for (i in 1:p) {
  rowNamesPhi <- c(rowNamesPhi, paste0("Lag", rep(i, k)))
}
if(!intercept) rowNamesPhi <- rowNamesPhi[-1]

dimnames(Phi_post_median) <- list(rowNamesPhi, varNames)
dimnames(Sigma_post_median) <- list(varNames, varNames)

# view results
Phi_post_median
Sigma_post_median


# ---- predict data ------------------------------------------------------------

# draw from posterior predictive distribution

start <- Sys.time()
Y_forecasts <- predict_bvar(
  Phi_store   = res_gibbs$Phi,
  Sigma_store = res_gibbs$Sigma,
  Yraw        = Yraw,
  p           = p,
  H           = h,
  draw_shocks = TRUE,
  intercept = intercept,
  n_cores = 2
)
pred <- Sys.time() - start

((pred + mcmc) * 10000) / (60*60)

# inspect results
alpha <- 0.05

Y_forecast_median <- apply(Y_forecasts, c(1, 2), median)
Y_forecast_ci_lower <- apply(Y_forecasts, c(1, 2), quantile, probs = alpha/2)
Y_forecast_ci_upper <- apply(Y_forecasts, c(1, 2), quantile, probs = (1-alpha/2))

Y_forecast_median
Y_forecast_ci_lower
Y_forecast_ci_upper


# ---- evaluate prediction -----------------------------------------------------

# --- root mean squared error (rmse)
# individual forecasts
pred_error_mat <- matrix(NA, nrow(Y_forecast_median), ncol(Y_forecast_median), 
                   dimnames = list(rownames(Y_forecast_median), paste0("Y", 1:k)))

for (row in 1:nrow(pred_error_mat)) {
  for (col in 1:ncol(pred_error_mat)) {
    pred_error_mat[row, col] <- sqrt(mean((Y_forecast_median[row, col] - Ypred[row, col])^2))
  }
}

# h-step ahead forecast
row_rmse <- matrix(rep(NA, h), ncol = 1, dimnames = list(paste0("T+", 1:h), "Y1 ... Yk"))
for (i in 1:h) {
  row_rmse[i,] <- sqrt(mean((Y_forecast_median[i,] - Ypred[i,])^2))
}

# variable forecast
col_rmse <- matrix(rep(NA, k), nrow = 1, dimnames = list("T1 ... Th", paste0("Y", 1:k)))
for (i in 1:k) {
  col_rmse[,i] <- sqrt(mean((Y_forecast_median[,i] - Ypred[,i])^2))
}

# total rmse
all_rmse <- sqrt(mean((Y_forecast_median - Ypred)^2))


pred_error_mat
row_rmse
col_rmse
all_rmse
