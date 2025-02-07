# test functions


# ---- import libraries --------------------------------------------------------

source(paste0(getwd(),"/lib.r"))


# ---- import functions --------------------------------------------------------
function_files <- c(
  "_distr_samplers.R",
  "_forecast.R",
  "_hierarch_fun.R",
  "_data_prep.R",
  "_mcmc.R",
  "_ssvs.R",
  "_prior_specification.R",
  "bvar.R",
  "simulate_data.R"
)
for (i in 1:length(function_files)) {
  source(paste0(getwd(),"/functions/", function_files[i]))
}
rm(i)
rm(function_files)


# ---- simulate data -----------------------------------------------------------
k <- 7
p <- 2
h <- 4
T_data <- 40 + h
data <- sim_y(k, p, T_data, 
              min_indiv_shocks = 0, 
              max_indiv_shocks = 0, 
              min_high_vol_periods = 0, 
              max_high_vol_periods = 0,
              min_exog_shocks = 0,
              max_exog_shocks = 0)
Ypred <- data[(T_data-h+1):T_data,]
Yraw <- data[1:(T_data-h),]


# ---- run Minnesota Gibbs Sampling --------------------------------------------
dummy_pars_ex <- list(mu = 1.0,
                      gamma = 1.0,
                      prior_mean = NULL)
intercept = TRUE
p_bvar <- 4

start <- Sys.time()
res_gibbs <- run_bvar_minnesota(
  Yraw         = Yraw,
  p            = p_bvar,
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
Sys.time() - start
Y_forecast_gibbs_median <- bvar_forecast(Yraw, 
                                         res_gibbs$Phi, 
                                         res_gibbs$Sigma, 
                                         p = p_bvar, h = 4, 
                                         alpha = 0.05, 
                                         intercept = intercept,
                                         draw_shocks = TRUE)$median
sqrt(mean((Y_forecast_gibbs_median - Ypred)^2))

# ---- run Metropolis Hastings sampler -----------------------------------------
data <- sim_y(k, p, T_data, 
              min_indiv_shocks = 0, 
              max_indiv_shocks = 0, 
              min_high_vol_periods = 0, 
              max_high_vol_periods = 0,
              min_exog_shocks = 0,
              max_exog_shocks = 0)
Yraw <- data[1:(T_data-h),]
Ypred <- data[(T_data-h+1):T_data,]

mh_params <- list(
  n_thin <- 1,
  scale_hess <- 0.06,
  adjust_burn <- 0.75,
  acc_lower <- 0.12,
  acc_upper <- 0.28,
  acc_change <- 0.05
)
hyper_params <- list(pi1_param = calc_gamma_pdf_params(m = 0.2, s = 0.4),
                     mu_param = calc_gamma_pdf_params(m = 1, s = 1),
                     gamma_param = calc_gamma_pdf_params(m = 1, s = 1),
                     s2_diag_param = list(a = 0.02^2, b = 0.02^2, 
                                  mode = hhelper_compute_sigma_vec(Yraw, p)))

start <- Sys.time()
res_mh <- run_bvar_hierarch(Yraw = Yraw, lag_mean = 1.0,
                            p = p_bvar, 
                            intercept = intercept,
                            hyper_params = hyper_params,
                            mh_params = mh_params)
Sys.time() - start
res_mh$acceptance_rate
plot(res_mh$hyper_parameters[,1], type = "l")


# ---- run SSVS ----------------------------------------------------------------
data <- sim_y(k, p, T_data, 
              min_indiv_shocks = 0, 
              max_indiv_shocks = 0, 
              min_high_vol_periods = 0, 
              max_high_vol_periods = 0,
              min_exog_shocks = 0,
              max_exog_shocks = 0)
Yraw <- data[1:(T_data-h),]
Ypred <- data[(T_data-h+1):T_data,]

start <- Sys.time()
res_ssvs <- run_bvar_ssvs(Yraw, p_bvar, 
                        intercept = intercept, 
                        use_dummies = FALSE,
                        dummy_pars = list(mu = 1.0, gamma = 1.0, prior_mean = NULL),
                        lag_mean = 1,
                        tau0 = 1/10,
                        tau1 = 10,
                        delta_prob = 0.8,
                        n_draws = 5000,
                        burnin = 1000)
Sys.time() - start
Y_forecast_ssvs_median <- bvar_forecast(Yraw, res_ssvs$Phi, res_ssvs$Sigma, p = p_bvar, h = 4, alpha = 0.05, intercept = intercept)$median
sqrt(mean((Y_forecast_ssvs_median - Ypred)^2))


# ---- inspect outputs ---------------------------------------------------------
# --- gibs
Phi_post_median   <- apply(res_gibbs$Phi, c(1,2), median)
Sigma_post_median <- apply(res_gibbs$Sigma, c(1,2), median)

# give names to results
varNames <- paste0("Y", 1:k)
rowNamesPhi <- c("Intercept")
for (i in 1:p_bvar) {
  rowNamesPhi <- c(rowNamesPhi, paste0("Lag", rep(i, k)))
}
if(!intercept) rowNamesPhi <- rowNamesPhi[-1]

dimnames(Phi_post_median) <- list(rowNamesPhi, varNames)
dimnames(Sigma_post_median) <- list(varNames, varNames)

# view results
Phi_post_median
Sigma_post_median


# --- mh
Phi_post_median   <- apply(res_mh$Phi, c(1,2), median)
Sigma_post_median <- apply(res_mh$Sigma, c(1,2), median)

# give names to results
varNames <- paste0("Y", 1:k)
rowNamesPhi <- c("Intercept")
for (i in 1:p_bvar) {
  rowNamesPhi <- c(rowNamesPhi, paste0("Lag", rep(i, k)))
}
if(!intercept) rowNamesPhi <- rowNamesPhi[-1]

dimnames(Phi_post_median) <- list(rowNamesPhi, varNames)
dimnames(Sigma_post_median) <- list(varNames, varNames)

# view results
Phi_post_median
Sigma_post_median


# --- ssvs
Phi_post_median   <- apply(res_ssvs$Phi, c(1,2), median)
Sigma_post_median <- apply(res_ssvs$Sigma, c(1,2), median)

# give names to results
varNames <- paste0("Y", 1:k)
rowNamesPhi <- c("Intercept")
for (i in 1:p_bvar) {
  rowNamesPhi <- c(rowNamesPhi, paste0("Lag", rep(i, k)))
}
if(!intercept) rowNamesPhi <- rowNamesPhi[-1]

dimnames(Phi_post_median) <- list(rowNamesPhi, varNames)
dimnames(Sigma_post_median) <- list(varNames, varNames)

# view results
Phi_post_median
Sigma_post_median
