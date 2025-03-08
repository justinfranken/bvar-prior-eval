# monte carlo simulating process


# ------------------------------------------------------------------------------
# import libraries
# ------------------------------------------------------------------------------

source(paste0(getwd(),"/lib.r"))


# ------------------------------------------------------------------------------
# import functions
# ------------------------------------------------------------------------------

function_files <- c(
  "_distr_samplers.R",
  "_monte_carlo_fun.R",
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


# ------------------------------------------------------------------------------
# Monte Carlo Simulation Setup
# ------------------------------------------------------------------------------

# ---- setting up simulation ---------------------------------------------------
# -- global simulation settings -----------------
n_iter <- 5000
n_draws <- 5000
burnin <- 1000
intercept <- TRUE
h <- 4
p_bvar <- 4
p <- 4
k <- 7
m <- k * p_bvar + as.integer(intercept)
lag_mean <- 1.0
alpha <- 0.05
dummy_pars <- list(mu = 1.0,
                   gamma = 1.0,
                   prior_mean = NULL)
mh_params <- list(
  n_thin <- 1,
  scale_hess <- 0.04,
  adjust_burn <- 0.75,
  acc_lower <- 0.12,
  acc_upper <- 0.28,
  acc_change <- 0.05
)


# -- data simulation settings -------------------
# - small data --------------
small_no_shocks_low_corr <- list(
  d_min = 0.30,
  d_max = 0.40,
  off_d_min = 0.02,
  off_d_max = 0.07,
  init_sd = 0.025,
  shock_diag_min = 0.02,
  shock_diag_max = 0.07,
  mean_vola = 0.035,
  sd_vola = 0.007,
  min_indiv_shocks = 0, 
  max_indiv_shocks = 0, 
  min_high_vol_periods = 0, 
  max_high_vol_periods = 0,
  min_exog_shocks = 0,
  max_exog_shocks = 0
)

small_with_shocks_low_corr <- list(
  d_min = 0.30,
  d_max = 0.40,
  off_d_min = 0.02,
  off_d_max = 0.07,
  init_sd = 0.025,
  shock_diag_min = 0.02,
  shock_diag_max = 0.07,
  mean_vola = 0.035,
  sd_vola = 0.007,
  min_indiv_shocks = 0, 
  max_indiv_shocks = 1, 
  min_high_vol_periods = 1, 
  max_high_vol_periods = 1,
  min_exog_shocks = 1,
  max_exog_shocks = 1
)

small_no_shocks_high_corr <- list(
  d_min = 0.30,
  d_max = 0.40,
  off_d_min = 0.05,
  off_d_max = 0.10,
  init_sd = 0.025,
  shock_diag_min = 0.40,
  shock_diag_max = 0.50,
  mean_vola = 0.035,
  sd_vola = 0.007,
  min_indiv_shocks = 0, 
  max_indiv_shocks = 0, 
  min_high_vol_periods = 0, 
  max_high_vol_periods = 0,
  min_exog_shocks = 0,
  max_exog_shocks = 0
)

small_with_shocks_high_corr <- list(
  d_min = 0.30,
  d_max = 0.40,
  off_d_min = 0.05,
  off_d_max = 0.10,
  init_sd = 0.025,
  shock_diag_min = 0.40,
  shock_diag_max = 0.50,
  mean_vola = 0.035,
  sd_vola = 0.007,
  min_indiv_shocks = 0, 
  max_indiv_shocks = 1, 
  min_high_vol_periods = 1, 
  max_high_vol_periods = 1,
  min_exog_shocks = 1,
  max_exog_shocks = 1
)


# - medium data --------------
medium_no_shocks_low_corr <- list(
  d_min = 0.30,
  d_max = 0.40,
  off_d_min = 0.02,
  off_d_max = 0.07,
  init_sd = 0.025,
  shock_diag_min = 0.02,
  shock_diag_max = 0.07,
  mean_vola = 0.035,
  sd_vola = 0.007,
  min_indiv_shocks = 0, 
  max_indiv_shocks = 0, 
  min_high_vol_periods = 0, 
  max_high_vol_periods = 0,
  min_exog_shocks = 0,
  max_exog_shocks = 0
)

medium_with_shocks_low_corr <- list(
  d_min = 0.30,
  d_max = 0.40,
  off_d_min = 0.02,
  off_d_max = 0.07,
  init_sd = 0.025,
  shock_diag_min = 0.02,
  shock_diag_max = 0.07,
  mean_vola = 0.035,
  sd_vola = 0.007,
  min_indiv_shocks = 1, 
  max_indiv_shocks = 1, 
  min_high_vol_periods = 2, 
  max_high_vol_periods = 2,
  min_exog_shocks = 2,
  max_exog_shocks = 2
)

medium_no_shocks_high_corr <- list(
  d_min = 0.30,
  d_max = 0.40,
  off_d_min = 0.05,
  off_d_max = 0.10,
  init_sd = 0.025,
  shock_diag_min = 0.40,
  shock_diag_max = 0.50,
  mean_vola = 0.035,
  sd_vola = 0.007,
  min_indiv_shocks = 0, 
  max_indiv_shocks = 0, 
  min_high_vol_periods = 0, 
  max_high_vol_periods = 0,
  min_exog_shocks = 0,
  max_exog_shocks = 0
)

medium_with_shocks_high_corr <- list(
  d_min = 0.30,
  d_max = 0.40,
  off_d_min = 0.05,
  off_d_max = 0.10,
  init_sd = 0.025,
  shock_diag_min = 0.40,
  shock_diag_max = 0.50,
  mean_vola = 0.035,
  sd_vola = 0.007,
  min_indiv_shocks = 1, 
  max_indiv_shocks = 1, 
  min_high_vol_periods = 2, 
  max_high_vol_periods = 2,
  min_exog_shocks = 2,
  max_exog_shocks = 2
)


# - large data --------------
large_no_shocks_low_corr <- list(
  d_min = 0.30,
  d_max = 0.40,
  off_d_min = 0.02,
  off_d_max = 0.07,
  init_sd = 0.025,
  shock_diag_min = 0.02,
  shock_diag_max = 0.07,
  mean_vola = 0.035,
  sd_vola = 0.007,
  min_indiv_shocks = 0, 
  max_indiv_shocks = 0, 
  min_high_vol_periods = 0, 
  max_high_vol_periods = 0,
  min_exog_shocks = 0,
  max_exog_shocks = 0
)

large_with_shocks_low_corr <- list(
  d_min = 0.30,
  d_max = 0.40,
  off_d_min = 0.02,
  off_d_max = 0.07,
  init_sd = 0.025,
  shock_diag_min = 0.02,
  shock_diag_max = 0.07,
  mean_vola = 0.035,
  sd_vola = 0.007,
  min_indiv_shocks = 1, 
  max_indiv_shocks = 2, 
  min_high_vol_periods = 3, 
  max_high_vol_periods = 3,
  min_exog_shocks = 3,
  max_exog_shocks = 3
)

large_no_shocks_high_corr <- list(
  d_min = 0.30,
  d_max = 0.40,
  off_d_min = 0.05,
  off_d_max = 0.10,
  init_sd = 0.025,
  shock_diag_min = 0.40,
  shock_diag_max = 0.50,
  mean_vola = 0.035,
  sd_vola = 0.007,
  min_indiv_shocks = 0, 
  max_indiv_shocks = 0, 
  min_high_vol_periods = 0, 
  max_high_vol_periods = 0,
  min_exog_shocks = 0,
  max_exog_shocks = 0
)

large_with_shocks_high_corr <- list(
  d_min = 0.30,
  d_max = 0.40,
  off_d_min = 0.05,
  off_d_max = 0.10,
  init_sd = 0.025,
  shock_diag_min = 0.40,
  shock_diag_max = 0.50,
  mean_vola = 0.035,
  sd_vola = 0.007,
  min_indiv_shocks = 1, 
  max_indiv_shocks = 2, 
  min_high_vol_periods = 3, 
  max_high_vol_periods = 3,
  min_exog_shocks = 3,
  max_exog_shocks = 3
)


# ------------------------------------------------------------------------------
# Monte Carlo Simulation
# ------------------------------------------------------------------------------

# ---- Monte Carlo Simulation small observations -------------------------------


n_obs <- 40

# DONE
small_out_no_shocks_low_corr <- simulation(n_iter = n_iter, 
                                           n_draws = n_draws, 
                                           burnin = burnin, 
                                           n_obs = n_obs, 
                                           k = k, p = p, m = m, h = h, 
                                           intercept = intercept, 
                                           p_bvar = p_bvar, 
                                           alpha = alpha, 
                                           lag_mean = lag_mean, 
                                           dummy_pars = dummy_pars, 
                                           mh_params = mh_params,  
                                           sim_params = small_no_shocks_low_corr)

# DONE
small_out_no_shocks_high_corr <- simulation(n_iter = n_iter, 
                                            n_draws = n_draws, 
                                            burnin = burnin, 
                                            n_obs = n_obs, 
                                            k = k, p = p, m = m, h = h, 
                                            intercept = intercept, 
                                            p_bvar = p_bvar, 
                                            alpha = alpha, 
                                            lag_mean = lag_mean, 
                                            dummy_pars = dummy_pars, 
                                            mh_params = mh_params,  
                                            sim_params = small_no_shocks_high_corr)

# DONE
small_out_with_shocks_low_corr <- simulation(n_iter = n_iter, 
                                             n_draws = n_draws, 
                                             burnin = burnin, 
                                             n_obs = n_obs, 
                                             k = k, p = p, m = m, h = h, 
                                             intercept = intercept, 
                                             p_bvar = p_bvar, 
                                             alpha = alpha, 
                                             lag_mean = lag_mean, 
                                             dummy_pars = dummy_pars, 
                                             mh_params = mh_params,  
                                             sim_params = small_with_shocks_low_corr)

# DONE
small_out_with_shocks_high_corr <- simulation(n_iter = n_iter, 
                                              n_draws = n_draws, 
                                              burnin = burnin, 
                                              n_obs = n_obs, 
                                              k = k, p = p, m = m, h = h, 
                                              intercept = intercept, 
                                              p_bvar = p_bvar, 
                                              alpha = alpha, 
                                              lag_mean = lag_mean, 
                                              dummy_pars = dummy_pars, 
                                              mh_params = mh_params,  
                                              sim_params = small_with_shocks_high_corr)


saveRDS(small_out_no_shocks_low_corr, file="small_out_no_shocks_low_corr.RData") # DONE
saveRDS(small_out_no_shocks_high_corr, file="small_out_no_shocks_high_corr.RData") # DONE
saveRDS(small_out_with_shocks_low_corr, file="small_out_with_shocks_low_corr.RData") # DONE
saveRDS(small_out_with_shocks_high_corr, file="small_out_with_shocks_high_corr.RData") # DONE

# ---- Monte Carlo Simulation medium observations ------------------------------

n_obs <- 80

# DONE
medium_out_no_shocks_low_corr <- simulation(n_iter = n_iter, 
                                            n_draws = n_draws, 
                                            burnin = burnin, 
                                            n_obs = n_obs, 
                                            k = k, p = p, m = m, h = h, 
                                            intercept = intercept, 
                                            p_bvar = p_bvar, 
                                            alpha = alpha, 
                                            lag_mean = lag_mean, 
                                            dummy_pars = dummy_pars, 
                                            mh_params = mh_params,  
                                            sim_params = medium_no_shocks_low_corr)

# DONE
medium_out_no_shocks_high_corr <- simulation(n_iter = n_iter, 
                                             n_draws = n_draws, 
                                             burnin = burnin, 
                                             n_obs = n_obs, 
                                             k = k, p = p, m = m, h = h, 
                                             intercept = intercept, 
                                             p_bvar = p_bvar, 
                                             alpha = alpha, 
                                             lag_mean = lag_mean, 
                                             dummy_pars = dummy_pars, 
                                             mh_params = mh_params,  
                                             sim_params = medium_no_shocks_high_corr)

# DONE
medium_out_with_shocks_low_corr <- simulation(n_iter = n_iter, 
                                              n_draws = n_draws, 
                                              burnin = burnin, 
                                              n_obs = n_obs, 
                                              k = k, p = p, m = m, h = h, 
                                              intercept = intercept, 
                                              p_bvar = p_bvar, 
                                              alpha = alpha, 
                                              lag_mean = lag_mean, 
                                              dummy_pars = dummy_pars, 
                                              mh_params = mh_params,  
                                              sim_params = medium_with_shocks_low_corr)

# DONE
medium_out_with_shocks_high_corr <- simulation(n_iter = n_iter, 
                                               n_draws = n_draws, 
                                               burnin = burnin, 
                                               n_obs = n_obs, 
                                               k = k, p = p, m = m, h = h, 
                                               intercept = intercept, 
                                               p_bvar = p_bvar, 
                                               alpha = alpha, 
                                               lag_mean = lag_mean, 
                                               dummy_pars = dummy_pars, 
                                               mh_params = mh_params,  
                                               sim_params = medium_with_shocks_high_corr)


saveRDS(medium_out_no_shocks_low_corr, file="medium_out_no_shocks_low_corr.RData") # DONE
saveRDS(medium_out_no_shocks_high_corr, file="medium_out_no_shocks_high_corr.RData") # DONE
saveRDS(medium_out_with_shocks_low_corr, file="medium_out_with_shocks_low_corr.RData") # DONE
saveRDS(medium_out_with_shocks_high_corr, file="medium_out_with_shocks_high_corr.RData") # DONE

# ---- Monte Carlo Simulation large observations -------------------------------

n_obs <- 200

# DONE
large_out_no_shocks_low_corr <- simulation(n_iter = n_iter, 
                                           n_draws = n_draws, 
                                           burnin = burnin, 
                                           n_obs = n_obs, 
                                           k = k, p = p, m = m, h = h, 
                                           intercept = intercept, 
                                           p_bvar = p_bvar, 
                                           alpha = alpha, 
                                           lag_mean = lag_mean, 
                                           dummy_pars = dummy_pars, 
                                           mh_params = mh_params,  
                                           sim_params = large_no_shocks_low_corr)
# DONE
large_out_no_shocks_high_corr <- simulation(n_iter = n_iter, 
                                            n_draws = n_draws, 
                                            burnin = burnin, 
                                            n_obs = n_obs, 
                                            k = k, p = p, m = m, h = h, 
                                            intercept = intercept, 
                                            p_bvar = p_bvar, 
                                            alpha = alpha, 
                                            lag_mean = lag_mean, 
                                            dummy_pars = dummy_pars, 
                                            mh_params = mh_params,  
                                            sim_params = large_no_shocks_high_corr)

# DONE
large_out_with_shocks_low_corr <- simulation(n_iter = n_iter, 
                                             n_draws = n_draws, 
                                             burnin = burnin, 
                                             n_obs = n_obs, 
                                             k = k, p = p, m = m, h = h, 
                                             intercept = intercept, 
                                             p_bvar = p_bvar, 
                                             alpha = alpha, 
                                             lag_mean = lag_mean, 
                                             dummy_pars = dummy_pars, 
                                             mh_params = mh_params,  
                                             sim_params = large_with_shocks_low_corr)

# DONE
large_out_with_shocks_high_corr <- simulation(n_iter = n_iter, 
                                              n_draws = n_draws, 
                                              burnin = burnin, 
                                              n_obs = n_obs, 
                                              k = k, p = p, m = m, h = h, 
                                              intercept = intercept, 
                                              p_bvar = p_bvar, 
                                              alpha = alpha, 
                                              lag_mean = lag_mean, 
                                              dummy_pars = dummy_pars, 
                                              mh_params = mh_params,  
                                              sim_params = large_with_shocks_high_corr)

saveRDS(large_out_no_shocks_low_corr, file="large_out_no_shocks_low_corr.RData") # DONE
saveRDS(large_out_no_shocks_high_corr, file="large_out_no_shocks_high_corr.RData") # DONE
saveRDS(large_out_with_shocks_low_corr, file="large_out_with_shocks_low_corr.RData") # DONE
saveRDS(large_out_with_shocks_high_corr, file="large_out_with_shocks_high_corr.RData") # DONE

# ------------------------------------------------------------------------------
# evaluate simulation results
# ------------------------------------------------------------------------------

# --- clarify models being used ------------------------------------------------
models <- c("classic_mn", "hierarch_mn", "ssvs", "flat_mn", "var")

# -- small ----------------------------
small_no_shock_low_corr <- evaluate_sim_res(small_out_no_shocks_low_corr, 
                                             models = models, 
                                             h = h, 
                                             n_iter = n_iter) # DONE
small_with_shock_low_corr <- evaluate_sim_res(small_out_with_shocks_low_corr, 
                                               models = models, 
                                               h = h, 
                                               n_iter = n_iter) # DONE
small_no_shock_high_corr <- evaluate_sim_res(small_out_no_shocks_high_corr, 
                                             models = models, 
                                             h = h, 
                                             n_iter = n_iter) # DONE
small_with_shock_high_corr <- evaluate_sim_res(small_out_with_shocks_high_corr, 
                                               models = models, 
                                               h = h, 
                                               n_iter = n_iter) # DONE

# -- medium ---------------------------
medium_no_shock_low_corr <- evaluate_sim_res(medium_out_no_shocks_low_corr, 
                                             models = models, 
                                             h = h, 
                                             n_iter = n_iter) # DONE
medium_with_shock_low_corr <- evaluate_sim_res(medium_out_with_shocks_low_corr, 
                                               models = models, 
                                               h = h, 
                                               n_iter = n_iter) # DONE
medium_no_shock_high_corr <- evaluate_sim_res(medium_out_no_shocks_high_corr, 
                                             models = models, 
                                             h = h, 
                                             n_iter = n_iter) # DONE
medium_with_shock_high_corr <- evaluate_sim_res(medium_out_with_shocks_high_corr, 
                                               models = models, 
                                               h = h, 
                                               n_iter = n_iter) # DONE

# -- large ----------------------------
large_no_shock_low_corr <- evaluate_sim_res(large_out_no_shocks_low_corr, 
                                             models = models, 
                                             h = h, 
                                             n_iter = n_iter) # DONE
large_with_shock_low_corr <- evaluate_sim_res(large_out_with_shocks_low_corr, 
                                               models = models, 
                                               h = h, 
                                               n_iter = n_iter)
large_no_shock_high_corr <- evaluate_sim_res(large_out_no_shocks_high_corr, 
                                             models = models, 
                                             h = h, 
                                             n_iter = n_iter) # DONE
large_with_shock_high_corr <- evaluate_sim_res(large_out_with_shocks_high_corr, 
                                               models = models, 
                                               h = h, 
                                               n_iter = n_iter)
