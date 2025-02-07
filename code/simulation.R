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
n_iter <- 10
n_draws <- 5000
burnin <- 500
intercept <- TRUE
h <- 4
p_bvar <- 4
p <- 3
k <- 7
m <- k * p + as.integer(intercept)
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


# -- sim_pars settings --------------------------
sim_params_no_shocks <- list(
  min_indiv_shocks = 0, 
  max_indiv_shocks = 0, 
  min_high_vol_periods = 0, 
  max_high_vol_periods = 0,
  min_exog_shocks = 0,
  max_exog_shocks = 0
)

sim_params_with_shocks <- list(
  min_indiv_shocks = 0, 
  max_indiv_shocks = 1, 
  min_high_vol_periods = 1, 
  max_high_vol_periods = 1,
  min_exog_shocks = 1,
  max_exog_shocks = 1
)


# ---- Monte Carlo Simulation small observations -------------------------------

n_obs <- 40

start <- Sys.time()
small_out_no_shock <- simulation(n_iter = n_iter, 
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
                                 sim_params = sim_params_no_shocks)
end <- Sys.time() - start

small_out_with_shock <- simulation(n_iter = n_iter, 
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
                                 sim_params = sim_params_with_shocks)


# ---- Monte Carlo Simulation medium observations ------------------------------

n_obs <- 60

medium_out_no_shock <- simulation(n_iter = n_iter, 
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
                                 sim_params = sim_params_no_shocks)


medium_out_with_shock <- simulation(n_iter = n_iter, 
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
                                   sim_params = sim_params_with_shocks)


# ---- Monte Carlo Simulation large observations -------------------------------

n_obs <- 80

large_out_no_shock <- simulation(n_iter = n_iter, 
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
                                  sim_params = sim_params_no_shocks)


large_out_with_shock <- simulation(n_iter = n_iter, 
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
                                    sim_params = sim_params_with_shocks)


# ------------------------------------------------------------------------------
# evaluate simulation results
# ------------------------------------------------------------------------------

# --- clarify models being used ------------------------------------------------
models <- c("classic_mn", "hierarch_mn", "ssvs", "flat_mn", "var")

small_no_shock <- evaluate_sim_res(small_out_no_shock, models = models, h = h, n_iter = n_iter)
small_with_shock <- evaluate_sim_res(small_out_with_shock, models = models, h = h, n_iter = n_iter)

medium_no_shock <- evaluate_sim_res(medium_out_no_shock, models = models, h = h, n_iter = n_iter)
medium_with_shock <- evaluate_sim_res(medium_out_with_shock, models = models, h = h, n_iter = n_iter)

large_no_shock <- evaluate_sim_res(large_out_no_shock, models = models, h = h, n_iter = n_iter)
large_with_shock <- evaluate_sim_res(large_out_with_shock, models = models, h = h, n_iter = n_iter)

