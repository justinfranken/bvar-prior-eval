# evaluation of mcmc convergence


# ------------------------------------------------------------------------------
# import libraries
# ------------------------------------------------------------------------------

source(paste0(getwd(),"/lib.r"))


# ------------------------------------------------------------------------------
# import functions
# ------------------------------------------------------------------------------

function_files <- c(
  "_distr_samplers.R",
  "_hierarch_fun.R",
  "_data_prep.R",
  "_mcmc.R",
  "_ssvs.R",
  "_prior_specification.R",
  "bvar.R",
  "simulate_data.R",
  "_mcmc_convergence.R"
)
for (i in 1:length(function_files)) {
  source(paste0(getwd(),"/functions/", function_files[i]))
}
rm(i)
rm(function_files)


# ------------------------------------------------------------------------------
# simulate data
# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------
# get posterior chain draws
# ------------------------------------------------------------------------------

# ---- set hyper parameters and setting ----------------------------------------
intercept = TRUE
p_bvar <- 4
dummy_pars_ex <- list(mu = 1.0,
                      gamma = 1.0,
                      prior_mean = NULL)
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


# ---- run Metropolis Hastings sampler -----------------------------------------
res_mh <- run_bvar_hierarch(Yraw = Yraw, 
                            pi4 = 1000, 
                            lag_mean = 1.0,
                            p = p_bvar, 
                            intercept = intercept,
                            hyper_params = hyper_params,
                            mh_params = mh_params,
                            n_draws = 5000, 
                            burnin = 1000)


# ---- run SSVS sampler --------------------------------------------------------
res_ssvs <- run_bvar_ssvs(Yraw, 
                          p_bvar, 
                          intercept = intercept, 
                          use_dummies = FALSE,
                          dummy_pars = dummy_pars_ex,
                          lag_mean = 1,
                          tau0 = 1/100,
                          tau1 = 1,
                          delta_prob = 0.8,
                          n_draws = 5000,
                          burnin = 1000)


# ------------------------------------------------------------------------------
# run convergence tests
# ------------------------------------------------------------------------------

# ---- run Metropolis Hastings sampler -----------------------------------------
geweke_coef(res_mh, 1, intercept = intercept, alpha = 0.05)
geweke_hyp(res_mh, alpha = 0.05)
brooks_plot(res_mh, B = 500, coef = 1, coef_row = 3, coef_col = 2)
traceplot_coef(res_mh, coef = 1, coef_row = 3, coef_col = 2)
traceplot_hyper(res_mh, hyp_para = 1)


# ---- run SSVS sampler --------------------------------------------------------
geweke_coef(res_ssvs, 1, intercept = intercept, alpha = 0.05)
brooks_plot(res_ssvs, B = 500, coef = 1, coef_row = 3, coef_col = 2)
traceplot_coef(res_ssvs, coef = 1, coef_row = 3, coef_col = 2)
traceplot_hyper(res_ssvs, ssvs = TRUE)

