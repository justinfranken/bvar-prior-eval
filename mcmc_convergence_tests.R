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
p <- 4
h <- 4
n_obs <- 40
T_data <- n_obs + h
data <- sim_y(k, p, T_data, 
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
              max_exog_shocks = 1)
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
# test different convergence tests
# ------------------------------------------------------------------------------

# ---- Metropolis Hastings sampler -----------------------------------------
geweke_coef(res_mh, 1, intercept = intercept, alpha = 0.05)
geweke_hyp(res_mh, alpha = 0.05)
brooks_plot(res_mh, B = 1000, coef = 1, coef_row = 3, coef_col = 2)
traceplot_coef(res_mh, coef = 1, coef_row = 3, coef_col = 2)
traceplot_hyper(res_mh, hyp_para = 1)


# ---- SSVS sampler --------------------------------------------------------
geweke_coef(res_ssvs, 1, intercept = intercept, alpha = 0.05)
brooks_plot(res_ssvs, B = 1000, coef = 1, coef_row = 3, coef_col = 2)
traceplot_coef(res_ssvs, coef = 1, coef_row = 3, coef_col = 2)
traceplot_hyper(res_ssvs, ssvs = TRUE)

# ---- plot for thesis ----------------------------------------------------
par(mar = c(2, 2, 0.1, 0.1),
    lwd = 2)
brooks_plot(res_mh, 
            B = 1000, 
            coef = 1, 
            coef_row = 2, 
            coef_col = 1, 
            main = "",
            cex_axis = 1.5)
geweke_coef(res_mh, 1, intercept = intercept, alpha = 0.05)
dev.off()
# ------------------------------------------------------------------------------
# create specific traceplots
# ------------------------------------------------------------------------------

# ---- SSVS sampler --------------------------------------------------------
par(mfrow = c(3,1), mar = c(2.5, 2.5, 2.5, 0.5))
plot(rowMeans(res_ssvs$delta_draws), 
     type = "l", 
     ylim = c(0,1),
     main = "Small Sample Size",
     cex.axis = 1.5,
     cex.main = 1.5)
abline(h = 7 / 29, col = "red", lty = 2)
abline(v = 1000, col = "red")

plot(rowMeans(res_ssvs$delta_draws), 
     type = "l", 
     ylim = c(0,1),
     main = "Medium Sample Size",
     cex.axis = 1.5,
     cex.main = 1.3)
abline(h = 7 / 29, col = "red", lty = 2)
abline(v = 1000, col = "red")

plot(rowMeans(res_ssvs$delta_draws), 
     type = "l", 
     ylim = c(0,1),
     main = "Large Sample Size",
     cex.axis = 1.5,
     cex.main = 1.3)
abline(h = 7 / 29, col = "red", lty = 2)
abline(v = 1000, col = "red")
dev.off()


# ---- Metropolis Hastings sampler -----------------------------------------
par(mfrow = c(4,1), mar = c(2.5, 2.5, 2.5, 0.5))
plot(res_mh$hyper_parameters[,1], 
     type = "l",
     main = "π1 Hyperparameter Traceplot for Large Sample Size")
abline(h = mean(res_mh$hyper_parameters[,1]), col = "red", lty = 2)
plot(res_mh$hyper_parameters[,2], 
     type = "l",
     main = "λ Hyperparameter Traceplot for Large Sample Size")
abline(h = mean(res_mh$hyper_parameters[,2]), col = "red", lty = 2)
plot(res_mh$hyper_parameters[,3], 
     type = "l",
     main = "ϑ Hyperparameter Traceplot for Large Sample Size")
abline(h = mean(res_mh$hyper_parameters[,3]), col = "red", lty = 2)
plot(res_mh$hyper_parameters[,4], 
     type = "l",
     main = "First Hyperparameter κ1 Traceplot for Large Sample Size")
abline(h = mean(res_mh$hyper_parameters[,4]), col = "red", lty = 2)

