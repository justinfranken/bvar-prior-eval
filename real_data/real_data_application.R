# evaluating real financial data


# ------------------------------------------------------------------------------
# prepare real data 
# ------------------------------------------------------------------------------

# ---- import and transform to matrix ------------------------------------------
load("data/chemicals.Rda")
load("data/building.Rda")
load("data/retail.Rda")
load("data/software.Rda")
load("data/banks.Rda")
load("data/oag.Rda")
chem_mat <- as.matrix(chemicals[-59,-1])
build_mat <- as.matrix(building[-59,-1])
retail_mat <- as.matrix(retail[-59,-1])
soft_mat <- as.matrix(software[-59,-1])
bank_mat <- as.matrix(banks[-59,-1])
oag_mat <- as.matrix(oag[,-1])
rm(list = c("chemicals", "building", "retail", "software", "banks", "oag"))


# ------------------------------------------------------------------------------
# import BVAR functions
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
  "bvar.R"
)
for (i in 1:length(function_files)) {
  source(file.path(dirname(getwd()), "functions", function_files[i]))
}
source(paste0(getwd(),"/real_data_funs.r"))
rm(i)
rm(function_files)


# ------------------------------------------------------------------------------
# running rolling window evaluations
# ------------------------------------------------------------------------------

# ---- setting up evaluations --------------------------------------------------
n_obs <- 40
h <- 4
k <- 7
p <- 4
m <- k * p + 1
alpha <- 0.05
n_draws <- 5000
burnin <- 1000
total <- length(chem_mat[,1]) - h - n_obs + 1


mh_params <- list(
  n_thin <- 1,
  scale_hess <- 0.04,
  adjust_burn <- 0.75,
  acc_lower <- 0.12,
  acc_upper <- 0.28,
  acc_change <- 0.05
)

dummy_pars <- list(mu = 1.0,
                   gamma = 1.0,
                   prior_mean = NULL)


# ---- run rolling windows -----------------------------------------------------
res_chem <- run_real_data_evaluation(data_mat = chem_mat, 
                                     n_obs, 
                                     h, k, p, 
                                     dummy_pars, 
                                     mh_params, 
                                     intercept = TRUE, 
                                     alpha = 0.05, 
                                     n_draws = 5000, 
                                     burnin = 1000)

res_build <- run_real_data_evaluation(data_mat = build_mat, 
                                      n_obs, 
                                      h, k, p, 
                                      dummy_pars, 
                                      mh_params, 
                                      intercept = TRUE, 
                                      alpha = 0.05, 
                                      n_draws = 5000, 
                                      burnin = 1000)

res_retail <- run_real_data_evaluation(data_mat = retail_mat, 
                                       n_obs, 
                                       h, k, p, 
                                       dummy_pars, 
                                       mh_params, 
                                       intercept = TRUE, 
                                       alpha = 0.05, 
                                       n_draws = 5000, 
                                       burnin = 1000)

res_soft <- run_real_data_evaluation(data_mat = soft_mat, 
                                     n_obs, 
                                     h, k, p, 
                                     dummy_pars, 
                                     mh_params, 
                                     intercept = TRUE, 
                                     alpha = 0.05, 
                                     n_draws = 5000, 
                                     burnin = 1000)

res_bank <- run_real_data_evaluation(data_mat = bank_mat, 
                                     n_obs, 
                                     h, k, p, 
                                     dummy_pars, 
                                     mh_params, 
                                     intercept = TRUE, 
                                     alpha = 0.05, 
                                     n_draws = 5000, 
                                     burnin = 1000)

res_oag <- run_real_data_evaluation(data_mat = oag_mat, 
                                    n_obs, 
                                    h, k, p, 
                                    dummy_pars, 
                                    mh_params, 
                                    intercept = TRUE, 
                                    alpha = 0.05, 
                                    n_draws = 5000, 
                                    burnin = 1000)

# ---- evaluate results --------------------------------------------------------
chem <- evaluate_sim_res(res_chem, c("niw", "hierarch", "ssvs", "flat", "var"), h, total)
build <- evaluate_sim_res(res_build, c("niw", "hierarch", "ssvs", "flat", "var"), h, total)
retail <- evaluate_sim_res(res_retail, c("niw", "hierarch", "ssvs", "flat", "var"), h, total)
soft <- evaluate_sim_res(res_soft, c("niw", "hierarch", "ssvs", "flat", "var"), h, total)
bank <- evaluate_sim_res(res_bank, c("niw", "hierarch", "ssvs", "flat", "var"), h, total)
oag <- evaluate_sim_res(res_oag, c("niw", "hierarch", "ssvs", "flat", "var"), h, total)
