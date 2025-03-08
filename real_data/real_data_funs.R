# functions to evaluate real data


run_one_evaluation <- function(train_data, Ypred,
                               n_draws, burnin, k, p, m, h, 
                               intercept, alpha, 
                               lag_mean, dummy_pars, mh_params){
  #' This function evaluates and fits several BVAR models using different prior 
  #' specifications and estimation methods (including classic Minnesota with inverse Wishart, 
  #' hierarchical Minnesota, SSVS, Minnesota with flat priors, and a frequentist VAR), generates 
  #' forecasts from each model, and computes forecast evaluation metrics (RMSE and prediction 
  #' interval accuracy). 
  #'
  #' Parameters:
  #' - train_data (matrix): The real data matrix the BVAR should learn from.
  #' - Ypred (matrix): The next h following real realizations of the data matrix.
  #' - n_draws (integer): Total number of posterior draws to be used in the Monte Carlo simulation.
  #' - burnin (integer): Number of initial posterior draws to discard as burn-in.
  #' - k (integer): Number of endogenous variables in the VAR model.
  #' - p (integer): Lag order to be used for the VAR estimation (may be equal to p).
  #' - m (integer): The total number of parameters in the VAR model.
  #' - h (integer): Forecast horizon (number of steps ahead to forecast).
  #' - intercept (logical): Indicator whether to include an intercept in the VAR model.
  #' - alpha (numeric): Significance level for forecast prediction intervals (e.g., 0.05 for 95% intervals).
  #' - lag_mean (numeric): Lag mean parameter used in the Minnesota prior.
  #' - dummy_pars (list): Parameters for generating dummy observations, including mu, gamma, and 
  #'   (optionally) a prior mean.
  #' - mh_params (list): Parameters controlling the Metropolis-Hastings sampler in the hierarchical 
  #'   Minnesota approach (e.g., thinning interval, scaling factors, burn-in adjustment, acceptance rate thresholds).
  #'
  #' Returns:
  #' - A list containing forecast evaluation results for different methods:
  #'    - gibbs: Results for classic Minnesota with inverse Wishart.
  #'    - mh: Results for the hierarchical Minnesota approach.
  #'    - ssvs: Results for the SSVS (stochastic search variable selection) approach.
  #'    - flat: Results for the classic Minnesota approach with flat priors.
  #'    - var: Results for the frequentist VAR.
  
  # ---- rigorous methods ------------------------------------------------------
  # -- classic minnesota with inverse wishart -------------
  # posterior coefficient draws
  res_gibbs <- run_bvar_minnesota(
    Yraw = train_data,
    p = p,
    intercept = intercept,
    use_dummies = FALSE,
    dummy_pars = dummy_pars,
    use_flat = FALSE,
    lag_mean = lag_mean,
    pi1 = 0.2,
    pi3 = 1,
    pi4 = 1000,
    sigma_vec = NULL, 
    n_draws = n_draws,
    burnin = burnin
  )
  
  # forecasts and evaluation
  fcst_gibbs <- bvar_forecast(Y = train_data, 
                              phi_draws = res_gibbs$Phi, 
                              Sigma_draws = res_gibbs$Sigma, 
                              p = p, 
                              h = h, 
                              alpha = alpha)
  rmse_gibbs <- fcst_rmse(fcst_obj = fcst_gibbs, Ypred = Ypred, h = h, k = k)
  pred_acc_gibbs <- fcst_pred_int_acc(fcst_obj = fcst_gibbs, Ypred = Ypred, h = h)
  
  
  # -- hierarchical minnesota -----------------------------
  # get hyper parameters for mh
  hyper_params <- list(pi1_param = calc_gamma_pdf_params(m = 0.2, s = 0.4),
                       mu_param = calc_gamma_pdf_params(m = 1, s = 1),
                       gamma_param = calc_gamma_pdf_params(m = 1, s = 1),
                       s2_diag_param = list(
                         a = 0.02^2, 
                         b = 0.02^2, 
                         mode = hhelper_compute_sigma_vec(train_data, p)
                       )
  )
  
  # posterior coefficient draws
  res_mh <- run_bvar_hierarch(
    Yraw = train_data, 
    p = p, 
    intercept = intercept,
    mh_params = mh_params,
    hyper_params = hyper_params,
    pi4 = 1000,
    lag_mean = lag_mean, 
    n_draws = n_draws, 
    burnin = burnin
  )
  
  # forecasts and evaluation
  fcst_mh <- bvar_forecast(Y = train_data, 
                           phi_draws = res_mh$Phi, 
                           Sigma_draws = res_mh$Sigma, 
                           p = p, 
                           h = h, 
                           alpha = alpha)
  rmse_mh <- fcst_rmse(fcst_obj = fcst_mh, Ypred = Ypred, h = h, k = k)
  pred_acc_mh <- fcst_pred_int_acc(fcst_obj = fcst_mh, Ypred = Ypred, h = h)
  
  
  # -- ssvs -----------------------------------------------
  # posterior coefficient draws
  res_ssvs <- run_bvar_ssvs(
    train_data, p, 
    intercept = intercept, 
    use_dummies = FALSE,
    dummy_pars = dummy_pars,
    lag_mean = lag_mean,
    tau0 = 1/100,
    tau1 = 1,
    delta_prob = 0.8,
    n_draws = n_draws,
    burnin = burnin
  )
  
  # forecasts and evaluation
  fcst_ssvs <- bvar_forecast(Y = train_data, 
                             phi_draws = res_ssvs$Phi, 
                             Sigma_draws = res_ssvs$Sigma, 
                             p = p, 
                             h = h, 
                             alpha = alpha)
  rmse_ssvs <- fcst_rmse(fcst_obj = fcst_ssvs, Ypred = Ypred, h = h, k = k)
  pred_acc_ssvs <- fcst_pred_int_acc(fcst_obj = fcst_ssvs, Ypred = Ypred, h = h)
  
  # ---- baseline methods ------------------------------------------------------
  # -- classic minnesota with flat priors -----------------
  res_flat <- run_bvar_minnesota(
    Yraw = train_data,
    p = p,
    intercept = intercept,
    use_dummies = FALSE,
    dummy_pars = dummy_pars,
    use_flat = TRUE,
    lag_mean = lag_mean,
    pi1 = 0.2,
    pi3 = 1,
    pi4 = 1000,
    sigma_vec = NULL, 
    n_draws = n_draws,
    burnin = burnin
  )
  
  # forecasts and evaluation
  fcst_flat <- bvar_forecast(Y = train_data, 
                             phi_draws = res_flat$Phi, 
                             Sigma_draws = res_flat$Sigma, 
                             p = p, 
                             h = h, 
                             alpha = alpha)
  rmse_flat <- fcst_rmse(fcst_obj = fcst_flat, Ypred = Ypred, h = h, k = k)
  pred_acc_flat <- fcst_pred_int_acc(fcst_obj = fcst_flat, Ypred = Ypred, h = h)
  
  
  # -- frequentist VAR ------------------------------------
  # compute coefficients
  colnames(train_data) <- paste0("Y", 1:k)
  res_var <- VAR(train_data, p = p, type = "const")
  
  # forecasts and evaluation
  temp_var <- predict(res_var, n.ahead = h, ci = 1 - alpha)
  fcst_var <- list(
    median = sapply(temp_var$fcst, function(x) x[,1]),
    lower  = sapply(temp_var$fcst, function(x) x[,2]),
    upper  = sapply(temp_var$fcst, function(x) x[,3])
  )
  rmse_var <- fcst_rmse(fcst_obj = fcst_var, Ypred = Ypred, h = h, k = k)
  pred_acc_var <- fcst_pred_int_acc(fcst_obj = fcst_var, Ypred = Ypred, h = h)
  
  
  # ---- return forecast evaluations -------------------------------------------
  return(list(
    gibbs = list(rmse = rmse_gibbs, pred_acc = pred_acc_gibbs),
    mh = list(rmse = rmse_mh, pred_acc = pred_acc_mh),
    ssvs = list(rmse = rmse_ssvs, pred_acc = pred_acc_ssvs),
    flat = list(rmse = rmse_flat, pred_acc = pred_acc_flat),
    var = list(rmse = rmse_var, pred_acc = pred_acc_var)
  ))
}


run_real_data_evaluation <- function(data_mat, n_obs, h, k, p, dummy_pars, mh_params, 
                                     intercept = TRUE, alpha = 0.05, 
                                     n_draws = 5000, burnin = 1000){
  #' Evaluates forecasting performance on real data using a rolling window approach.
  #'
  #' Parameters:
  #'   - data_mat: The data matrix containing the real data.
  #'   - n_obs: The number of observations to be used for training in each rolling window.
  #'   - h: The forecast horizon.
  #'   - n_draws: Number of posterior draws for each evaluation.
  #'   - burnin: Number of initial draws to discard as burn-in.
  #'   - k, p: Model dimensions (k = number of variables, p = lag order).
  #'   - alpha: Significance level for forecast prediction intervals.
  #'   - dummy_pars: A list of parameters for dummy observations.
  #'   - mh_params: A list of parameters for the Metropolis-Hastings sampler.
  #'
  #' Returns:
  #' - A list containing:
  #'     - results: A list of evaluation outputs for each rolling window.
  
  total <- length(data_mat[,1]) - h - n_obs + 1
  m <- k * p + as.integer(intercept)
  
  result <- list()
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  for (iter in seq(total)) {
    # get rolling window 
    chem_train <- data_mat[iter:(n_obs+iter-1),]
    chem_true <- data_mat[(n_obs+iter):(n_obs+h+iter-1),]
    
    # run evaluation
    out <- run_one_evaluation(chem_train, 
                              chem_true, 
                              n_draws, 
                              burnin, 
                              k, p, m, h, 
                              intercept = intercept, 
                              alpha = alpha, 
                              lag_mean = 1, 
                              dummy_pars = dummy_pars, 
                              mh_params = mh_params)
    
    # store results
    result[[iter]] <- out
    setTxtProgressBar(pb, iter)
  }
  res <- list(results = result)
  return(res)
}
