# function which do monte carlo simulations and evaluate its results.


simulation <- function(n_iter,
                       n_draws, burnin, n_obs, k, p, m, h,
                       intercept, p_bvar, alpha, lag_mean,
                       dummy_pars, mh_params, sim_params,
                       max_tries = 5) {
  #' Runs multiple Monte Carlo simulation iterations in parallel.
  #'
  #' Parameters:
  #' - n_iter (integer): The total number of simulation iterations to run.
  #' - n_draws (integer): Number of posterior draws to use in each simulation iteration.
  #' - burnin (integer): Number of initial draws to discard as burn-in.
  #' - n_obs (integer): Number of observations for estimation (excluding the forecast horizon).
  #' - k (integer): Number of endogenous variables in the VAR model.
  #' - p (integer): Lag order used for data simulation.
  #' - m (integer): Total number of parameters in the VAR model.
  #' - h (integer): Forecast horizon (number of steps ahead to forecast).
  #' - intercept (logical): Indicates whether an intercept is included in the VAR model.
  #' - p_bvar (integer): Lag order used in the VAR estimation.
  #' - alpha (numeric): Significance level for forecast prediction intervals (e.g., 0.05 for 95% intervals).
  #' - lag_mean (numeric): Lag mean parameter for the Minnesota prior.
  #' - dummy_pars (list): Parameters for dummy observations (e.g., mu, gamma, and optionally prior_mean).
  #' - mh_params (list): Parameters controlling the Metropolis-Hastings sampler in the hierarchical VAR.
  #' - sim_params (list): Parameters for the data simulation (e.g., shock intensities, volatility periods).
  #' - max_tries (integer, default = 5): Maximum number of attempts allowed for each iteration if errors occur.
  #'
  #' Returns:
  #' - A list containing:
  #'     - results: A list of simulation results from all successful iterations.
  #'     - total_error_count: The total number of errors encountered (i.e., the sum of errors across iterations).
  #'     - n_errors: The number of iterations that ended in an error after max_tries.
  #'     - error_messages: A list of error messages from the iterations that failed.
  
  # set up a parallel plan
  plan(multisession)
  
  # handler
  handlers(handler_progress(format = ":current/:total [:bar] :percent  [Elapsed time: :elapsed | ETA: :eta]"))
  
  # run 'n_iter' times in parallel
  results_list <- with_progress({
    prog <- progressor(steps = n_iter + 1)
    future_lapply(
      X = seq_len(n_iter),
      FUN = function(i) {
        prog()
        run_one_iteration(
          n_draws    = n_draws,
          burnin     = burnin,
          n_obs      = n_obs,
          k          = k,
          p          = p,
          m          = m,
          h          = h,
          intercept  = intercept,
          p_bvar     = p_bvar,
          alpha      = alpha,
          lag_mean   = lag_mean,
          dummy_pars = dummy_pars,
          mh_params  = mh_params,
          sim_params = sim_params,
          max_tries  = max_tries
        )
      },
      future.seed = TRUE
    )
  })
  
  # separate successes from errors
  successful_runs <- Filter(function(x) x$status == "success", results_list)
  error_runs      <- Filter(function(x) x$status == "error",   results_list)
  
  # summarize total error count from all runs
  total_error_count <- sum(
    vapply(successful_runs, function(x) x$error_count, numeric(1)),
    vapply(error_runs,      function(x) x$error_count, numeric(1))
  )
  
  # collect results from successful runs
  final_results <- lapply(successful_runs, function(x) x$result)
  
  # return summary
  list(
    results           = final_results,
    total_error_count = total_error_count,
    n_errors          = length(error_runs),
    error_messages    = lapply(error_runs, function(x) x$error_message)
  )
}




run_one_iteration <- function(n_draws, burnin, n_obs, k, p, m, h,
                              intercept, p_bvar, alpha, lag_mean,
                              dummy_pars, mh_params, sim_params,
                              max_tries = 5) {
  #' Performs one iteration of the Monte Carlo simulation with error handling.
  #'
  #' Parameters:
  #' - n_draws (integer): Number of posterior draws to use in the simulation.
  #' - burnin (integer): Number of initial draws to discard as burn-in.
  #' - n_obs (integer): Number of observations used for estimation (excluding the forecast horizon).
  #' - k (integer): Number of endogenous variables in the VAR model.
  #' - p (integer): Lag order used for data simulation.
  #' - m (integer): Total number of parameters.
  #' - h (integer): Forecast horizon (number of steps ahead).
  #' - intercept (logical): Whether an intercept is included in the VAR model.
  #' - p_bvar (integer): Lag order used in VAR estimation.
  #' - alpha (numeric): Significance level for forecast prediction intervals.
  #' - lag_mean (numeric): Lag mean parameter used in the Minnesota prior.
  #' - dummy_pars (list): Parameters for dummy observations (e.g., mu, gamma, prior_mean).
  #' - mh_params (list): Parameters controlling the Metropolis-Hastings sampler in the hierarchical VAR.
  #' - sim_params (list): Parameters for data simulation (e.g., shock intensities, volatility periods).
  #' - max_tries (integer, default = 5): Maximum number of attempts to run one simulation iteration.
  #'
  #' Returns:
  #' - A list containing:
  #'    - status: "success" if the simulation iteration succeeded, "error" otherwise.
  #'    - result: The simulation result (only if status is "success").
  #'    - error_message: Error message (only if status is "error").
  #'    - error_count: The number of errors that occurred before success (or total errors if failed).
  
  local_error_count <- 0
  
  for (attempt in seq_len(max_tries)) {
    this_try <- tryCatch(
      {
        # attempt one simulation iteration
        res <- monte_carlo_simulation(
          n_draws    = n_draws,
          burnin     = burnin,
          n_obs      = n_obs,
          k          = k,
          p          = p,
          m          = m,
          h          = h,
          intercept  = intercept,
          p_bvar     = p_bvar,
          alpha      = alpha,
          lag_mean   = lag_mean,
          dummy_pars = dummy_pars,
          mh_params  = mh_params,
          sim_params = sim_params
        )
        
        return(list(
          status = "success",
          result = res,
          error_count = local_error_count
        ))
      },
      error = function(e) {
        # if error, we increase the local_error_count
        local_error_count <<- local_error_count + 1
        
        # if not reached max_tries, return NULL. If reached the last attempt, return an error object
        if (attempt < max_tries) {
          return(NULL)
        } else {
          return(list(
            status = "error",
            error_message = e$message,
            error_count   = local_error_count
          ))
        }
      }
    )
    
    # we stop trying further and return whatever it gave us.
    if (!is.null(this_try)) {
      return(this_try)
    }
  }
  
  # safeguard
  return(list(
    status = "error",
    error_message = "Unknown error after max_tries",
    error_count = local_error_count
  ))
}


monte_carlo_simulation <- function(n_draws, burnin, n_obs, k, p, m, h, intercept, 
                                   p_bvar, alpha, 
                                   lag_mean, dummy_pars, mh_params, sim_params){
  #' This function simulates time series data, fits several BVAR models using different prior 
  #' specifications and estimation methods (including classic Minnesota with inverse Wishart, 
  #' hierarchical Minnesota, SSVS, Minnesota with flat priors, and a frequentist VAR), generates 
  #' forecasts from each model, and computes forecast evaluation metrics (RMSE and prediction 
  #' interval accuracy). The simulation is run for a specified number of observations, forecast 
  #' horizon, and posterior draws.
  #'
  #' Parameters:
  #' - n_draws (integer): Total number of posterior draws to be used in the Monte Carlo simulation.
  #' - burnin (integer): Number of initial posterior draws to discard as burn-in.
  #' - n_obs (integer): Number of observations in the simulated data (excluding forecast horizon).
  #' - k (integer): Number of endogenous variables in the VAR model.
  #' - p (integer): Lag order for data generation.
  #' - m (integer): The total number of parameters in the VAR model.
  #' - h (integer): Forecast horizon (number of steps ahead to forecast).
  #' - intercept (logical): Indicator whether to include an intercept in the VAR model.
  #' - p_bvar (integer): Lag order to be used for the VAR estimation (may be equal to p).
  #' - alpha (numeric): Significance level for forecast prediction intervals (e.g., 0.05 for 95% intervals).
  #' - lag_mean (numeric): Lag mean parameter used in the Minnesota prior.
  #' - dummy_pars (list): Parameters for generating dummy observations, including mu, gamma, and 
  #'   (optionally) a prior mean.
  #' - mh_params (list): Parameters controlling the Metropolis-Hastings sampler in the hierarchical 
  #'   Minnesota approach (e.g., thinning interval, scaling factors, burn-in adjustment, acceptance rate thresholds).
  #' - sim_params (list): Parameters for data simulation, including minimum and maximum individual shocks, 
  #'   high volatility periods, and exogenous shocks.
  #'
  #' Returns:
  #' - A list containing forecast evaluation results for different methods:
  #'    - gibbs: Results for classic Minnesota with inverse Wishart.
  #'    - mh: Results for the hierarchical Minnesota approach.
  #'    - ssvs: Results for the SSVS (stochastic search variable selection) approach.
  #'    - flat: Results for the classic Minnesota approach with flat priors.
  #'    - var: Results for the frequentist VAR.

  # ----------------------------------------------------------------------------
  # generating data
  # ----------------------------------------------------------------------------
  
  T_data <- n_obs + h
  data <- sim_y(k, p, T_data, 
                d_min = sim_params$d_min,
                d_max = sim_params$d_max,
                off_d_min = sim_params$off_d_min,
                off_d_max = sim_params$off_d_max,
                init_sd = sim_params$init_sd,
                shock_diag_min = sim_params$shock_diag_min,
                shock_diag_max = sim_params$shock_diag_max,
                mean_vola = sim_params$mean_vola,
                sd_vola = sim_params$sd_vola,
                min_indiv_shocks = sim_params$min_indiv_shocks,
                max_indiv_shocks = sim_params$max_indiv_shocks, 
                min_high_vol_periods = sim_params$min_high_vol_periods,
                max_high_vol_periods = sim_params$max_high_vol_periods,
                min_exog_shocks = sim_params$min_exog_shocks,
                max_exog_shocks = sim_params$max_exog_shocks
                )
  Ypred <- data[(T_data-h+1):T_data,]
  Yraw <- data[1:(T_data-h),]
  
  # ----------------------------------------------------------------------------
  # draw from posterior predictive distribution n_draws times
  # ----------------------------------------------------------------------------
  
  # ---- rigorous methods ------------------------------------------------------
  # -- classic minnesota with inverse wishart -------------
  # posterior coefficient draws
  res_gibbs <- run_bvar_minnesota(
    Yraw = Yraw,
    p = p_bvar,
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
  fcst_gibbs <- bvar_forecast(Y = Yraw, 
                              phi_draws = res_gibbs$Phi, 
                              Sigma_draws = res_gibbs$Sigma, 
                              p = p_bvar, 
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
                         mode = hhelper_compute_sigma_vec(Yraw, p_bvar)
                       )
  )
  
  # posterior coefficient draws
  res_mh <- run_bvar_hierarch(
    Yraw = Yraw, 
    p = p_bvar, 
    intercept = intercept,
    mh_params = mh_params,
    hyper_params = hyper_params,
    pi4 = 1000,
    lag_mean = lag_mean, 
    n_draws = n_draws, 
    burnin = burnin
  )
  
  # forecasts and evaluation
  fcst_mh <- bvar_forecast(Y = Yraw, 
                           phi_draws = res_mh$Phi, 
                           Sigma_draws = res_mh$Sigma, 
                           p = p_bvar, 
                           h = h, 
                           alpha = alpha)
  rmse_mh <- fcst_rmse(fcst_obj = fcst_mh, Ypred = Ypred, h = h, k = k)
  pred_acc_mh <- fcst_pred_int_acc(fcst_obj = fcst_mh, Ypred = Ypred, h = h)
  
  
  # -- ssvs -----------------------------------------------
  # posterior coefficient draws
  res_ssvs <- run_bvar_ssvs(
    Yraw, p_bvar, 
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
  fcst_ssvs <- bvar_forecast(Y = Yraw, 
                             phi_draws = res_ssvs$Phi, 
                             Sigma_draws = res_ssvs$Sigma, 
                             p = p_bvar, 
                             h = h, 
                             alpha = alpha)
  rmse_ssvs <- fcst_rmse(fcst_obj = fcst_ssvs, Ypred = Ypred, h = h, k = k)
  pred_acc_ssvs <- fcst_pred_int_acc(fcst_obj = fcst_ssvs, Ypred = Ypred, h = h)
  
  # ---- baseline methods ------------------------------------------------------
  # -- classic minnesota with flat priors -----------------
  res_flat <- run_bvar_minnesota(
    Yraw = Yraw,
    p = p_bvar,
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
  fcst_flat <- bvar_forecast(Y = Yraw, 
                             phi_draws = res_flat$Phi, 
                             Sigma_draws = res_flat$Sigma, 
                             p = p_bvar, 
                             h = h, 
                             alpha = alpha)
  rmse_flat <- fcst_rmse(fcst_obj = fcst_flat, Ypred = Ypred, h = h, k = k)
  pred_acc_flat <- fcst_pred_int_acc(fcst_obj = fcst_flat, Ypred = Ypred, h = h)
  
  
  # -- frequentist VAR ------------------------------------
  # compute coefficients
  colnames(Yraw) <- paste0("Y", 1:k)
  res_var <- VAR(Yraw, p = p_bvar, type = "const")
  
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


evaluate_sim_res <- function(result_obj, models, h, n_iter){
  #' Evaluates simulation results across multiple models and iterations.
  #'
  #' Parameters:
  #' - result_obj (list): A list of length n_iter, where each element is itself a list containing
  #'   the evaluation results for each model. Each model's result should include:
  #'     - rmse: A list with elements 'all_rmse' (overall RMSE) and 'row_rmse' (h-step ahead RMSE).
  #'     - pred_acc: A matrix of prediction interval accuracy (dimensions h x 1).
  #' - models (character vector): A vector of model names corresponding to the different models evaluated.
  #' - h (integer): The forecast horizon (number of steps ahead).
  #' - n_iter (integer): The number of simulation iterations.
  #'
  #' Returns:
  #' - A list containing:
  #'     - all_rmse: A 1 x total_model matrix of overall RMSE means for each model.
  #'     - h_fcst_rmse: An h x total_model matrix of h-step ahead RMSE means for each model.
  #'     - pred_acc: An h x total_model matrix of mean prediction interval accuracy for each model.
  
  total_model <- length(models)
  all_rmse <- matrix(0, nrow = total_model, ncol = n_iter)
  h_fcst_rmse <- matrix(0, nrow = h * total_model, ncol = n_iter)
  pred_acc <- matrix(0, nrow = (h+1) * total_model, ncol = n_iter)
  
  # get results in matrix form
  for (model in seq_len(total_model)) {
    for (iter in seq_len(n_iter)) {
      all_rmse[model,iter] <- result_obj$results[[iter]][[model]]$rmse$all_rmse
      h_fcst_rmse[((model - 1) * h + 1):(model * h), iter] <- result_obj$results[[iter]][[model]]$rmse$row_rmse
      pred_acc[((model - 1) * (h+1) + 1):(model * (h+1)), iter] <- result_obj$results[[iter]][[model]]$pred_acc
    }
  }
  
  # overall mean
  all_rmse_mean <- matrix(rowMeans(all_rmse), nrow = 1)
  colnames(all_rmse_mean) <- models
  
  # T+h means
  h_fcst_rmse_mean <- matrix(
    rowMeans(h_fcst_rmse), 
    nrow = h, 
    ncol = total_model, 
    dimnames = list(paste0("T+", 1:h), models))
  
  # prediction interval accuracy
  pred_acc_mean <- matrix(
    rowMeans(pred_acc), 
    nrow = h+1, 
    ncol = total_model, 
    dimnames = list(c(paste0("T+", 1:h), "overall_acc"), models))
  
  return(list(
    all_rmse = drop(all_rmse_mean),
    h_fcst_rmse = h_fcst_rmse_mean,
    pred_acc = pred_acc_mean
  ))
}
