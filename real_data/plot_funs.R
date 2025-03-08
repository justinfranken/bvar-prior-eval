# Functions for plotting real data.


plot_gr <- function(industry_gr, 
                    adj_cex = 0.7,
                    title = "Growth Rates of Companies Over Time", 
                    y_lim = c(-35, 35), 
                    legend_off = FALSE){
  #' Plot industry growth rates.
  #' 
  #' Parameters:
  #' - industry_gr (data.frame): Output of get_industry_growth_rates().
  #' - adj_cex (scalar): Adjust size of legend in plot.
  #' - title (string): Title of the plot.
  #' - legend_off (boolean): If TRUE, legend won't be displayed.
  #' 
  #' Returns:
  #' - Plot of the input growth rates.
  if(title == "off"){
    matplot(
      industry_gr[,1], industry_gr[,-1], type = "l", lty = 1, 
      xlab = "Year", ylab = "Growth Rate",
      col = 1:ncol(industry_gr[,-1]),
      ylim = y_lim,
      cex.axis = 1.5
    )  
  } else{
    matplot(
      industry_gr[,1], industry_gr[,-1], type = "l", lty = 1, 
      col = 1:ncol(industry_gr[,-1]),
      xlab = "Year", ylab = "Growth Rate", 
      main = title,
      ylim = y_lim,
      cex.axis = 1.5
    ) 
  }
  if (legend_off == FALSE) {
    legend(
      "topleft", 
      legend = colnames(industry_gr[,-1]), 
      col = 1:ncol(industry_gr[,-1]), 
      lty = 1, 
      cex = adj_cex
    )
  }
  abline(h = 0, col = "black", lty = 2)
}



plot_pred_intervals <- function(data_mat, data_df, p, h, alpha, industry, y_lim, main = TRUE){
  #' Plots forecast prediction intervals for the first column firm and industry, training with the initial 40
  #' observations and predicting h-steps ahead.
  #'
  #' Parameters:
  #' - data_mat (matrix): A numeric data matrix with dimensions T x k. The first column corresponds to the focal firm (ticker).
  #' - data_df (data.frame): A data frame that contains a date column named "date" used for plotting the x-axis.
  #' - p (integer): The lag order used in the BVAR models.
  #' - h (integer): The forecast horizon.
  #' - alpha (numeric): The significance level for constructing prediction intervals (e.g., 0.05 for 95% intervals).
  #' - industry (character): A string specifying the industry name, used in the plot title and legend.
  #' - y_lim (numeric vector): A two-element vector specifying the y-axis limits for the plot.
  #' - main (logical): If set to FALSE, no plot title is used.
  #'
  #' Returns:
  #' - A plot is generated showing the actual time series and the forecasted median values along with the
  #'   upper and lower prediction intervals for each of the following methods: NIW (classic Minnesota with
  #'   inverse Wishart), Hierarchical Minnesota, SSVS, Flat Minnesota, and a frequentist VAR.
  
  k <- ncol(data_mat)
  ticker <- colnames(data_mat)[1]
  if(main){
    main <- paste0(industry, " Industry Prediction of ", ticker)
  } else{
    main <- ""
  }
  
  # ---- NIW ---------------------------------------------------------------------
  res_niw <- run_bvar_minnesota(data_mat[1:40,], 
                                p, 
                                intercept = TRUE, 
                                use_dummies = FALSE, 
                                n_draws = 5000, 
                                burnin = 1000)
  fcst_niw <- bvar_forecast(data_mat[1:40,], 
                            res_niw$Phi, 
                            res_niw$Sigma, 
                            p, h, alpha)
  
  # ---- Hierarch ---------------------------------------------------------------------
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
                                            mode = hhelper_compute_sigma_vec(data_mat[1:40,], p)))
  res_hierarch <- run_bvar_hierarch(data_mat[1:40,], 
                                    p, 
                                    intercept = TRUE, 
                                    hyper_params = hyper_params,
                                    mh_params = mh_params, 
                                    n_draws = 5000, 
                                    burnin = 1000)
  fcst_hierarch <- bvar_forecast(data_mat[1:40,], 
                                 res_hierarch$Phi, 
                                 res_hierarch$Sigma, 
                                 p, h, alpha)
  
  # ---- SSVS ---------------------------------------------------------------------
  res_ssvs <- run_bvar_ssvs(data_mat[1:40,], 
                            p, 
                            intercept = TRUE, 
                            use_dummies = FALSE,
                            tau0 = 1/100,
                            tau1 = 1,
                            delta_prob = 0.8,
                            n_draws = 5000, 
                            burnin = 1000)
  fcst_ssvs <- bvar_forecast(data_mat[1:40,], 
                             res_ssvs$Phi, 
                             res_ssvs$Sigma, 
                             p, h, alpha)
  
  # ---- Flat ---------------------------------------------------------------------
  res_flat <- run_bvar_minnesota(data_mat[1:40,], 
                                 use_flat = TRUE,
                                 p, 
                                 intercept = TRUE, 
                                 use_dummies = FALSE, 
                                 n_draws = 5000, 
                                 burnin = 1000)
  fcst_flat <- bvar_forecast(data_mat[1:40,], 
                             res_flat$Phi, 
                             res_flat$Sigma, 
                             p, h, alpha)
  
  # ---- VAR ---------------------------------------------------------------------
  res_var <- VAR(data_mat[1:40,], p = p, type = "const")
  temp_var <- predict(res_var, n.ahead = h, ci = 1 - alpha)
  fcst_var <- list(
    median = sapply(temp_var$fcst, function(x) x[,1]),
    lower  = sapply(temp_var$fcst, function(x) x[,2]),
    upper  = sapply(temp_var$fcst, function(x) x[,3])
  )
  
  
  # ---- create prediction plot - Microsoft --------------------------------------
  matplot(data_df$date[35:44], data_mat[35:44, ], type = "l", lty = 1,
          col = c("black", rep("grey", 6)),
          lwd = c(2, rep(1, 6)),
          main = main,
          ylim = y_lim,
          xlab = "Date", ylab = "Growth Rates",
          cex.axis = 1.5)
  abline(h = 0, col = "black", lty = 2)
  
  # add niw lines
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_niw$median[,1]), col = "red", lty = 1, lwd = 2)
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_niw$upper[,1]), col = "red", lty = 2, lwd = 2)
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_niw$lower[,1]), col = "red", lty = 2, lwd = 2)
  # add hierarch lines
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_hierarch$median[,1]), col = "blue", lty = 1, lwd = 2)
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_hierarch$upper[,1]), col = "blue", lty = 2, lwd = 2)
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_hierarch$lower[,1]), col = "blue", lty = 2, lwd = 2)
  # add ssvs lines
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_ssvs$median[,1]), col = "purple", lty = 1, lwd = 2)
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_ssvs$upper[,1]), col = "purple", lty = 2, lwd = 2)
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_ssvs$lower[,1]), col = "purple", lty = 2, lwd = 2)
  # add flat lines
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_flat$median[,1]), col = "orange", lty = 1, lwd = 2)
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_flat$upper[,1]), col = "orange", lty = 2, lwd = 2)
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_flat$lower[,1]), col = "orange", lty = 2, lwd = 2)
  # add var lines
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_var$median[,1]), col = "turquoise", lty = 1, lwd = 2)
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_var$upper[,1]), col = "turquoise", lty = 2, lwd = 2)
  lines(data_df$date[40:44], c(data_mat[40,1],fcst_var$lower[,1]), col = "turquoise", lty = 2, lwd = 2)
  
  legend(
    "topleft", 
    legend = c(ticker, paste0("Other ", industry, " firms"), "NIW", "Hierarch", "SSVS", "Flat", "VAR"), 
    col = c("black", "grey", "red", "blue", "purple", "orange", "turquoise"), 
    lty = rep(1,7),
    lwd = c(2,1,2,2,2,2,2),
    cex = 1
  )
}
