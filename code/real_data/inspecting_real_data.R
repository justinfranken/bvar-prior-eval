# testing the behavior of real financial revenue growth for adjusting MC simulations.

# import sources
source(paste0(getwd(),"/downloading_real_data.r"))

# functions for analyzing real revenue data
get_industry_growth_rates <- function(tickers, corp_names, type){
  #' Loop over input tickers to download industry growth rates
  #' 
  #' Parameters:
  #' - tickers (vector): Containing the tickers of the corporates to be downloaded in upper cases.
  #' - corp_names (vector): Containing the names of the corporates in lower cases.
  #' - type (scalar): Number indicating which growth rate is downloaded:
  #'  - 2: quarterly growth rates.
  #'  - 3: 12 months quarterly growth rates.
  #'  - 4: TTM quarterly growth rates.
  #'  - 5: TTM 12 months growth rates.
  #' 
  #' Returns:
  #' - Data.frame of industry growth rates of determined type.
  
  # download revenue data and store in list
  df_list <- list()
  for (i in 1:length(tickers)) {
    df_list[[i]] <- get_growth_rates(download_revenue_data(
      tickers[i], 
      corp_names[i]))[,c(1,type)]
    colnames(df_list[[i]]) = c("date", tickers[i])
  }
  
  # merge elements
  out_df <- Reduce(function(x, y) merge(x, y, by = "date"), df_list)
  return(out_df)
}


plot_gr <- function(industry_gr){
  #' Plot industry growth rates.
  #' 
  #' Parameters:
  #' - industry_gr (data.frame): Output of get_industry_growth_rates().
  #' 
  #' Returns:
  #' - Plot of the input growth rates.
  matplot(
    industry_gr[,1], industry_gr[,-1], type = "l", lty = 1, 
    xlab = "Year", ylab = "Growth Rate", 
    main = "Growth Rates of Companies Over Time"
  )
  legend(
    "topleft", 
    legend = colnames(industry_gr[,-1]), 
    col = 1:ncol(industry_gr[,-1]), 
    lty = 1, 
    cex = 0.7
  )
  abline(h = 0, col = "black", lty = 2)
}


# --- large cap pharmaceutical -------------------------------------------------
pharma_tickers <- c("LLY", "NVO", "JNJ", "ABBV", "MRK", "AZN", "NVS", "PFE")
pharma_corp_names <- c("eli-lilly", 
                       "novo-nordisk", 
                       "johnson_johnson", 
                       "abbvie", "merck", 
                       "astrazeneca", 
                       "novartis-ag", 
                       "pfizer")
pharma <- get_industry_growth_rates(pharma_tickers, pharma_corp_names, 4)
plot_gr(pharma)

# summary statistics
summary(pharma[,-1])
mean(colMeans(as.matrix(pharma[,-1]), na.rm = TRUE))
abline(h = mean(colMeans(as.matrix(pharma[,-1]), na.rm = TRUE)), col = "red", lty = 2)

# check time series
pharma_ts <- ts(na.omit(pharma[,-1]), start = c(2012, 1), frequency = 4)
colnames(pharma_ts) <- pharma_tickers
# stationarity
pharma_adf <- matrix(NA, ncol = length(pharma_tickers), nrow = 1)
for (i in 1:length(pharma_tickers)) {
  pharma_adf[1,i] <- adf.test(pharma_ts[,i])$p.value
}
colnames(pharma_adf) <- pharma_tickers
pharma_adf
# lags
acf(pharma_ts)
VARselect(pharma_ts, type = "const", lag.max = 10)
# VAR coefficients
pharma_var <- VAR(pharma_ts[,c("LLY", "JNJ", "MRK", "AZN", "NVS")], p = 5, type = "const")
summary(pharma_var)


