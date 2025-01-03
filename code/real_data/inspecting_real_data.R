# testing the behavior of real financial revenue growth for adjusting MC simulations.

# import sources
source(paste0(getwd(),"/downloading_real_data.r"))

# functions for analyzing real revenue data
get_industry_growth_rates <- function(tickers, corp_names){
  #' Loop over input tickers to download adjusted industry growth rates.
  #' 
  #' Parameters:
  #' - tickers (vector): Containing the tickers of the corporates to be downloaded in upper cases.
  #' - corp_names (vector): Containing the names of the corporates in lower cases.
  #' 
  #' Returns:
  #' - Data.frame of adjusted industry growth rates.
  
  # download revenue data and store in list
  df_list <- list()
  for (i in 1:length(tickers)) {
    df_list[[i]] <- transform_to_stationary(download_revenue_data(
      tickers[i], 
      corp_names[i]))
    colnames(df_list[[i]]) = c("date", tickers[i])
    cat(paste0("Downloaded ", tickers[i], ", which is of row ", nrow(df_list[[i]]), " and col ", ncol(df_list[[i]]), "\n"))
  }
  
  # merge elements
  out_df <- Reduce(function(x, y) merge(x, y, by = "date"), df_list)
  return(out_df)
}


plot_gr <- function(industry_gr, adj_cex = 0.7, title = "Growth Rates of Companies Over Time", y_lim = c(-0.35, 0.35), legend_off = FALSE){
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
  matplot(
    industry_gr[,1], industry_gr[,-1], type = "l", lty = 1, 
    col = 1:ncol(industry_gr[,-1]),
    ylim = y_lim,
    xlab = "Year", ylab = "Growth Rate", 
    main = title
  )
  abline(h = mean(colMeans(as.matrix(industry_gr[,-1]), na.rm = TRUE)), col = "red", lty = 2)
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


# --- Utility - Electrical Power Distribution ----------------------------------
utilities_tickers <- c("NEE", "SO", "DUK", "PCG", "AEP", "D", "PEG")
utilities_corp_names <- c("nextera-energy", 
                       "southern", 
                       "duke-energy", 
                       "pacific-gas-electric",
                       "american-electric-power", 
                       "dominion-energy", 
                       "public-service-enterprise-group")
utilities <- get_industry_growth_rates(utilities_tickers, utilities_corp_names)
plot_gr(utilities, title = "Utilities")

# summary statistics
summary(utilities[,-1])
mean(colMeans(as.matrix(utilities[,-1]), na.rm = TRUE))
abline(h = mean(colMeans(as.matrix(utilities[,-1]), na.rm = TRUE)), col = "red", lty = 2)

# check time series
utilities_ts <- ts(na.omit(utilities[,-1]), start = c(2009, 2), frequency = 4)
colnames(utilities_ts) <- utilities_tickers
# stationarity
utilities_adf <- matrix(NA, ncol = length(utilities_tickers), nrow = 1)
for (i in 1:length(utilities_tickers)) {
  utilities_adf[1,i] <- adf.test(utilities_ts[,i])$p.value
}
colnames(utilities_adf) <- utilities_tickers
utilities_adf
# lags
acf(utilities_ts)
VARselect(utilities_ts, type = "none", lag.max = 10)
# VAR coefficients
utilities_var <- VAR(utilities_ts[,c("SO", "NEE", "AEP", "PCG", "PEG")], p = 8, type = "none")
summary(utilities_var)


# --- computer software --------------------------------------------------------
soft_tickers <- c("MSFT", "SAP", "CDNS", "MSTR", "ANSS", "PTC", "SSNC", "MANH", "PEGA")
soft_corp_names <- c("microsoft",
                       "sap-se",
                       "cadence-design-systems",
                       "microstrategy",
                       "ansys",
                       "ptc",
                       "ss-c-technologies-holdings",
                       "manhattan-associates",
                       "pegasystems"
                       )
soft <- get_industry_growth_rates(soft_tickers, soft_corp_names)
plot_gr(soft, title = "Software")

# summary statistics
summary(soft[,-1])
mean(colMeans(as.matrix(soft[,-1]), na.rm = TRUE))
abline(h = mean(colMeans(as.matrix(soft[,-1]), na.rm = TRUE)), col = "red", lty = 2)

# check time series
soft_ts <- ts(na.omit(soft[,-1]), start = c(2010, 1), frequency = 4)
colnames(soft_ts) <- soft_tickers
# stationarity
soft_adf <- matrix(NA, ncol = length(soft_tickers), nrow = 1)
for (i in 1:length(soft_tickers)) {
  soft_adf[1,i] <- adf.test(soft_ts[,i])$p.value
}
colnames(soft_adf) <- soft_tickers
soft_adf
# lags
acf(soft_ts)
VARselect(soft_ts, type = "none", lag.max = 10)
# VAR coefficients
soft_var <- VAR(soft_ts[,c("MSFT", "SAP", "MSTR", "SSNC", "PEGA")], p = 6, type = "none")
summary(soft_var)


# --- oil & gas --------------------------------------------------------
oag_tickers <- c("XOM", "CVX", "SHEL", "BP", "E", "YPF")
oag_corp_names <- c("exxon",
                     "chevron",
                     "shell",
                     "bp",
                     "eni-spa",
                     "ypf-sociedad-anonima"
)
oag <- get_industry_growth_rates(oag_tickers, oag_corp_names)
plot_gr(oag, title = "Oil and Gas")

# summary statistics
summary(oag[,-1])
mean(colMeans(as.matrix(oag[,-1]), na.rm = TRUE))
abline(h = mean(colMeans(as.matrix(oag[,-1]), na.rm = TRUE)), col = "red", lty = 2)

# check time series
oag_ts <- ts(na.omit(oag[,-1]), start = c(2010, 1), frequency = 4)
colnames(oag_ts) <- oag_tickers
# stationarity
oag_adf <- matrix(NA, ncol = length(oag_tickers), nrow = 1)
for (i in 1:length(oag_tickers)) {
  oag_adf[1,i] <- adf.test(oag_ts[,i])$p.value
}
colnames(oag_adf) <- oag_tickers
oag_adf
# lags
acf(oag_ts)
VARselect(oag_ts, type = "none", lag.max = 10)
# VAR coefficients
oag_var <- VAR(oag_ts[,c("XOM", "CVX", "SHEL", "BP", "YPF")], p = 9, type = "none")
summary(oag_var)


# --- Medical Products Manufacturers -------------------------------------------
med_prod_tickers <- c("ABT", "SYK", "BSX", "RMD", "ZBH", "PODD")
med_prod_corp_names <- c("abbott-laboratories",
                     "stryker",
                     "boston-scientific",
                     "resmed",
                     "zimmer-biomet-holdings",
                     "insulet"
)
med_prod <- get_industry_growth_rates(med_prod_tickers, med_prod_corp_names)
plot_gr(med_prod, title = "Medical Products_Manufacturers")

# summary statistics
summary(med_prod[,-1])
mean(colMeans(as.matrix(med_prod[,-1]), na.rm = TRUE))
abline(h = mean(colMeans(as.matrix(med_prod[,-1]), na.rm = TRUE)), col = "red", lty = 2)

# check time series
med_prod_ts <- ts(na.omit(med_prod[,-1]), start = c(2010, 1), frequency = 4)
colnames(med_prod_ts) <- med_prod_tickers
# stationarity
med_prod_adf <- matrix(NA, ncol = length(med_prod_tickers), nrow = 1)
for (i in 1:length(med_prod_tickers)) {
  med_prod_adf[1,i] <- adf.test(med_prod_ts[,i])$p.value
}
colnames(med_prod_adf) <- med_prod_tickers
med_prod_adf
# lags
acf(med_prod_ts)
VARselect(med_prod_ts, type = "none", lag.max = 10)
# VAR coefficients
med_prod_var <- VAR(med_prod_ts[,c("ABT", "SYK", "BSX", "RMD", "ZBH")], p = 9, type = "none")
summary(med_prod_var)


# --- consumer products --------------------------------------------------------
cp_tickers <- c("PG", "CL", "KMB", "CHD", "CLX", "NWL", "ENR")
cp_corp_names <- c("procter-gamble",
                     "colgate-palmolive",
                     "kimberly-clark",
                     "church-dwight",
                     "clorox",
                     "newell-brands",
                     "energizer-holdings"
)
cp <- get_industry_growth_rates(cp_tickers, cp_corp_names)
plot_gr(cp[,c(-7, -8)], title = "Commercial Products Manufacturers")

# summary statistics
summary(cp[,-1])
mean(colMeans(as.matrix(cp[,-1]), na.rm = TRUE))
abline(h = mean(colMeans(as.matrix(cp[,-1]), na.rm = TRUE)), col = "red", lty = 2)

# check time series
cp_ts <- ts(na.omit(cp[,-1]), start = c(2010, 1), frequency = 4)
colnames(cp_ts) <- cp_tickers
# stationarity
cp_adf <- matrix(NA, ncol = length(cp_tickers), nrow = 1)
for (i in 1:length(cp_tickers)) {
  cp_adf[1,i] <- adf.test(cp_ts[,i])$p.value
}
colnames(cp_adf) <- cp_tickers
cp_adf
# lags
acf(cp_ts)
VARselect(cp_ts, type = "none", lag.max = 10)
# VAR coefficients
cp_var <- VAR(cp_ts[,c("PG", "KMB", "CHD", "CLX", "NWL")], p = 8, type = "none")
summary(cp_var)


# --- general machinery --------------------------------------------------------
machinary_tickers <- c("PH", "ITW", "ATLKY", "DOV", "TRMB", "IEX", "GGG", "RRX")
machinary_corp_names <- c("parker-hannifin",
                   "illinois-tool-works",
                   "atlas-copco-ab",
                   "dover",
                   "trimble",
                   "idex",
                   "graco",
                   "regal-rexnord"
)
machinary <- get_industry_growth_rates(machinary_tickers, machinary_corp_names)
plot_gr(machinary, title = "Machinery Producers")

# summary statistics
summary(machinary[,-1])
mean(colMeans(as.matrix(machinary[,-1]), na.rm = TRUE))
abline(h = mean(colMeans(as.matrix(machinary[,-1]), na.rm = TRUE)), col = "red", lty = 2)

# check time series
machinary_ts <- ts(na.omit(machinary[,-1]), start = c(2010, 1), frequency = 4)
colnames(machinary_ts) <- machinary_tickers
# stationarity
machinary_adf <- matrix(NA, ncol = length(machinary_tickers), nrow = 1)
for (i in 1:length(machinary_tickers)) {
  machinary_adf[1,i] <- adf.test(machinary_ts[,i])$p.value
}
colnames(machinary_adf) <- machinary_tickers
machinary_adf
# lags
acf(machinary_ts)
VARselect(machinary_ts, type = "none", lag.max = 10)
# VAR coefficients
machinary_var <- VAR(machinary_ts[,c("ITW", "ATLKY", "DOV", "TRMB", "GGG")], p = 7, type = "none")
summary(machinary_var)


# --- american banks -----------------------------------------------------------
banks_tickers <- c("JPM", "BAC", "WFC", "MS", "GS", "C", "SCHW")
banks_corp_names <- c("jpmorgan-chase",
                          "bank-of-america",
                          "wells-fargo",
                          "morgan-stanley",
                          "goldman-sachs",
                          "citigroup",
                          "charles-schwab"
)
banks <- get_industry_growth_rates(banks_tickers, banks_corp_names)
plot_gr(banks, title = "Banks")

# summary statistics
summary(banks[,-1])
mean(colMeans(as.matrix(banks[,-1]), na.rm = TRUE))
abline(h = mean(colMeans(as.matrix(banks[,-1]), na.rm = TRUE)), col = "red", lty = 2)

# check time series
banks_ts <- ts(na.omit(banks[,-1]), start = c(2010, 1), frequency = 4)
colnames(banks_ts) <- banks_tickers
# stationarity
banks_adf <- matrix(NA, ncol = length(banks_tickers), nrow = 1)
for (i in 1:length(banks_tickers)) {
  banks_adf[1,i] <- adf.test(banks_ts[,i])$p.value
}
colnames(banks_adf) <- banks_tickers
banks_adf
# lags
acf(banks_ts)
VARselect(banks_ts, type = "none", lag.max = 10)
# VAR coefficients
banks_var <- VAR(banks_ts[,c("JPM", "WFC", "C", "BAC", "MS")], p = 8, type = "none")
summary(banks_var)


# --- get overall average mean and its standard deviation ----------------------
list_all_industries <- list(med_prod, banks, cp, machinary, oag, utilities, soft)
merged_all_industry <- Reduce(function(x, y) merge(x, y, by = "date"), list_all_industries)
mean(colMeans(as.matrix(merged_all_industry[,-1]), na.rm = TRUE))
sd(colMeans(as.matrix(merged_all_industry[,-1]), na.rm = TRUE))


# --- plot all industry growth rates -------------------------------------------
par(mfrow = c(4,2))
indust_list <- list(utilities, soft[,-8], oag, med_prod, cp[,c(-7, -8)], machinary, banks)
indust_names <- c("utilities growth rates",
                  "software growth rates",
                  "oil and gas growth rates",
                  "med_prod manufacturers growth rates",
                  "comercial products growth rates",
                  "machinery growth rates",
                  "banks growth rates")
for (i in 1:length(indust_list)) {
  plot_gr(indust_list[[i]], title = indust_names[i], legend_off = TRUE)
}
par(mfrow = c(1,1))




y <- abbvie$quarterly_gr

min_y <- min(y, na.rm = TRUE)  # Find the minimum value
shift <- ifelse(min_y <= 0, abs(min_y) + 1, 0)  # Determine the shift amount
shifted_y <- y + shift
log_y <- log(shifted_y)
seasonally_differenced <- diff(log_y, lag = 4)
stationary_series <- diff(seasonally_differenced)

abbvie$log_qrt_gr <- log_y
abbvie$log_diffed_qrt_gr <- c(rep(NA, 5), stationary_series)
