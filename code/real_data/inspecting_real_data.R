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
    cat(paste0("Downloaded ", tickers[i], ", which is of row ", nrow(df_list[[i]]), " and col ", ncol(df_list[[i]]), "\n"))
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
    col = 1:ncol(industry_gr[,-1]),
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
soft <- get_industry_growth_rates(soft_tickers, soft_corp_names, 4)
plot_gr(soft)

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
VARselect(soft_ts, type = "const", lag.max = 10)
# VAR coefficients
soft_var <- VAR(soft_ts[,c("MSFT", "SAP", "MSTR", "SSNC", "PEGA")], p = 6, type = "const")
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
oag <- get_industry_growth_rates(oag_tickers, oag_corp_names, 4)
plot_gr(oag)

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
VARselect(oag_ts, type = "const", lag.max = 10)
# VAR coefficients
oag_var <- VAR(oag_ts[,c("XOM", "CVX", "SHEL", "BP", "YPF")], p = 8, type = "const")
summary(oag_var)


# --- automotive ---------------------------------------------------------------
auto_tickers <- c("F", "GM", "TM", "BAMXF", "VWAGY", "MBGYY", "TSLA", "HOG")
auto_corp_names <- c("ford-motor",
                     "general-motors",
                     "toyota",
                     "bmw",
                     "volkswagen-ag",
                     "mercedes-benz-group-ag",
                     "tesla",
                     "harley-davidson"
)
auto <- get_industry_growth_rates(auto_tickers, auto_corp_names, 4)
plot_gr(auto[,-8])

# summary statistics
summary(auto[,-1])
mean(colMeans(as.matrix(auto[,-1]), na.rm = TRUE))
abline(h = mean(colMeans(as.matrix(auto[,-1]), na.rm = TRUE)), col = "red", lty = 2)

# check time series
auto_ts <- ts(na.omit(auto[,-1]), start = c(2010, 1), frequency = 4)
colnames(auto_ts) <- auto_tickers
# stationarity
auto_adf <- matrix(NA, ncol = length(auto_tickers), nrow = 1)
for (i in 1:length(auto_tickers)) {
  auto_adf[1,i] <- adf.test(auto_ts[,i])$p.value
}
colnames(auto_adf) <- auto_tickers
auto_adf
# lags
acf(auto_ts)
VARselect(auto_ts, type = "const", lag.max = 10)
# VAR coefficients
auto_var <- VAR(auto_ts[,c("GM", "TM", "BAMXF", "VWAGY", "MBGYY")], p = 6, type = "const")
summary(auto_var)


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
cp <- get_industry_growth_rates(cp_tickers, cp_corp_names, 4)
plot_gr(cp[,-7])

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
VARselect(cp_ts, type = "const", lag.max = 10)
# VAR coefficients
cp_var <- VAR(cp_ts[,c("PG", "KMB", "CHD", "CLX", "NWL")], p = 7, type = "const")
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
machinary <- get_industry_growth_rates(machinary_tickers, machinary_corp_names, 4)
plot_gr(machinary)

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
VARselect(machinary_ts, type = "const", lag.max = 10)
# VAR coefficients
machinary_var <- VAR(machinary_ts[,c("ITW", "ATLKY", "DOV", "TRMB", "GGG")], p = 6, type = "const")
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
banks <- get_industry_growth_rates(banks_tickers, banks_corp_names, 4)
plot_gr(banks)

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
VARselect(banks_ts, type = "const", lag.max = 10)
# VAR coefficients
banks_var <- VAR(banks_ts[,c("JPM", "WFC", "C", "BAC", "MS")], p = 7, type = "const")
summary(banks_var)


# --- get overall average mean and its standard deviation ----------------------
list_all_industries <- list(auto, banks, cp, machinary, oag, pharma, soft)
merged_all_industry <- Reduce(function(x, y) merge(x, y, by = "date"), list_all_industries)
mean(colMeans(as.matrix(merged_all_industry[,-1]), na.rm = TRUE))
sd(colMeans(as.matrix(merged_all_industry[,-1]), na.rm = TRUE))


