# functions to download real financial data from macrotrends.net and return stationary growth rates.


# ------------------------------------------------------------------------------
# functions to download real data
# ------------------------------------------------------------------------------

download_revenue_data <- function(ticker, corp_name){
  #' Downloads quarterly real financial revenue data from listed equities from macrotrends.
  #' 
  #' Parameters:
  #' - ticker (string): Ticker of the equity to be downloaded in upper case. E.g. 'AAPL'.
  #' - corp_name (string): Corporate name in lower case of that company. E.g. 'apple'.
  #' 
  #' Returns:
  #' - a data.frame containing 2 columns:
  #'  - date: string of the format 'YYYY-MM-DD'.
  #'  - revenue: Numerical revenues in mio USD.
  
  # get the data
  url <- paste0(
    "https://macrotrends.net/stocks/charts/", ticker, "/", corp_name, "/revenue"
    )
  URL <- GET(url,user_agent("Justin"))
  webpage <- read_html(URL)
  
  temp <- webpage %>%
    html_nodes(xpath = "//table[contains(@class, 'historical_data_table')]")
  
  quarterly_table <- temp[2] %>% 
    html_table(fill = TRUE)
  
  # shape downloaded data
  colnames(quarterly_table[[1]]) <- c(
    "date", 
    "revenue"
    )
  cleaned_data <- quarterly_table[[1]]
  
  cleaned_data$date <- as.Date(cleaned_data$date, format = "%Y-%m-%d")
  cleaned_data$revenue <- as.numeric(
    gsub("[$,]", "", cleaned_data$revenue)
    )
  cleaned_data <- cleaned_data[order(cleaned_data$date), ]

  # small break to prevent overuse of downloading capacity
  Sys.sleep(5)
  
  out <- as.data.frame(cleaned_data)
  return(out)
}


transform_to_stationary <- function(data){
  data$date <- as.Date(data$date)
  
  # Sort the data by date (just in case it's not already sorted)
  data <- data[order(data$date), ]
  
  # Log Transformation
  data$log_revenue <- log(data$revenue)

  # First-order Differencing (subtract the previous quarter's value)
  data$log_diff <- c(NA, diff(data$log_revenue, lag = 1))
  
  # Log + Seasonal Differencing
  data$log_seasonal_diff <- c(rep(NA, 4), diff(data$log_diff, lag = 4)) * 100
  
  # Remove rows with NA values (from differencing)
  data_clean <- na.omit(data)
  return(data_clean[,c(1,5)])
}


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


# ------------------------------------------------------------------------------
# download real data
# ------------------------------------------------------------------------------

# ---- Chemicals ---------------------------------------------------------------
tickers_chem <- c("SHW", "PPG", "EMN", "CE", "NEU", "SXT", "SCL")
corp_names_chem <- c("sherwin-williams",
                     "ppg-industries",
                     "eastman-chemical",
                     "celanese",
                     "newmarket", 
                     "sensient-technologies",
                     "stepan")

chemicals <- get_industry_growth_rates(tickers = tickers_chem, corp_names = corp_names_chem)
round(cor(chemicals[,-1]), digits = 4)

# ---- Building ----------------------------------------------------------------
tickers_building <- c("URI", "DHI", "VMC", "MLM", "ROL", "NVR", "EME")
corp_names_building <- c("united-rentals",
                         "dr-horton",
                         "vulcan-materials",
                         "martin-marietta-materials",
                         "rollins",
                         "nvr",
                         "emcor")

building <- get_industry_growth_rates(tickers = tickers_building, corp_names = corp_names_building)
round(cor(building[,-1]), digits = 4)

# ---- Retail - Supermarkets ---------------------------------------------------
tickers_retail <- c("WMT", "KR", "TGT", "DG", "VLGEA", "TJX", "M")
corp_names_retail <- c("walmart",
                       "kroger",
                       "target",
                       "dollar-general",
                       "village-super-market",
                       "tjx",
                       "macys")

retail <- get_industry_growth_rates(tickers = tickers_retail, corp_names = corp_names_retail)
round(cor(retail[,-1]), digits = 4)

# ---- Banks -------------------------------------------------------------------
tickers_banks <- c("JPM", "BAC", "WFC", "PNC", "C", "SCHW", "RJF")
corp_names_banks <- c("jpmorgan-chase",
                      "bank-of-america",
                      "wells-fargo",
                      "pnc-financial-services",
                      "citigroup",
                      "charles-schwab",
                      "raymond-james-financial")

banks <- get_industry_growth_rates(tickers = tickers_banks, corp_names = corp_names_banks)
round(cor(banks[,-1]), digits = 4)

# ---- Software ----------------------------------------------------------------
tickers_software <- c("MSFT", "SAP", "CDNS", "MSTR", "ANSS", "PTC", "MANH")
corp_names_software <- c("microsoft",
                         "sap-se",
                         "cadence-design-systems",
                         "microstrategy",
                         "ansys",
                         "ptc",
                         "manhattan-associates")

software <- get_industry_growth_rates(tickers = tickers_software, corp_names = corp_names_software)
round(cor(software[,-1]), digits = 4)

# ---- oil and gas -------------------------------------------------------------
tickers_oag <- c("XOM", "CVX", "SHEL", "BP", "WMB", "YPF", "ENB")
corp_names_oag <- c("exxon",
                    "chevron",
                    "shell",
                    "bp",
                    "williams",
                    "ypf-sociedad-anonima",
                    "enbridge-inc")

oag <- get_industry_growth_rates(tickers = tickers_oag, corp_names = corp_names_oag)
round(cor(oag[,-1]), digits = 4)


# ------------------------------------------------------------------------------
# store data
# ------------------------------------------------------------------------------

save(chemicals, file = "chemicals.Rda")
save(building, file = "building.Rda")
save(retail, file = "retail.Rda")
save(banks, file = "banks.Rda")
save(software, file = "software.Rda")
save(oag, file = "oag.Rda")

