# function to download real financial data from macrotrends.net


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
  
  out = as.data.frame(cleaned_data)
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






#' transform_to_stationary <- function(data) {
#'   #' Computes log differences for growth rates and adjusts quarterly input data 
#'   #' by trends and seasonality.
#'   #' 
#'   #' Parameters:
#'   #' - data (data.frame): Data.frame from download_revenue_data.
#'   #' 
#'   #' Returns:
#'   #' - Trend and seasonal adjusted logarithmic matrix of growth rates.
#' 
#'   # take the natural logarithm of the data
#'   log_data <- log(data$revenue)
#'   
#'   # de-trending and de-seasonalizing adjustment
#'   ts_log_data <- ts(log_data, frequency = 4)
#'   decomposed <- stl(ts_log_data, s.window = "periodic")
#'   seasonally_adjusted <- decomposed$time.series[, "remainder"]
#'   
#'   # calculate the first differences of the seasonally adjusted data
#'   diff_log_data <- diff(seasonally_adjusted)
#'   
#'   # convert result to a data.frame
#'   out <- data.frame(matrix(NA, dim(data)[1], 2))
#'   out[,1] <- data$date
#'   out[,2] <- matrix(c(NA, diff_log_data), ncol = 1)
#'   
#'   return(out)
#' }
