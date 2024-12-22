# function to download real financial data from macrotrends.net


download_revenue_data <- function(ticker, corp_name){
  #' Downloads quarterly real financial revenue data from listed equities from macrotrends.
  #' 
  #' Parameters:
  #' - ticker (string): Ticker of the equity to be downloaded in upper case. E.g. 'AAPL'.
  #' - corp_name (string): Corporate name in lower case of that company. E.g. 'apple'.
  #' 
  #' Returns:
  #' - a data.frame containing 3 columns:
  #'  - date: string of the format 'YYYY-MM-DD'.
  #'  - revenue: Numerical revenues in mio USD.
  #'  - TTM: Trailing 12 months, which shows the revenue of the past 12 months of the quarterly data.
  
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
  
  # computing TTM
  cleaned_data$TTM <- NA
  for (i in 4:nrow(cleaned_data)) {
    cleaned_data$TTM[i] <- sum(
      cleaned_data$revenue[(i-3):i]
      )
  }
  
  out = as.data.frame(cleaned_data)
  return(out)
}


get_growth_rates <- function(input){
  #' Compute the growth rates for quarterly revenue data
  #' 
  #' Parameters:
  #' - input (data.frame): The output of download_revenue_data.
  #' 
  #' Returns:
  #' - A data frame containing 5 columns with computed growth rates:
  #'  - date: Column indicating the date.
  #'  - quarterly_gr: Growth rate of raw quarterly data each quarter.
  #'  - twelve_months_gr: Growth rate of quarterly data compared to previous 12 month revenue.
  #'  - TTM_quarterly_gr: Growth rate of TTM data each quarter.
  #'  - TTM_twelve_months_gr: Growth rate of TTM data compared to previous 12 month TTM revenue.
  
  # initiate matrix
  growth_rates <- data.frame(
    matrix(NA, dim(input)[1], 5)
    )
  colnames(growth_rates) <- c(
    "date", 
    "quarterly_gr", 
    "twelve_months_gr", 
    "TTM_quarterly_gr", 
    "TTM_twelve_months_gr"
    )
  
  # date
  growth_rates[,1] <- input$date
  
  # quarterly data
  growth_rates[,2] <- c(NA, diff(input[,2]) / input[-length(input$TTM), 2])
  growth_rates[,3] <- c(
    rep(NA, 4), 
    (input[-(1:4), 2] - input[1:(nrow(input) - 4), 2])/input[1:(nrow(input) - 4), 2]
  )
  
  # TTM data
  growth_rates[,4] <- c(NA, diff(input[,3]) / input[-length(input$TTM), 3])
  growth_rates[,5] <- c(
    rep(NA, 4), 
    (input[-(1:4), 3] - input[1:(nrow(input) - 4), 3])/input[1:(nrow(input) - 4), 3]
  )
  
  # small break to prevent overuse of downloading capacity
  Sys.sleep(5)
  
  return(growth_rates)
}
