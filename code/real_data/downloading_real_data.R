# Load required libraries
library(rvest)
library(dplyr)

# Define the URL
url <- "https://macrotrends.net/stocks/charts/AAPL/apple/revenue"

# Read the webpage
webpage <- read_html(url)

# Extract the "Apple Quarterly Revenue" table
tables <- webpage %>%
  html_nodes(xpath = "//table[contains(@class, 'historical_data_table')]")

# Select the second table (Apple Quarterly Revenue)
quarterly_table <- tables[2] %>% html_table(fill = TRUE)

# Clean and format the data
cleaned_data <- table_data %>%
  rename(Date = 1, Revenue = 2) %>%
  mutate(
    Revenue = gsub("[$,]", "", Revenue), # Remove dollar sign and commas
    Revenue = as.numeric(Revenue),      # Convert to numeric
    Date = as.Date(Date, format = "%Y-%m-%d") # Convert to Date format
  )

# View the cleaned data
print(cleaned_data)

# Optionally save the data
write.csv(cleaned_data, "apple_quarterly_revenue.csv", row.names = FALSE)
