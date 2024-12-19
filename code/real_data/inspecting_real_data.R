# testing the behavior of real financial revenue growth for adjusting MC simulations.

# import sources
source(paste0(getwd(),"/downloading_real_data.r"))

# --- large cap pharmaceutical -------------------------------------------------
# download growth rates
eli_lilly <- get_growth_rates(download_revenue_data('LLY', 'eli-lilly'))
novo_nordisk <- get_growth_rates(download_revenue_data('NVO', 'novo-nordisk'))
johnson_johnson <- get_growth_rates(download_revenue_data('JNJ', 'johnson_johnson'))
abbvie <- get_growth_rates(download_revenue_data('ABBV', 'abbvie'))
merck <- get_growth_rates(download_revenue_data('MRK', 'merck'))
astrazeneca <- get_growth_rates(download_revenue_data('AZN', 'astrazeneca'))
novartis <- get_growth_rates(download_revenue_data('NVS', 'novartis-ag'))
pfizer <- get_growth_rates(download_revenue_data('PFE', 'pfizer'))

# merge growth rates
merged_TTM_quarterly_gr = data.frame(
  year = eli_lilly$date,
  eli_lilly = eli_lilly$TTM_quarterly_gr,
  novo_nordisk = novo_nordisk$TTM_quarterly_gr,
  johnson_johnson = johnson_johnson$TTM_quarterly_gr,
  merck = merck$TTM_quarterly_gr,
  astrazeneca = astrazeneca$TTM_quarterly_gr,
  novartis = novartis$TTM_quarterly_gr,
  pfizer = pfizer$TTM_quarterly_gr
)

# plot growth rates
growth_rates <- merged_TTM_quarterly_gr[, -1]
matplot(
  merged_TTM_quarterly_gr$year, growth_rates, type = "l", lty = 1, 
  xlab = "Year", ylab = "Growth Rate", 
  main = "Growth Rates of Companies Over Time"
)
legend(
  "topright", legend = colnames(growth_rates), 
  col = 1:ncol(growth_rates), lty = 1, cex = 0.8
)
