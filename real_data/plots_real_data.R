# functions which create plots of real quarterly revenue data analysis for different industries.


# ------------------------------------------------------------------------------
# import BVAR functions
# ------------------------------------------------------------------------------

function_files <- c(
  "_distr_samplers.R",
  "_monte_carlo_fun.R",
  "_forecast.R",
  "_hierarch_fun.R",
  "_data_prep.R",
  "_mcmc.R",
  "_ssvs.R",
  "_prior_specification.R",
  "bvar.R"
)
for (i in 1:length(function_files)) {
  source(file.path(dirname(getwd()), "functions", function_files[i]))
}
source(paste0(getwd(),"/plot_funs.r"))
rm(i)
rm(function_files)

# ------------------------------------------------------------------------------
# prepare real data 
# ------------------------------------------------------------------------------

# ---- import and transform to matrix ------------------------------------------
load("data/chemicals.Rda")
load("data/building.Rda")
load("data/retail.Rda")
load("data/software.Rda")
load("data/banks.Rda")
load("data/oag.Rda")
chem_mat <- as.matrix(chemicals[-59,-1])
build_mat <- as.matrix(building[-59,-1])
retail_mat <- as.matrix(retail[-59,-1])
soft_mat <- as.matrix(software[-59,-1])
bank_mat <- as.matrix(banks[-59,-1])
oag_mat <- as.matrix(oag[,-1])


# ------------------------------------------------------------------------------
# plotting growth rates
# ------------------------------------------------------------------------------

# ---- get correlations --------------------------------------------------------
round(cor(chemicals[,-1]), digits = 4)
round(cor(building[,-1]), digits = 4)
round(cor(retail[,-1]), digits = 4)
round(cor(banks[,-1]), digits = 4)
round(cor(software[,-1]), digits = 4)
round(cor(oag[,-1]), digits = 4)


# ---- store plots -------------------------------------------------------------
pdf(file = "chemicals.pdf", width = 8, height = 5.5)
par(mar = c(2.5, 2.5, 0, 0.5))
plot_gr(chemicals, legend_off = FALSE, adj_cex = 0.9, title = "off")
dev.off()

pdf(file = "building.pdf", width = 8, height = 5.5)
par(mar = c(2.5, 2.5, 0, 0.5))
plot_gr(building, legend_off = FALSE, adj_cex = 0.9, title = "off")
dev.off()

pdf(file = "retail.pdf", width = 8, height = 5.5)
par(mar = c(2.5, 2.5, 0, 0.5))
plot_gr(retail, legend_off = FALSE, adj_cex = 0.9, title = "off")
dev.off()

pdf(file = "software.pdf", width = 8, height = 5.5)
par(mar = c(2.5, 2.5, 0, 0.5))
plot_gr(software, legend_off = FALSE, adj_cex = 0.9, title = "off")
dev.off()

pdf(file = "banks.pdf", width = 8, height = 5.5)
par(mar = c(2.5, 2.5, 0, 0.5))
plot_gr(banks, legend_off = FALSE, y_lim = c(-90,90), adj_cex = 0.9, title = "off")
dev.off()

pdf(file = "oag.pdf", width = 8, height = 5.5)
par(mar = c(2.5, 2.5, 0, 0.5))
plot_gr(oag, legend_off = FALSE, y_lim = c(-90,90), adj_cex = 0.9, title = "off")
dev.off()


# ------------------------------------------------------------------------------
# plot predictive intervals
# ------------------------------------------------------------------------------
p <- 4
h <- 4
alpha <- 0.05

par(mar = c(2.2, 2.2, 0.1, 0.1))
plot_pred_intervals(soft_mat, software, p, h, alpha, "Software", c(-35,45), main = FALSE)
plot_pred_intervals(chem_mat, chemicals, p, h, alpha, "Chemicals", c(-50,40), main = FALSE)
plot_pred_intervals(build_mat, building, p, h, alpha, "Construction", c(-40,45), main = FALSE)
plot_pred_intervals(retail_mat, retail, p, h, alpha, "Retail", c(-50,60), main = FALSE)
plot_pred_intervals(bank_mat, banks, p, h, alpha, "Investment Banking", c(-60,75), main = FALSE)
plot_pred_intervals(oag_mat, oag, p, h, alpha, "Oil and Gas", c(-200,200), main = FALSE)


