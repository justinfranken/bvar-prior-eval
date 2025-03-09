# Evaluating Prior Distributions in Bayesian Vector Autoregressions: A Comparative Analysis

This project has been developed for my master's thesis *"Evaluating Prior Distributions in Bayesian Vector Autoregressions: A Comparative Analysis"* at the University of Bonn. The goal of this work is to compare the predictive power of different *Minnesota* style priors in Bayesian vector autoregressive models, especially in settings with high volatily time series and limited amount of observations at a time. Comparative analysis is done on simulated data and real world quarterly revenue data sourced from [macrotrends](https://www.macrotrends.net/). Below you can find an overview of the repositories structure and the accommodating content.



## Overview
```
bvar-prior-eval
│   LICENSE
│   README.md
│   lib.R
│   mcmc_convergence_tests.R
│   simulation.R
│   temp_execute.R
│
└───functions
│   │   _data_prep.R
│   │   _distr_sampler.R
│   │   _forecast.R
│   │   _hierarch_fun.R
│   │   _mcmc.R
│   │   _mcmc_convergence.R
│   │   _monte_carlo_fun.R
│   │   _prior_specification.R
│   │   _ssvs.R
│   │   bvar.R
│   │   simulate_data.R
│
└───plots
│   (20 .pdf and .png figures)
│
└───real_data
│   │   application_evaluations.txt
│   │   downloading_real_data.R
│   │   plots_funs.R
│   │   plots_real_data.R
│   │   real_data_application.R
│   │   real_data_funs.R
│   │
│   └───data
│       │   banks.Rda
│       │   building.Rda
│       │   chemicals.Rda
│       │   oag.Rda
│       │   retail.Rda
│       │   software.Rda
│   
└───simulations_results
    │   simulation_results.txt
    │   (12 .RData files storing the simulation
    │   results for each simulation setting)
```

## Notable Directories and Files

* `/lib.R`: This script loads all required R packages necessary to run the BVAR algorithms and do predictions with evaluations as well to download and transform real world quarterly revenue data.
* `/simulation.R`: Contains the Monte Carlo simulation process. It reads all BVAR functions in from /functions.
* `/functions/`: Contains all functions used to perform BVAR analysis, evaluations and perform monte carlo simulations with it. Key scripts in this file are:
    * `_forecast.R`: Includes functions to predict posterior coefficient draws and meassure RMSFE and predictive interval accuracy.
    * `_hierarch_fun.R`: Includes functions for establishing hierarchical priors, following Giannone et al. (2015) and the R-package [bvar](https://github.com/nk027/bvar).
    * `_ssvs.R`: Includes functions needed for SSVS bvar analysis, following Koop (2013), where rank-1 Sherman-Morrison updates are used to improve computational efficiency.
    * `_mcmc.R`: Functions to sample from posterior predictive distributions utilizing Makrov Chain Monte Carlo Simulations.
    * `_prior_specification.R`: Includes functions which establish Minnesota prior beliefs.
    * `simulate_data.R`: Includes functions for our data generating process, needed for the monte carlo simulation.
* `/real_data/`: Containts functions and corresponding analysis results/plots of real world data applications. Key scripts and files in this directory are:
    * `application_evaluations.txt`: File containing the rolling window BVAR analysis of real world quarterly revenue data for 6 different industries.
    * `downloading_real_data.R`: Functions which source and transform real world quarterly revenue data from [macrotrends](https://www.macrotrends.net/) to stationary time series.
    * `real_data_application.R`: Functions performing our real data application evaluation of 6 different industries.
* `/simulation_results/`: Contains the results of the monte carlo simulation for all simulation settings (small, medium, large data generating processes / high and low correlation / with and without shocks).



## Getting started

You can either download the repository directly via Github or clone the repository to you system by running the following in your console:

```
git clone https://github.com/justinfranken/bvar-prior-eval.git
```

Every R file which requires functions from other R files has in each preliminary part a source function. The R file `simulation.R` for example contains:

```R
function_files <- c(
  "_distr_samplers.R",
  "_monte_carlo_fun.R",
  "_forecast.R",
  "_hierarch_fun.R",
  "_data_prep.R",
  "_mcmc.R",
  "_ssvs.R",
  "_prior_specification.R",
  "bvar.R",
  "simulate_data.R"
)
for (i in 1:length(function_files)) {
  source(paste0(getwd(),"/functions/", function_files[i]))
}
```

## References
Main papers and books used for this R project are:
* Giannone, D., M. Lenza, and G. Primiceri (2015). “Prior Selection for Vector Autoregressions”. In: The Review of Economics and Statistics 97.2, pp. 436–451. url: https://direct.mit.edu/rest/article-abstract/97/2/436/58236/Prior-Selection-for-Vector-Autoregressions?redirectedFrom=fulltext.
* Karlsson, S. (2013). *“Forecasting with Bayesian vector autoregressions”*. In: Handbook of Economic Forecasting. Ed. by G. Elliott and A. Timmermann. Vol. 2. Amsterdam: Elsevier, pp. 791–897.
* Kuschnig, N. and L. Vashold (2021). “BVAR: Bayesian Vector Autoregressions with Hierarchical Prior Selection in R”. In: Journal of Statistical Software 100.14, pp. 1–27. doi: 10.18637/jss.v100.i14. url: https://doi.org/10.18637/jss.v100.i14.
* Koop, G. 2013). “Forecasting with Medium and Large Bayesian VARs”. In: Journal of Applied Econometrics 28, pp. 177–203. doi: 10.1002/jae.1270. url: https://doi.org/10.
* Rockova, Veronika (Nov. 2013). “Bayesian Variable Selection in High-dimensional Applications”. PhD thesis. url: http://hdl.handle.net/1765/51587.

For the hierarchical prior estimation, inspirations are taken from the R-package [bvar](https://github.com/nk027/bvar).