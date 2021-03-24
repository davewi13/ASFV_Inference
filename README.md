This is the code to carry out parameter inference on the African Swine Fever data, as described in the paper "Exact Bayesian inference of epidemiological parameters from mortality data: application to African swine fever virus".

The code is run from the runChains.R file.
This file calls either the MCMCfunc_combinedherds.R or the MCMCfunc_separateherds.R file depending on whether or not data is being pooled across herds to estimate global latent and infectious period parameters.
These functions require the C++ loglikelihood_rjmcmc_.cpp files (for calculating the log likelihood, updating proposed exposure/infection times and updating parameter values), the initFunc.R file (for proposing initial sets of exposure/infection times), the importHerds.R file (if loading in farm mortality data) and the simulate_data.R and simulate_SEIR.cpp files (if simulating data).

The code to produce Figures 3-6 in the paper is in the produce_figures.R file.

The file ASF_RJMCMC.yml contains the conda environment required to re-run these analyses.
