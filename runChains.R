source('./MCMCfunc_combinedherds.R')

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
chain_no <- as.numeric(slurm_arrayid)
mcmc.size <- 1000000
burn_in <- 0
thinning <- 250
t_upd <- 3
simulated <- F
resultsFilepath <- "../Results/"

MCMCfunc(chain_no,mcmc.size,burn_in,thinning,t_upd,simulated,resultsFilepath)
