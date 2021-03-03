MCMCfunc <- function(chain_no,mcmc.size,burn_in,thinning,t_upd,simulated,resultsFilepath){
# Function to run RJ MCMC chains
# Inputs: chain_no - the chain number to set the seed and ensure different (and reproducible) chains
# mcmc.size - length of the chain
# burn_in - length of any burn-in period
# thinning - number of observations to thin over
# t_upd - the number of individuals to propose an updated set of times for on each iteration
# simulated - true or false stating whether to simulate data or import information about given herds
# resultsFilepath - where to store results
# Outputs: samps - sampled parameter values
# times - final set of proposed times
# init_vals - values used to initialise the chain
# chain_length - number of iterations
# burn_in - number of iterations discarded as burn in
# thinning - value of any thinning parameter
# priors - information about the priors used
    
  set.seed(chain_no)
  # Load function to propose an initial set of times
  source('./initFunc.R')
  # Load Rcpp library
  library(Rcpp)
  # Load c++ functions to calculate log likelihoods, update parameter values and update times
  sourceCpp("./loglikelihood_rjmcmc_combinedherds.cpp")
  # Simulating values or using data from studied herds
  if(simulated){
    # Set herd sizes and cull times
    M1 <- M2 <- M3 <- M4 <- M5 <- M6 <-  M7 <- M8 <- M9 <- 500
    Tfinal1 <- Tfinal2 <- Tfinal3 <- Tfinal4 <- Tfinal5 <- Tfinal6 <-  Tfinal7 <- Tfinal8 <- Tfinal9 <- 28
    # Load R and C++ code to simulate from SEIR model
    sourceCpp("./simulate_SEIR.cpp")
    source('./simulate_data.R')
    # simulate data and load into global environment
    list2env(simulate_data(N=M1, Tfinal=Tfinal1),globalenv())
  } else {
    # If using real data, just load it in
    source('./importHerds.R')
  }
  
  # Reserve space for the MCMC output
  mcmc.samples           <- matrix(0,(mcmc.size-burn_in)/thinning+1,28) 
  mcmc.samples           <- as.data.frame(mcmc.samples)
  colnames(mcmc.samples) <- c("beta1_ext","beta2_ext","beta3_ext","beta4_ext","beta5_ext","beta6_ext","beta7_ext","beta8_ext","beta9_ext","beta1_int","beta2_int","beta3_int","beta4_int","beta5_int","beta6_int","beta7_int","beta8_int","beta9_int","exp_mean","exp_shape","inf_mean","inf_shape","accept_exp_mean","accept_exp_shape","accept_inf_mean","accept_inf_shape","accept_beta_int","accept_beta_ext")
  
  # Initialize the model unknowns
  beta_ext_vec <- rep(0,9)
  beta_int_vec <- rep(1,9)
  beta_int_shape <- 2
  beta_ext_shape <- 1
  beta_int_rate <- 2/2
  beta_ext_rate <- 0.001
  exp_mean <- rgamma(1,10,10/6.25)
  exp_mean_rate <- 10/6.25
  exp_mean_shape <- 10
  exp_shape <- rgamma(1,5,5/19.39)
  exp_shape_rate <- 5/19.39
  exp_shape_shape <- 5
  inf_mean <- rgamma(1,10,10/9.12)
  inf_mean_rate <- 10/9.12
  inf_mean_shape <- 10
  inf_shape <- rgamma(1,5,5/22)
  inf_shape_rate <- 5/22
  inf_shape_shape <- 5
  priors_list <- list(beta_int_rate = beta_int_rate,
                      beta_ext_rate = beta_ext_rate,
                      beta_int_shape = beta_int_shape,
                      beta_ext_shape = beta_ext_shape,
                      exp_mean_rate = exp_mean_rate, 
                      exp_mean_shape = exp_mean_shape,
                      exp_shape_rate = exp_shape_rate,
                      exp_shape_shape = exp_shape_shape,
                      inf_mean_rate = inf_mean_rate,
                      inf_mean_shape = inf_mean_shape,
                      inf_shape_rate = inf_shape_rate,
                      inf_shape_shape = inf_shape_shape)
  # Set proposal variances and initialise acceptance rates
  beta_int_sigma <- 0.5
  beta_ext_sigma <- 0.005
  exp_mean_sigma <- 0.5
  exp_shape_sigma <- 2
  inf_mean_sigma <- 1
  inf_shape_sigma <- 10
  accept_exp_mean <- 0
  accept_exp_shape <- 0
  accept_inf_mean <- 0
  accept_inf_shape <- 0
  accept_beta_int <- 0
  accept_beta_ext <- 0
  accept_times <- 0
  init_vals <- list(exp_mean = exp_mean, exp_shape=exp_shape, inf_mean = inf_mean, inf_shape = inf_shape)
  
  # Store herd sizes and cull times in vectors
  M.vec <- c(M1,M2,M3,M4,M5,M6,M7,M8,M9)
  Tfinal.vec <- c(Tfinal1,Tfinal2,Tfinal3,Tfinal4,Tfinal5,Tfinal6,Tfinal7,Tfinal8,Tfinal9)

  # Put removal times for each herd in a list
  remtimes.ls <- list(remtimes1=remtimes1,
                      remtimes2=remtimes2,
                      remtimes3=remtimes3,
                      remtimes4=remtimes4,
                      remtimes5=remtimes5,
                      remtimes6=remtimes6,
                      remtimes7=remtimes7,
                      remtimes8=remtimes8,
                      remtimes9=remtimes9)
  
  # Propose intial sets of times for each herd
  for(i in 1:9){
    assign(paste0("completedata", i), initFunc(remtimes.ls[[i]],beta_ext_vec[i],beta_int_vec[i],exp_mean,exp_shape,inf_mean,inf_shape,M.vec[i],Tfinal.vec[i]))
  }
  
  # Store data frames of proposed times in list
  completedata.ls <- list(completedata1=completedata1,
                          completedata2=completedata2,
                          completedata3=completedata3,
                          completedata4=completedata4,
                          completedata5=completedata5,
                          completedata6=completedata6,
                          completedata7=completedata7,
                          completedata8=completedata8,
                          completedata9=completedata9)
  
  # Store the initial values
  mcmc.samples[1,1:9] <- beta_ext_vec
  mcmc.samples[1,10:18] <- beta_int_vec
  mcmc.samples[1,]$exp_mean <- exp_mean
  mcmc.samples[1,]$exp_shape <- exp_shape
  mcmc.samples[1,]$inf_mean <- inf_mean
  mcmc.samples[1,]$inf_shape <- inf_shape
  mcmc.samples[1,]$accept_exp_mean <- accept_exp_mean
  mcmc.samples[1,]$accept_exp_shape <- accept_exp_shape
  mcmc.samples[1,]$accept_inf_mean <- accept_inf_mean
  mcmc.samples[1,]$accept_inf_shape <- accept_inf_shape
  mcmc.samples[1,]$accept_beta_int <- accept_beta_int
  mcmc.samples[1,]$accept_beta_ext <- accept_beta_ext


  # The MCMC iterations
  for (i in 2:mcmc.size){
    
    # Update the parameters and acceptance rates
    # If external transmission is set to zero then best to comment the first for loop out
    for(j in 1:9){
      beta_ext_output <- update_beta_ext(completedata.ls[[j]],beta_ext_vec[j],beta_ext_sigma,beta_int_vec[j],beta_ext_rate,beta_ext_shape, M.vec[j],Tfinal.vec[j], accept_beta_ext)
      beta_ext_vec[j] <- beta_ext_output[[1]]
      accept_beta_ext <- beta_ext_output[[2]]
    }
    for(j in 1:9){
       beta_int_output <- update_beta_int(completedata.ls[[j]],beta_int_vec[j],beta_int_sigma,beta_ext_vec[j],beta_int_rate,beta_int_shape, M.vec[j],Tfinal.vec[j], accept_beta_int)
       beta_int_vec[j] <- beta_int_output[[1]]
       accept_beta_int <- beta_int_output[[2]]
    }
    output_exp_mean <- update_exp_mean(completedata.ls,exp_mean,exp_mean_sigma,exp_shape,exp_mean_rate,exp_mean_shape,Tfinal.vec,accept_exp_mean)
    exp_mean <- output_exp_mean[[1]]
    accept_exp_mean <- output_exp_mean[[2]]
    output_exp_shape <- update_exp_shape(completedata.ls,exp_mean,exp_shape,exp_shape_sigma,exp_shape_rate,exp_shape_shape,Tfinal.vec,accept_exp_shape)
    exp_shape <- output_exp_shape[[1]]
    accept_exp_shape <- output_exp_shape[[2]]
    output_inf_mean <- update_inf_mean(completedata.ls,inf_mean,inf_mean_sigma,inf_shape,inf_mean_rate,inf_mean_shape,Tfinal.vec,accept_inf_mean)
    inf_mean <- output_inf_mean[[1]]
    accept_inf_mean <- output_inf_mean[[2]]
    output_inf_shape <- update_inf_shape(completedata.ls,inf_mean,inf_shape,inf_shape_sigma,inf_shape_rate,inf_shape_shape,Tfinal.vec,accept_inf_shape)
    inf_shape <- output_inf_shape[[1]]
    accept_inf_shape <- output_inf_shape[[2]]

    # Update t_upd of the exposure and infection times for each herd
    for(j in 1:9){
      for(k in 1:t_upd){
        completedata.ls[[j]] <- update_times(completedata.ls[[j]],beta_ext_vec[j],beta_int_vec[j],exp_mean,exp_shape,inf_mean,inf_shape,M.vec[j],Tfinal.vec[j])
      }
    }
    
    # Store the current iterates for parameters and acceptance rates
    if((i >= burn_in) & ((i %% thinning) == 0)){
      update_row <- (i%/%thinning)-(burn_in/thinning-1)
      mcmc.samples[update_row,1:9] <- beta_ext_vec
      mcmc.samples[update_row,10:18] <- beta_int_vec
      mcmc.samples[update_row,]$exp_mean <- exp_mean
      mcmc.samples[update_row,]$exp_shape <- exp_shape
      mcmc.samples[update_row,]$inf_mean <- inf_mean
      mcmc.samples[update_row,]$inf_shape <- inf_shape
      mcmc.samples[update_row,]$accept_exp_mean <- accept_exp_mean/i
      mcmc.samples[update_row,]$accept_exp_shape <- accept_exp_shape/i
      mcmc.samples[update_row,]$accept_inf_mean <- accept_inf_mean/i
      mcmc.samples[update_row,]$accept_inf_shape <- accept_inf_shape/i
      mcmc.samples[update_row,]$accept_beta_ext <- accept_beta_ext/(9*i)
      mcmc.samples[update_row,]$accept_beta_int <- accept_beta_int/(9*i)
    }
    # Output message and results to keep track through several iterations
    if((i %% 10000) == 0){
      print(paste(i, "completed iterations of run no", chain_no))
      results <- list(samps = mcmc.samples, times = completedata.ls, init_vals = init_vals, chain_length = mcmc.size, burn_in = burn_in, thinning = thinning, priors = priors_list)
      save(results,file=paste(resultsFilepath,'chain',chain_no,'.Rdata',sep=''))
    }
  }
  
  # Output
  results <- list(samps = mcmc.samples, times = completedata.ls, init_vals = init_vals, chain_length = mcmc.size, burn_in = burn_in, thinning = thinning, priors = priors_list)
  save(results,file=paste(resultsFilepath,'chain',chain_no,'.Rdata',sep=''))
  # Return the MCMC sample of the two model parameters
  return(results)
  
}