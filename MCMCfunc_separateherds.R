MCMCfunc <- function(chain_no,mcmc.size,burn_in,thinning,t_upd,simulated,resultsFilepath){
  # Function to run RJ MCMC chains fitting latent and infectious parameters for each herd
  # Code is very similar to the code for combined herds, so only comment where there are differences  
  
  set.seed(chain_no)
  source('./initFunc.R')
  library(Rcpp)
  sourceCpp("./loglikelihood_rjmcmc_separateherds_twobetas_freq_dep.cpp")
  if(simulated){
    M1 <- M2 <- M3 <- M4 <- M5 <- M6 <-  M7 <- M8 <- M9 <- 500
    Tfinal1 <- Tfinal2 <- Tfinal3 <- Tfinal4 <- Tfinal5 <- Tfinal6 <-  Tfinal7 <- Tfinal8 <- Tfinal9 <- 28
    sourceCpp("./simulate_SEIR.cpp")
    source('./simulate_data.R')
    list2env(simulate_data(N=M1, Tfinal=Tfinal1),globalenv())
  } else {
    source('./importHerds.R')
  }
  
  # Reserve space for the MCMC output (more space here as more parameters)
  mcmc.samples           <- matrix(0,(mcmc.size-burn_in)/thinning+1,60) 
  mcmc.samples           <- as.data.frame(mcmc.samples)
  colnames(mcmc.samples) <- c("beta1_ext","beta2_ext","beta3_ext","beta4_ext","beta5_ext","beta6_ext","beta7_ext","beta8_ext","beta9_ext","beta1_int","beta2_int","beta3_int","beta4_int","beta5_int","beta6_int","beta7_int","beta8_int","beta9_int","exp_mean1","exp_mean2","exp_mean3","exp_mean4","exp_mean5","exp_mean6","exp_mean7","exp_mean8","exp_mean9","exp_shape1","exp_shape2","exp_shape3","exp_shape4","exp_shape5","exp_shape6","exp_shape7","exp_shape8","exp_shape9","inf_mean1","inf_mean2","inf_mean3","inf_mean4","inf_mean5","inf_mean6","inf_mean7","inf_mean8","inf_mean9","inf_shape1","inf_shape2","inf_shape3","inf_shape4","inf_shape5","inf_shape6","inf_shape7","inf_shape8","inf_shape9","accept_exp_mean","accept_exp_shape","accept_inf_mean","accept_inf_shape","accept_beta_int","accept_beta_ext")
  
  # Now have vectors for latent and infectious mean parameters instead of single values
  beta_ext_vec <- rep(0.001,9)
  beta_int_vec <- rep(1,9)
  beta_int_shape <- 2
  beta_ext_shape <- 1
  beta_int_rate <- 2/2
  beta_ext_rate <- 0.001
  exp_mean_vec <- rgamma(9,10,10/6.25)
  exp_mean_rate <- 10/6.25
  exp_mean_shape <- 10
  exp_shape_vec <- rgamma(9,5,5/19.39)
  exp_shape_rate <- 5/19.39
  exp_shape_shape <- 5
  inf_mean_vec <- rgamma(9,10,10/9.12)
  inf_mean_rate <- 10/9.12
  inf_mean_shape <- 10
  inf_shape_vec <- rgamma(9,5,5/22) #22
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
  init_vals <- list(exp_mean_vec = exp_mean_vec, exp_shape_vec = exp_shape_vec, inf_mean_vec = inf_mean_vec, inf_shape_vec = inf_shape_vec)
  
  M.vec <- c(M1,M2,M3,M4,M5,M6,M7,M8,M9)
  Tfinal.vec <- c(Tfinal1,Tfinal2,Tfinal3,Tfinal4,Tfinal5,Tfinal6,Tfinal7,Tfinal8,Tfinal9)

  remtimes.ls <- list(remtimes1=remtimes1,
                      remtimes2=remtimes2,
                      remtimes3=remtimes3,
                      remtimes4=remtimes4,
                      remtimes5=remtimes5,
                      remtimes6=remtimes6,
                      remtimes7=remtimes7,
                      remtimes8=remtimes8,
                      remtimes9=remtimes9)
  
  for(i in 1:9){
    assign(paste0("completedata", i), initFunc(remtimes.ls[[i]],beta_ext_vec[i],beta_int_vec[i],exp_mean_vec[i],exp_shape_vec[i],inf_mean_vec[i],inf_shape_vec[i],M.vec[i],Tfinal.vec[i]))
  }
  
  completedata.ls <- list(completedata1=completedata1,
                          completedata2=completedata2,
                          completedata3=completedata3,
                          completedata4=completedata4,
                          completedata5=completedata5,
                          completedata6=completedata6,
                          completedata7=completedata7,
                          completedata8=completedata8,
                          completedata9=completedata9)
  
  # Store the initial values - again more values to store
  mcmc.samples[1,1:9] <- beta_ext_vec
  mcmc.samples[1,10:18] <- beta_int_vec
  mcmc.samples[1,19:27] <- exp_mean_vec
  mcmc.samples[1,28:36] <- exp_shape_vec
  mcmc.samples[1,37:45] <- inf_mean_vec
  mcmc.samples[1,46:54] <- inf_shape_vec
  mcmc.samples[1,]$accept_exp_mean <- accept_exp_mean
  mcmc.samples[1,]$accept_exp_shape <- accept_exp_shape
  mcmc.samples[1,]$accept_inf_mean <- accept_inf_mean
  mcmc.samples[1,]$accept_inf_shape <- accept_inf_shape
  mcmc.samples[1,]$accept_beta_int <- accept_beta_int
  mcmc.samples[1,]$accept_beta_ext <- accept_beta_ext


  for (i in 2:mcmc.size){
    
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
    # Separate parameter values for each herd now
    for(j in 1:9){
      output_exp_mean <- update_exp_mean(completedata.ls[[j]],exp_mean_vec[j],exp_mean_sigma,exp_shape_vec[j],exp_mean_rate,exp_mean_shape,Tfinal.vec[j],accept_exp_mean)
    exp_mean_vec[j] <- output_exp_mean[[1]]
    accept_exp_mean <- output_exp_mean[[2]]
    }
    for(j in 1:9){
      output_exp_shape <- update_exp_shape(completedata.ls[[j]],exp_mean_vec[j],exp_shape_vec[j],exp_shape_sigma,exp_shape_rate,exp_shape_shape,Tfinal.vec[j],accept_exp_shape)
    exp_shape_vec[j] <- output_exp_shape[[1]]
    accept_exp_shape <- output_exp_shape[[2]]
    }
    for(j in 1:9){
      output_inf_mean <- update_inf_mean(completedata.ls[[j]],inf_mean_vec[j],inf_mean_sigma,inf_shape_vec[j],inf_mean_rate,inf_mean_shape,Tfinal.vec[j],accept_inf_mean)
    inf_mean_vec[j] <- output_inf_mean[[1]]
    accept_inf_mean <- output_inf_mean[[2]]
    }
    for(j in 1:9){
      output_inf_shape <- update_inf_shape(completedata.ls[[j]],inf_mean_vec[j],inf_shape_vec[j],inf_shape_sigma,inf_shape_rate,inf_shape_shape,Tfinal.vec[j],accept_inf_shape)
    inf_shape_vec[j] <- output_inf_shape[[1]]
    accept_inf_shape <- output_inf_shape[[2]]
    }

    for(j in 1:9){
      for(k in 1:t_upd){
        completedata.ls[[j]] <- update_times(completedata.ls[[j]],beta_ext_vec[j],beta_int_vec[j],exp_mean_vec[j],exp_shape_vec[j],inf_mean_vec[j],inf_shape_vec[j],M.vec[j],Tfinal.vec[j])
      }
    }
    
    # Store the current iterates for parameters and acceptance rates - more to store
    if((i >= burn_in) & ((i %% thinning) == 0)){
      update_row <- (i%/%thinning)-(burn_in/thinning-1)
      mcmc.samples[update_row,1:9] <- beta_ext_vec
      mcmc.samples[update_row,10:18] <- beta_int_vec
      mcmc.samples[update_row,19:27] <- exp_mean_vec
      mcmc.samples[update_row,28:36] <- exp_shape_vec
      mcmc.samples[update_row,37:45] <- inf_mean_vec
      mcmc.samples[update_row,46:54] <- inf_shape_vec
      mcmc.samples[update_row,]$accept_exp_mean <- accept_exp_mean/(9*i)
      mcmc.samples[update_row,]$accept_exp_shape <- accept_exp_shape/(9*i)
      mcmc.samples[update_row,]$accept_inf_mean <- accept_inf_mean/(9*i)
      mcmc.samples[update_row,]$accept_inf_shape <- accept_inf_shape/(9*i)
      mcmc.samples[update_row,]$accept_beta_ext <- accept_beta_ext/(9*i)
      mcmc.samples[update_row,]$accept_beta_int <- accept_beta_int/(9*i)
    }
    if((i %% 10000) == 0){
      print(paste(i, "completed iterations of run no", chain_no))
      results <- list(samps = mcmc.samples, times = completedata.ls, init_vals = init_vals, chain_length = mcmc.size, burn_in = burn_in, thinning = thinning, priors = priors_list)
      save(results,file=paste(resultsFilepath,'chain',chain_no,'.Rdata',sep=''))
    }
  }
  
  results <- list(samps = mcmc.samples, times = completedata.ls, init_vals = init_vals, chain_length = mcmc.size, burn_in = burn_in, thinning = thinning, priors = priors_list)
  save(results,file=paste(resultsFilepath,'chain',chain_no,'.Rdata',sep=''))
  # Return the MCMC sample of the two model parameters
  return(results)
  
}
