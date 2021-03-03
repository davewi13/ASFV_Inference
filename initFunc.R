initFunc <- function(remtimes,beta_ext,beta_int,exp_mean,exp_shape,inf_mean,inf_shape,M,Tfinal){
  # This function creates an initial data frame with estimates of the exposure and infection times
  # Inputs:
  # remtimes - a vector of removal (death) times for a given herd
  # beta_ext, beta_int - starting parameter values for the external and internal (within-herd) transmission rates
  # exp_mean, exp_shape - mean and shape of the exposed (latent) period
  # inf_mean, inf_shape - mean and shape of the infectious period
  # M - Herd size
  # Tfinal - time of cull
  # Outputs:
  # completedata - a data frame with proposed exposure and infection times for each removal
  
  # calculate the rate parameters of the gamma distributed latent and infectious periods
  exp_rate <- exp_shape/exp_mean
  inf_rate <- inf_shape/inf_mean
  
  # If there is no external transmission rate (set to zero) then it is possible to propose sets of times
  # where multiple individuals are exposed before the first infection, for which the likelihood is NA
  # Don't allow such sets of times
  viable <- FALSE
  iter <- 0
  while(viable == FALSE){
    # Determine the number of infections and create a data frame
    n            = length(remtimes)
    completedata = matrix(0,nrow=n,ncol=3)
    completedata = as.data.frame(completedata)
    colnames(completedata) = c("Exposure","Infection","Removal")
    
    # Add proposed exposure and infection times by sampling from gamma distribution with initial parameter estimates
    completedata$Removal <- remtimes
    completedata$Infection <- remtimes-rgamma(n,inf_shape,inf_rate)
    completedata$Exposure <- completedata$Infection-rgamma(n,exp_shape,exp_rate)
    completedata <- completedata[order(completedata$Exposure),]
    
    # If we have tried this lots of times and we don't have a viable set of times (one with a likelihood)
    # then propose exposure times in a systematic way to ensure something for which we can calculate a likelihood
    if((iter > 1000) & (beta_ext == 0)){
      completedata <- completedata[order(completedata$Exposure),]
      for(i in 2:n){
        if(completedata$Exposure[i] < min(completedata$Infection)){
          completedata$Exposure[i] <- (min(completedata$Infection)+completedata$Infection[i])/2
        }
      }
    }

    # Check that there is a likelihood and leave the loop if so
    if((loglikelihood(completedata,beta_ext,beta_int,exp_mean,exp_shape,inf_mean,inf_shape,M,Tfinal) != -Inf)){(viable <- TRUE)}
    iter <- iter+1
  } 
  # Order the set of proposed times by the removal time and output
  completedata <- completedata[order(completedata$Removal),]
  
  return(completedata)
  
}
