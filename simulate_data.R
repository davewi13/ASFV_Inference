simulate_data <- function(N, Tfinal){
# Function to simulate data from SEIR model (C++ function does most of the work, this just calls it)
# Inputs:
# N - number of animals in herd
# Tfinal - time of cull
  
  # Set up parameter values for simulation
  set.seed(10)
  beta_int_vec <- c(0.4,0.3,0.25,0.25,0.35,0.45,0.3,0.35,0.4)
  beta_ext_vec <- c(0.0002,0.0002,0.0001,0.0002,0.0001,0.0002,0.0001,0.0002,0.0001)
  exp_shape <- 19.39
  exp_rate <- 19.39/6.25
  inf_shape <- 22
  inf_rate <- 22/9.12
  # Loop through all herds
  for(i in 1:9){
    # Only accept simulated datasets with 15 or more removals (due to stochasticity sometimes the epidemic dies immediately)
    invalid <- T
    while(invalid){
      results <- SEIR_IBM(N=N,I0=1,exp_shape=exp_shape,exp_rate=exp_rate,inf_shape=inf_shape,inf_rate=inf_rate,beta_int=beta_int_vec[i],beta_ext=beta_ext_vec[i],t0=0,tfinal=Tfinal)
      if(tail(results[[5]],1) >= 15){
        invalid <- F
      }
    }
    # Determine removal times, round to nearest day and output in a vector for each herd
    remtimes <- numeric(0)
    counter <- 0
    removals_vec <- results[[5]]
    times_vec <- results[[1]]
    for(j in 2:length(results[[5]])){
      if(removals_vec[j] > removals_vec[j-1]){
        remtimes[counter] <- times_vec[j]
        counter <- counter + 1
      }
    }
    assign(paste0("remtimes",i), round(remtimes))
  }
  
  # Output list of removal times for each herd
  return(list(remtimes1=remtimes1,
              remtimes2=remtimes2,
              remtimes3=remtimes3,
              remtimes4=remtimes4,
              remtimes5=remtimes5,
              remtimes6=remtimes6,
              remtimes7=remtimes7,
              remtimes8=remtimes8,
              remtimes9=remtimes9))
}
