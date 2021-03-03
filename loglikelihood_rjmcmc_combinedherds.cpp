#include <Rcpp.h>
#include <cstdlib>
#include <cmath> 
#include <chrono>
using namespace Rcpp;

double inf = std::numeric_limits<double>::infinity();
double minusinf = -inf;

// [[Rcpp::export]]
NumericVector Rcpp_sort(NumericVector x, NumericVector y) {
    // Order the elements of x by sorting y
    // First create a vector of indices
    IntegerVector idx = seq_along(x) - 1;
    // Then sort that vector by the values of y
    std::sort(idx.begin(), idx.end(), [&](int i, int j){return y[i] < y[j];});
    // And return x in that order
    return x[idx];
}

// [[Rcpp::export]]
double totaltime_in_stage(const NumericVector &v1, const NumericVector &v2){

// Calculate the total time spent in a stage based on the vectors of times in and times out.
// Inputs: v1, v2 - times of entry to stage and times of exit from stage
// Outputs: the sum of the times spent in the stage

  NumericVector times_vec = v2-v1;
  for (int k = 0; k < times_vec.size(); ++k){
	  if(v2[k] == inf){
		  times_vec[k] = 0;
	  }
  }
  return sum(times_vec);
  
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
double totaltime_in_stage_log(const NumericVector &v1, const NumericVector &v2){
  
// Calculate the log of the total time spent in a stage based on the times in and times out.
// Inputs: v1, v2 - times of entry to stage and times of exit from stage
// Outputs: the sum of the log of the total time spent in the stage
  
  NumericVector times_vec = v2-v1;
  NumericVector log_times_vec = log(times_vec);
  for (int k = 0; k < log_times_vec.size(); ++k){
	  if(v2[k] == inf){
		  log_times_vec[k] = 0;
	  }
  }
  return sum(log_times_vec);
  
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
double totaltime_infpressure(Rcpp::DataFrame &eventtimes, double M, double Tfinal){

// Function to calculate the total time of infection pressure applied by infectious individuals on susceptibles.
// Inputs: eventtimes - a data frame with all proposed exposure, infection and removal times
// M - the total size of the herd
// Tfinal - the time of cull
// Outputs: the total time of infectious pressure

  // Extract the exposure times, infection times and removal times
  NumericVector exptimes = eventtimes["Exposure"];
  NumericVector inftimes = eventtimes["Infection"];
  NumericVector remtimes = eventtimes["Removal"];
  
  // Make new vectors without inf values.  Inf values are used for cryptic infections/exposures which will not have, for example, a removal time
  NumericVector inftimes_noinf = inftimes[inftimes != inf];
  NumericVector remtimes_noinf = remtimes[remtimes != inf];
  
  // Calculate the number of each type of event (exposure, infection, removal)
  int Nexp = exptimes.size();
  int Ninf = inftimes_noinf.size();
  int Nrem = remtimes_noinf.size();
  int no_events = Nexp+Ninf+Nrem;

  // Create some vectors for storage
  NumericVector alltimes (no_events);
  NumericVector allprocesses (no_events);
  
  // Create a vector with all the event times in it (alltimes) and an index vector (allprocesses) stating what type of event each is
  for (int i = 0; i < allprocesses.size(); ++i){
	if (i < Nexp) {
		alltimes[i] = exptimes[i];
	} else if ((i >= Nexp) & (i < (Nexp+Ninf))){
		alltimes[i] = inftimes_noinf[i-Nexp];
		allprocesses[i] = 1;
	} else if (i >= (Nexp+Ninf)){
		alltimes[i] = remtimes_noinf[i-(Nexp+Ninf)];
		allprocesses[i] = 2;
	}
  }
  
  // Order the alltimes and allprocesses vectors by the order in which events occur
  allprocesses = Rcpp_sort(allprocesses, alltimes);
  std::sort(alltimes.begin(), alltimes.end());
  
  // Set initial time of infectious pressure and initial numbers of susceptiples, exposed, infectious and removed
  double totaltime = 0;
  double CurrSus = M;
  double CurrExp = 0;
  double CurrInf = 0;
  double CurrTot = M;

  
  // Calculate the contribution to the total time of infectious pressure from those that acquire infection during the outbreak
  for (int j = 1; j < alltimes.size(); ++j){
	if (allprocesses[j-1] == 0){
		CurrSus += -1;
		CurrExp += 1;
	} else if (allprocesses[j-1] == 1){
		CurrExp += -1;
		CurrInf += 1;
	} else {
		CurrInf += -1;
		CurrTot += -1;
	}
	totaltime += (CurrSus*CurrInf/CurrTot)*(alltimes[j]-alltimes[j-1]);
  }
  // This last section deals with any gap between last time and Tfinal
  if (allprocesses[no_events-1] == 0){
		CurrSus += -1;
		CurrExp += 1;
	} else if (allprocesses[no_events-1] == 1){
		CurrExp += -1;
		CurrInf += 1;
	} else {
		CurrInf += -1;
		CurrTot += -1;
	}
	totaltime += (CurrSus*CurrInf/CurrTot)*(Tfinal - alltimes[no_events-1]);
  
  // The output
  return totaltime;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
double loglikelihood(Rcpp::DataFrame &eventtimes, double beta_ext, double beta_int, double exp_mean, double exp_shape, double inf_mean, double inf_shape, int M, double Tfinal){

// Calculates the log likelihood
// Inputs: eventtimes - a data frame with all proposed exposure, infection and removal times
// beta_ext, beta_int - transmission rates
// exp_mean, exp_shape - latent mean and shape
// inf_mean, inf_shape - infectious mean and shape
// M - herd size 
// Tfinal - time of cull 
// Outputs: the calculated log likelihood
  
  // Extract columns of event times data frame as vectors
  NumericVector exptimes = eventtimes["Exposure"];
  NumericVector inftimes = eventtimes["Infection"];
  NumericVector remtimes = eventtimes["Removal"];
  
  // Calculate the number of exposures, infections, removals and individuals which remain susceptible
  int Nexp = exptimes.size() - sum(exptimes == inf);
  int Ninf = inftimes.size() - sum(inftimes == inf);
  int Nrem = remtimes.size() - sum(remtimes == inf);
  int Nsusc = M - Nexp;
  
  // Reserve space: the numbers of infected individuals at each of the exposures
  NumericVector no_infected;
  
  // Store the time of the earliest exposure
  double earliest_exposure = min(exptimes);
  
  // Determine the number of infected individuals at the time of each exposure (except the first exposure which must come from an external source)
  for (int i = 0; i < exptimes.size(); ++i){
    
	if(exptimes[i] > earliest_exposure){
		double ref = exptimes[i];
		double ninf_int = std::count_if(inftimes.begin(), inftimes.end(), [ref](double j) { return j < ref; }); 
		double nrem_int = std::count_if(remtimes.begin(), remtimes.end(), [ref](double j) { return j < ref; }); 
		
		double pop_size = M-nrem_int;
		no_infected.push_back((ninf_int-nrem_int)/pop_size);
    }
  }

  
  //#########################################################################
  //# The log-likelihood contributions
  //#########################################################################

  // The log-likelihood contribution from the exposure
  double c1 = Ninf*exp_shape*log(exp_shape/exp_mean) - Ninf*log(tgamma(exp_shape)) + (exp_shape-1)*totaltime_in_stage_log(exptimes, inftimes);

  // The log-likelihood contribution from the removals
  double c2 = Nrem*inf_shape*log(inf_shape/inf_mean) - Nrem*log(tgamma(inf_shape)) + (inf_shape-1)*totaltime_in_stage_log(inftimes, remtimes);
  
  // The log-likelihood contribution from the infections
  double c3 = sum(log(beta_int*no_infected + beta_ext));

  // Calculate the contribution of the integral term 
  double c4a = - (exp_shape/exp_mean)*totaltime_in_stage(exptimes,inftimes);
  double c4b = - (inf_shape/inf_mean)*totaltime_in_stage(inftimes,remtimes);
  double c4c_a = - beta_int*totaltime_infpressure(eventtimes,M,Tfinal);
  double c4c_b =  - beta_ext*(Nsusc*(Tfinal-min(exptimes)) + sum(exptimes-min(exptimes)));
  
  double c5 = 0;
  // Calculate cryptic infection contributions
  for (int i = 0; i < exptimes.size(); ++i){
	  if((inftimes[i] != inf) && (remtimes[i] == inf)){
		  NumericVector crypt_inf_time = NumericVector::create(Tfinal-inftimes[i]);
		  c5 += log(1-pgamma(crypt_inf_time, inf_shape, inf_mean/inf_shape)[0]);
	  }
  }
  
  double c6 = 0;
  // Calculate cryptic exposure contributions
  for (int i = 0; i < exptimes.size(); ++i){
	  if(inftimes[i] == inf){
		  NumericVector crypt_exp_time = NumericVector::create(Tfinal-exptimes[i]);
		  c6 += log(1-pgamma(crypt_exp_time, exp_shape, exp_mean/exp_shape)[0]);
	  }
  }
  
  
  // Output: the log-likelihood
  return c1+c2+c3+c4a+c4b+c4c_a+c4c_b+c5+c6;
  
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
double loglikelihood_beta(Rcpp::DataFrame &eventtimes, double beta_ext, double beta_int, int M, double Tfinal){

// This function calculates parts of the loglikelihood relevant to proposed updates of transmission rates
// Some parts are unnecessary and can be dropped to decrease computation time
    
  // Extract columns of event times data frame as vectors
  NumericVector exptimes = eventtimes["Exposure"];
  NumericVector inftimes = eventtimes["Infection"];
  NumericVector remtimes = eventtimes["Removal"];
  int Nexp = exptimes.size() - sum(exptimes == inf);
  int Nsusc = M - Nexp;

  // Reserve space: the numbers of infected at each of the exposures
  NumericVector no_infected;
  
  double earliest_exposure = min(exptimes);
  
  for (int i = 0; i < exptimes.size(); ++i){
    
	if(exptimes[i] > earliest_exposure){
		double ref = exptimes[i];
		double ninf_int = std::count_if(inftimes.begin(), inftimes.end(), [ref](double j) { return j < ref; }); 
		double nrem_int = std::count_if(remtimes.begin(), remtimes.end(), [ref](double j) { return j < ref; }); 
		double pop_size = M-nrem_int;

		no_infected.push_back((ninf_int-nrem_int)/pop_size);
    }
  }

  
  //#########################################################################
  //# The log-likelihood contributions
  //#########################################################################

    // The log-likelihood contribution from the infections
  double c3 = sum(log(beta_int*no_infected + beta_ext));

  // Calculate the contribution of the integral term 
  double c4c = - beta_int*totaltime_infpressure(eventtimes,M,Tfinal) - beta_ext*(Nsusc*(Tfinal-min(exptimes)) + sum(exptimes-min(exptimes)));
  
  // Output: the log-likelihood
  return c3+c4c;
  
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
double loglikelihood_exp_mean(Rcpp::DataFrame &eventtimes, double exp_mean, double exp_shape, double Tfinal){

// This function calculates parts of the loglikelihood relevant to proposed updates of the latent mean
// Some parts are unnecessary and can be dropped to decrease computation time
  
  int Nexp = eventtimes.nrow();
  NumericVector inftimes = eventtimes["Infection"];
  NumericVector exptimes = eventtimes["Exposure"];
  NumericVector remtimes = eventtimes["Removal"];
  int Ninf = inftimes.size() - sum(inftimes == inf);

    // The log-likelihood contribution from the exposure
  double c1 = Ninf*exp_shape*log(exp_shape/exp_mean);

  // Calculate the contribution of the integral term 
  double c4 = -(exp_shape/exp_mean)*totaltime_in_stage(exptimes,inftimes);
  
  double c6 = 0;
  // Calculate cryptic exposure contributions
  for (int i = 0; i < Nexp; ++i){
	if(inftimes[i] == inf){
		NumericVector crypt_exp_time = NumericVector::create(Tfinal-exptimes[i]);
		c6 += log(1-pgamma(crypt_exp_time, exp_shape, exp_mean/exp_shape)[0]);
	}
  }
  
  // Output: the log-likelihood
  return c1+c4+c6;
  
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
double loglikelihood_exp_shape(Rcpp::DataFrame &eventtimes, double exp_mean, double exp_shape, double Tfinal){

// This function calculates parts of the loglikelihood relevant to proposed updates of latent shape
// Some parts are unnecessary and can be dropped to decrease computation time  

  int Nexp = eventtimes.nrow();
  NumericVector inftimes = eventtimes["Infection"];
  NumericVector exptimes = eventtimes["Exposure"];
  NumericVector remtimes = eventtimes["Removal"];
  int Ninf = inftimes.size() - sum(inftimes == inf);

  //#########################################################################
  //# The log-likelihood contributions
  //#########################################################################

  // The log-likelihood contribution from the exposure
  double c1 = Ninf*exp_shape*log(exp_shape/exp_mean) - Ninf*log(tgamma(exp_shape)) + (exp_shape-1)*totaltime_in_stage_log(exptimes,inftimes);
  
  double c4 = -(exp_shape/exp_mean)*totaltime_in_stage(exptimes,inftimes);

  double c6 = 0;
  // Calculate cryptic exposure contributions
  for (int i = 0; i < Nexp; ++i){
	  if(inftimes[i] == inf){
		  NumericVector crypt_exp_time = NumericVector::create(Tfinal-exptimes[i]);
		  c6 += log(1-pgamma(crypt_exp_time, exp_shape, exp_mean/exp_shape)[0]);
	  }
  }
  
  // Output: the log-likelihood
  return c1+c4+c6;
  
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
double loglikelihood_inf_mean(Rcpp::DataFrame &eventtimes, double inf_mean, double inf_shape, double Tfinal){

// This function calculates parts of the loglikelihood relevant to proposed updates of the infectious mean
// Some parts are unnecessary and can be dropped to decrease computation time  

  int Nexp = eventtimes.nrow();
  NumericVector inftimes = eventtimes["Infection"];
  NumericVector exptimes = eventtimes["Exposure"];
  NumericVector remtimes = eventtimes["Removal"];
  int Nrem = remtimes.size() - sum(remtimes == inf);

    // The log-likelihood contribution from the exposure
  double c2 = Nrem*inf_shape*log(inf_shape/inf_mean);

  // Calculate the contribution of the integral term 
  double c4 = -(inf_shape/inf_mean)*totaltime_in_stage(inftimes,remtimes);
  
  double c5 = 0;
  // Calculate cryptic exposure contributions
  for (int i = 0; i < Nexp; ++i){
	  if((inftimes[i] != inf) && (remtimes[i] == inf)){
		  NumericVector crypt_inf_time = NumericVector::create(Tfinal-inftimes[i]);
		  c5 += log(1-pgamma(crypt_inf_time, inf_shape, inf_mean/inf_shape)[0]);
	  }
  }
  
  // Output: the log-likelihood
  return c2+c4+c5;
  
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
double loglikelihood_inf_shape(Rcpp::DataFrame &eventtimes, double inf_mean, double inf_shape, double Tfinal){

// This function calculates parts of the loglikelihood relevant to proposed updates of the infectious shape
// Some parts are unnecessary and can be dropped to decrease computation time  

  int Nexp = eventtimes.nrow();
  NumericVector inftimes = eventtimes["Infection"];
  NumericVector exptimes = eventtimes["Exposure"];
  NumericVector remtimes = eventtimes["Removal"];
  int Nrem = remtimes.size() - sum(remtimes == inf);

  //#########################################################################
  //# The log-likelihood contributions
  //#########################################################################

  // The log-likelihood contribution from the exposure
  double c2 = Nrem*inf_shape*log(inf_shape/inf_mean) - Nrem*log(tgamma(inf_shape)) + (inf_shape-1)*totaltime_in_stage_log(inftimes,remtimes);
  
  double c4 = -(inf_shape/inf_mean)*totaltime_in_stage(inftimes,remtimes);

  double c5 = 0;
  // Calculate cryptic infection contributions
  for (int i = 0; i < Nexp; ++i){
	  if((inftimes[i] != inf) && (remtimes[i] == inf)){
		  NumericVector crypt_inf_time = NumericVector::create(Tfinal-inftimes[i]);
		  c5 += log(1-pgamma(crypt_inf_time, inf_shape, inf_mean/inf_shape)[0]);
	  }
  }
  
  // Output: the log-likelihood
  return c2+c4+c5;
  
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
List update_beta_ext(Rcpp::DataFrame &completedata, double beta_ext, double beta_ext_sigma, double beta_int, double beta_ext_rate, double beta_ext_shape, int M, double Tfinal, double accept_beta_ext){
  
// This function proposes and accepts or rejects a new value for the external transmission rate
// Inputs: completedata - data frame with all proposed exposure, infection and removal times
// beta_ext, beta_int - external and within-herd transmission rates
// beta_ext_sigma - variance of the proposal distribution for beta_ext
// beta_ext_rate, beta_ext_shape - parameters of the prior distribution for beta_ext
// M - herd size
// Tfinal - time of cull
// accept_beta_ext - current number of acceptances for updates to beta_ext
// Outputs: beta_ext - new value of beta_ext
// accept_beta_ext - new number of acceptances  

	// Propose a new value for beta_ext and create variables for new and old log likelihood
	double prop_beta_ext = rnorm(1,beta_ext,beta_ext_sigma)[0];
	double ll_new;
	double ll_old;
  
    // if proposed transmission rate is negative, make it positive
	if(prop_beta_ext < 0){
		prop_beta_ext = -prop_beta_ext;
	}
	// calculate new and old log likelihood
	ll_old = loglikelihood_beta(completedata,beta_ext,beta_int, M, Tfinal)+(beta_ext_shape-1)*log(beta_ext)-beta_ext_rate*beta_ext;
	ll_new = loglikelihood_beta(completedata,prop_beta_ext,beta_int, M, Tfinal)+(beta_ext_shape-1)*log(prop_beta_ext)-beta_ext_rate*prop_beta_ext;
	
	// Choose to accept or reject proposal and update acceptances as appropriate
	if(ll_new != minusinf){
		if(log(runif(1,0,1)[0]) < (ll_new-ll_old)){
			beta_ext = prop_beta_ext;
			accept_beta_ext += 1;
		}
	}
	
  // return new parameter value and number of acceptances	
  return List::create(beta_ext, accept_beta_ext);  
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
List update_beta_int(Rcpp::DataFrame &completedata, double beta_int, double beta_int_sigma, double beta_ext, double beta_int_rate, double beta_int_shape, int M, double Tfinal, double accept_beta_int){
  
// This function proposes and accepts or rejects a new value for beta_int
// It works the same way as the function for beta_ext above
  
    double prop_beta_int = rnorm(1,beta_int,beta_int_sigma)[0];
	double ll_new;
	double ll_old;
  
	if(prop_beta_int < 0){
		prop_beta_int = -prop_beta_int;
	}
	ll_old = loglikelihood_beta(completedata,beta_ext,beta_int, M, Tfinal)+(beta_int_shape-1)*log(beta_int)-beta_int_rate*beta_int;
	ll_new = loglikelihood_beta(completedata,beta_ext,prop_beta_int, M, Tfinal)+(beta_int_shape-1)*log(prop_beta_int)-beta_int_rate*prop_beta_int;	
	
	if(ll_new != minusinf){
		if(log(runif(1,0,1)[0]) < (ll_new-ll_old)){
			beta_int = prop_beta_int;
			accept_beta_int += 1;
		}
	}
	
  return List::create(beta_int, accept_beta_int);  
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
List update_exp_mean(Rcpp::List &completedata, double exp_mean, double exp_mean_sigma, double exp_shape, double exp_mean_rate, double exp_mean_shape, Rcpp::NumericVector &Tfinal, double accept_exp_mean){
  
// This function proposes and accepts or rejects a new value for exp_mean
// It works the same way as the function for beta_ext above except the log likelihood is calculated across all the herds
  
  double prop_exp_mean = rnorm(1,exp_mean,exp_mean_sigma)[0];
  NumericVector ll_old (9);
  NumericVector ll_new (9);
  
	if(prop_exp_mean < 0){
		prop_exp_mean = -prop_exp_mean;
	}
	for(int i = 0; i < 9; ++i){
		DataFrame completedata_temp = completedata[i];
		double Tfinal_temp = Tfinal[i];
		ll_old[i] = loglikelihood_exp_mean(completedata_temp,exp_mean,exp_shape, Tfinal_temp)+(exp_mean_shape-1)*log(exp_mean)-exp_mean_rate*exp_mean;
		ll_new[i] = loglikelihood_exp_mean(completedata_temp,prop_exp_mean,exp_shape, Tfinal_temp)+(exp_mean_shape-1)*log(prop_exp_mean)-exp_mean_rate*prop_exp_mean;
	}	
	
	double tot_ll_old = sum(ll_old);
	double tot_ll_new = sum(ll_new);
	if(tot_ll_new != minusinf){
		if(log(runif(1,0,1)[0]) < (tot_ll_new-tot_ll_old)){
			exp_mean = prop_exp_mean;
			accept_exp_mean += 1;
		}
	}
	
  return List::create(exp_mean, accept_exp_mean);  
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
List update_exp_shape(Rcpp::List &completedata, double exp_mean, double exp_shape, double exp_shape_sigma, double exp_shape_rate, double exp_shape_shape, Rcpp::NumericVector &Tfinal, double accept_exp_shape){

// This function proposes and accepts or rejects a new value for exp_shape
// It works the same way as the function for beta_ext above except the log likelihood is calculated across all the herds
  
  double prop_exp_shape = rnorm(1,exp_shape,exp_shape_sigma)[0];
  NumericVector ll_old (9);
  NumericVector ll_new (9);
  
  if(prop_exp_shape < 0){
	prop_exp_shape = -prop_exp_shape;
  }
  
  for(int i = 0; i < 9; ++i){
	DataFrame completedata_temp = completedata[i];
	double Tfinal_temp = Tfinal[i];
	ll_old[i] = loglikelihood_exp_shape(completedata_temp, exp_mean, exp_shape, Tfinal_temp)+(exp_shape_shape-1)*log(exp_shape)-exp_shape_rate*exp_shape;
	ll_new[i] = loglikelihood_exp_shape(completedata_temp, exp_mean, prop_exp_shape, Tfinal_temp)+(exp_shape_shape-1)*log(prop_exp_shape)-exp_shape_rate*prop_exp_shape;
  }
  
  double tot_ll_old = sum(ll_old);
  double tot_ll_new = sum(ll_new);
  if(tot_ll_new != minusinf){
	if(log(runif(1,0,1)[0]) < (tot_ll_new-tot_ll_old)){
		exp_shape = prop_exp_shape;
		accept_exp_shape += 1;
	}
  }
  return List::create(exp_shape, accept_exp_shape);  

}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
List update_inf_mean(Rcpp::List &completedata, double inf_mean, double inf_mean_sigma, double inf_shape, double inf_mean_rate, double inf_mean_shape, Rcpp::NumericVector &Tfinal, double accept_inf_mean){
  
// This function proposes and accepts or rejects a new value for inf_mean
// It works the same way as the function for beta_ext above except the log likelihood is calculated across all the herds 
  
	double prop_inf_mean = rnorm(1,inf_mean,inf_mean_sigma)[0];
	NumericVector ll_old (9);
	NumericVector ll_new (9);
  
	if(prop_inf_mean < 0){
		prop_inf_mean = -prop_inf_mean;
	}
	// Inverse gamma prior
	for(int i = 0; i < 9; ++i){
		DataFrame completedata_temp = completedata[i];
		double Tfinal_temp = Tfinal[i];
		ll_old[i] = loglikelihood_inf_mean(completedata_temp,inf_mean,inf_shape,Tfinal_temp)+(inf_mean_shape-1)*log(inf_mean)-inf_mean_rate*inf_mean;
		ll_new[i] = loglikelihood_inf_mean(completedata_temp,prop_inf_mean,inf_shape,Tfinal_temp)+(inf_mean_shape-1)*log(prop_inf_mean)-inf_mean_rate*prop_inf_mean;
	}
  
	double tot_ll_old = sum(ll_old);
	double tot_ll_new = sum(ll_new);	
	if(tot_ll_new != minusinf){
		if(log(runif(1,0,1)[0]) < (tot_ll_new-tot_ll_old)){
			inf_mean = prop_inf_mean;
			accept_inf_mean += 1;
		}
	}  
  return List::create(inf_mean, accept_inf_mean);  
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
List update_inf_shape(Rcpp::List &completedata, double inf_mean, double inf_shape, double inf_shape_sigma, double inf_shape_rate, double inf_shape_shape, Rcpp::NumericVector &Tfinal, double accept_inf_shape){
  
// This function proposes and accepts or rejects a new value for inf_shape
// It works the same way as the function for beta_ext above except the log likelihood is calculated across all the herds
  
	double prop_inf_shape = rnorm(1,inf_shape,inf_shape_sigma)[0];
	NumericVector ll_old (9);
	NumericVector ll_new (9);
  
	if(prop_inf_shape < 0){
		prop_inf_shape = -prop_inf_shape;
	}
	// Inverse gamma prior
	for(int i = 0; i < 9; ++i){
		DataFrame completedata_temp = completedata[i];
		double Tfinal_temp = Tfinal[i];	
		ll_old[i] = loglikelihood_inf_shape(completedata_temp,inf_mean,inf_shape,Tfinal_temp)+(inf_shape_shape-1)*log(inf_shape)-inf_shape_rate*inf_shape;
		ll_new[i] = loglikelihood_inf_shape(completedata_temp,inf_mean,prop_inf_shape,Tfinal_temp)+(inf_shape_shape-1)*log(prop_inf_shape)-inf_shape_rate*prop_inf_shape;
	}
  
	double tot_ll_old = sum(ll_old);
	double tot_ll_new = sum(ll_new);	
	if(tot_ll_new != minusinf){
		if(log(runif(1,0,1)[0]) < (tot_ll_new-tot_ll_old)){
			inf_shape = prop_inf_shape;
			accept_inf_shape += 1;
		}
	}  
  return List::create(inf_shape, accept_inf_shape);  
}

// [[Rcpp::export]]
DataFrame update_times(Rcpp::DataFrame& completedata, double beta_ext, double beta_int, double exp_mean, double exp_shape, double inf_mean, double inf_shape, int M, double Tfinal){

// This function selects an individual in the population and proposes an update to its exposure/infection time(s) as appropriate
// Inputs: completedata - a data frame with all proposed exposure, infection and removal times
// beta_ext, beta_int - transmission rates
// exp_mean, exp_shape - latent mean and shape
// inf_mean, inf_shape - infectious mean and shape
// M - herd size 
// Tfinal - time of cull 
// Outputs: an updated set of exposure, infection and removal times

  // Store exposure, infection and removal times as vectors and create clones of these for proposing new times
  NumericVector remtimes = completedata["Removal"];
  NumericVector exptimes = completedata["Exposure"];
  NumericVector inftimes = completedata["Infection"];
  NumericVector new_remtimes = clone(remtimes);
  NumericVector new_exptimes = clone(exptimes);
  NumericVector new_inftimes = clone(inftimes);
  
  // Store the number of exposures
  int Nexp = exptimes.size();
  
  // Calculate the log likelihood before updating any times
  double loglikelihood_old = loglikelihood(completedata,beta_ext,beta_int,exp_mean,exp_shape,inf_mean,inf_shape,M,Tfinal);

  // Choose an individual at random to update
  IntegerVector prop_update = sample(M,1)-1;
  int j = prop_update[0];
  
  // Create these variables to calculate acceptance probabilities
  double scale1 = 1;
  double scale2 = 1;
  
  if (j > (Nexp-1)){ 
      // What to do if we pick susceptible individual
	  int prob_exp = rbinom(1, 1, 0.5)[0];
	  if (prob_exp == 1){ 
	      // Make it a cryptic exposure
		  new_exptimes.push_back(runif(1,min(inftimes),Tfinal)[0]);
		  new_inftimes.push_back(inf);
		  new_remtimes.push_back(inf);
		  scale1 = 2*(Tfinal-min(inftimes));
		  scale2 = 3;
	  } else { 
	      // Make it a cryptic infection
		  new_remtimes.push_back(inf);
		  new_inftimes.push_back(runif(1,min(inftimes),Tfinal)[0]);
		  new_exptimes.push_back(new_inftimes[Nexp]-rgamma(1,exp_shape,exp_mean/exp_shape)[0]);
		  NumericVector exp_dur_new = NumericVector::create(new_inftimes[Nexp]-new_exptimes[Nexp]);
		  scale1 = 2*(Tfinal-min(inftimes));
		  scale2 = 3*dgamma(exp_dur_new,exp_shape,exp_mean/exp_shape)[0];
	  }
  } else if (remtimes[j] != inf) { 
      // What to do if we pick removed individual
	  int prob_exp = rbinom(1, 1, 0.5)[0];
	  if (prob_exp == 1){ 
	      // Update the exposure time
		  new_exptimes[j] = inftimes[j]-rgamma(1,exp_shape,exp_mean/exp_shape)[0];
		  NumericVector exp_dur_old = NumericVector::create(inftimes[j]-exptimes[j]);
		  NumericVector exp_dur_new = NumericVector::create(inftimes[j]-new_exptimes[j]);
		  scale1 = dgamma(exp_dur_old,exp_shape,exp_mean/exp_shape)[0];
		  scale2 = dgamma(exp_dur_new,exp_shape,exp_mean/exp_shape)[0];
	  } else { 
	      // Update the infection time
		  new_inftimes[j] = remtimes[j] - rgamma(1,inf_shape,inf_mean/inf_shape)[0];
		  NumericVector inf_dur_old = NumericVector::create(remtimes[j]-inftimes[j]);
		  NumericVector inf_dur_new = NumericVector::create(remtimes[j]-new_inftimes[j]);
		  scale1 = dgamma(inf_dur_old,inf_shape,inf_mean/inf_shape)[0];
		  scale2 = dgamma(inf_dur_new,inf_shape,inf_mean/inf_shape)[0];
	  }
  } else if (inftimes[j] != inf) { 
      // What to do if we pick a cryptic infection
	  IntegerVector probs = seq_len(3);
	  int index = sample(probs, 1)[0];
	  if (index == 1){ 
	      // delete infection and exposure
		  NumericVector exp_dur_old = NumericVector::create(inftimes[j]-exptimes[j]);
		  scale1 = 3*dgamma(exp_dur_old,exp_shape,exp_mean/exp_shape)[0];
		  scale2 = 2*(Tfinal-min(inftimes));
		  new_exptimes.erase(j);
		  new_inftimes.erase(j);
		  new_remtimes.erase(j);
	  } else if (index == 2){ 
	      // change to cryptic exposure
	      new_inftimes[j] = inf;
		  scale1 = 1;
		  scale2 = (Tfinal - exptimes[j]);
	  } else {
		  // update exposure and infection time
		  new_inftimes[j] = runif(1, min(inftimes), Tfinal)[0];
		  new_exptimes[j] = new_inftimes[j]-rgamma(1,exp_shape,exp_mean/exp_shape)[0];
		  NumericVector exp_dur_old = NumericVector::create(inftimes[j]-exptimes[j]);
		  NumericVector exp_dur_new = NumericVector::create(new_inftimes[j]-new_exptimes[j]);
		  scale1 = dgamma(exp_dur_old,exp_shape,exp_mean/exp_shape)[0];
		  scale2 = dgamma(exp_dur_new,exp_shape,exp_mean/exp_shape)[0];
	  }
  } else {
	  // What if we pick a cryptic exposure
	  IntegerVector probs = seq_len(3);
	  int index = sample(probs, 1)[0];
	  if (index == 1){
		  // delete exposure
		  new_exptimes.erase(j);
		  new_inftimes.erase(j);
		  new_remtimes.erase(j);
		  scale1 = 3;
		  scale2 = 2*(Tfinal-min(inftimes));
	  } else if (index == 2){
		  // change to cryptic infection
		  new_inftimes[j] = runif(1, exptimes[j], Tfinal)[0];
		  scale1 = (Tfinal - exptimes[j]);
		  scale2 = 1;
	  } else {
		  // update exposure time
		  new_exptimes[j] = runif(1, min(exptimes), Tfinal)[0];
	  }
  }
	  
  // Create a data frame with the proposed set of times
  DataFrame completedata_new = DataFrame::create(Named("Exposure")= new_exptimes, Named("Infection") = new_inftimes, Named("Removal") = new_remtimes);
  
  // Calculate the log likelihood of the proposed set of times
  double loglikelihood_new = loglikelihood(completedata_new,beta_ext,beta_int,exp_mean,exp_shape,inf_mean,inf_shape,M,Tfinal);
  
  // Determine whether to accept or reject the new set of times
  if (loglikelihood_new != minusinf){
	if (log(runif(1,0,1)[0]) < (loglikelihood_new+log(scale1)-loglikelihood_old-log(scale2))){
		completedata = completedata_new;
	}
  }
  
  // Output the updated set of times
  return completedata;
}