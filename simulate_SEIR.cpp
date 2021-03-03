#include <Rcpp.h>
#include <cstdlib>
#include <cmath> 
using namespace Rcpp;

// [[Rcpp::export]]
List SEIR_IBM(int N, int I0, double exp_shape, double exp_rate, double inf_shape, double inf_rate, double beta_int, double beta_ext, double t0, double tfinal){
// Function to simulate from an SEIR individual based model
// Inputs: N - herd size
// I0 - number of initial infections from outside the herd
// exp_shape, exp_rate - parameters of the latent distribution
// inf_shape, inf_rate - parameters of the infectious distribution
// beta_int, beta_ext - transmission rates
// t0 - start time for the outbreak
// tfinal - time of cull
// Outputs: vectors with exposure, infection and removal times
	
	// Set up variables for current time and current numbers of susc, exp, inf and rem
	double current_time = t0;
	int current_S = N-I0;
	int current_E = I0;
	int current_I = 0;
	int current_R = 0;
	
	// Create vectors showing the composition of the population each time a transition occurs and initialise with starting state
	std::vector<double> vec_time;
	std::vector<int> vec_S, vec_E, vec_I, vec_R, vec_Ind;
	vec_time.push_back(t0);
	vec_S.push_back(N-I0);
	vec_E.push_back(I0);
	vec_I.push_back(0);
	vec_R.push_back(0);
	
	// Create vector giving state of each individual in the population
	std::vector<int> state_vec(N, 0);
	int Inf_index = rand() % N;
	state_vec[Inf_index] = 1;
	vec_Ind.push_back(Inf_index);
	
	// Create vectors of exp, inf and rem times and set set all values to after the end of the epidemic
	std::vector<double> exposure_times(N, tfinal+1);
	std::vector<double> infection_times(N, tfinal+1);
	std::vector<double> removal_times(N, tfinal+1);
	
	// Initialise the exposure, infection and removal/removal time vectors 
	for(int i = 0; i < N; ++i){
		if(state_vec[i] == 0){
			if(current_I == 0){
				exposure_times[i] = tfinal+1;
			} else {
				double dt = rexp(1, (beta_int*current_I)/(N-current_R) + beta_ext)[0];
				exposure_times[i] = current_time + dt;
			}
		} else if(state_vec[i] == 1){
			double dt = rgamma(1, exp_shape, 1/exp_rate)[0];
			infection_times[i] = current_time + dt;
		} else if(state_vec[i] == 2){
			double dt = rgamma(1, inf_shape, 1/inf_rate)[0];
			removal_times[i] = current_time + dt;
		}
	}
	
	while((current_time < tfinal) && ((current_I+current_E) > 0)){
		// Determine what the next event will be.
		int minPos_exp = 0;
		for (unsigned i = 0; i < exposure_times.size(); ++i){
			if (exposure_times[i] < exposure_times[minPos_exp]){
				minPos_exp = i;
			}
		}
		int minPos_inf = 0;
		for (unsigned i = 0; i < infection_times.size(); ++i){
			if (infection_times[i] < infection_times[minPos_inf]){
				minPos_inf = i;
			}
		}
		int minPos_rec = 0;
		for (unsigned i = 0; i < removal_times.size(); ++i){
			if (removal_times[i] < removal_times[minPos_rec]){
				minPos_rec = i;
			}
		}
		std::vector<double> min_positions;
		min_positions.push_back(exposure_times[minPos_exp]);
		min_positions.push_back(infection_times[minPos_inf]);
		min_positions.push_back(removal_times[minPos_rec]);
		int minPos_combined = 0;
		for (unsigned i = 0; i < min_positions.size(); ++i){
			if (min_positions[i] < min_positions[minPos_combined]){
				minPos_combined = i;
			}
		}

		// Update the current time
		current_time = min_positions[minPos_combined];

		// Depending on next event, update state vectors and event times
		if(minPos_combined == 0){
			state_vec[minPos_exp] = 1;
			exposure_times[minPos_exp] = tfinal+1;
			infection_times[minPos_exp] = current_time + rgamma(1, exp_shape, 1/exp_rate)[0];
		} else if(minPos_combined == 1){
			state_vec[minPos_inf] = 2;
			infection_times[minPos_inf] = tfinal+1;
			removal_times[minPos_inf] = current_time + rgamma(1, inf_shape, 1/inf_rate)[0];
		} else if(minPos_combined == 2){
			state_vec[minPos_rec] = 3;
			removal_times[minPos_rec] = tfinal+1;
		}
		
		// Update counts of numbers in each state for each population
		current_S = std::count (state_vec.begin(), state_vec.end(), 0);
		current_E = std::count (state_vec.begin(), state_vec.end(), 1);
		current_I = std::count (state_vec.begin(), state_vec.end(), 2);
		current_R = std::count (state_vec.begin(), state_vec.end(), 3);
		
		vec_S.push_back(current_S);
		vec_E.push_back(current_E);
		vec_I.push_back(current_I);
		vec_R.push_back(current_R);
		if(minPos_combined == 0){
			vec_Ind.push_back(minPos_exp);
		} else if(minPos_combined == 1){
			vec_Ind.push_back(minPos_inf);
		} else if(minPos_combined == 2){
			vec_Ind.push_back(minPos_rec);
		}
		vec_time.push_back(current_time);
		
		if((current_I+current_E) == 0){
			vec_S.push_back(current_S);
			vec_E.push_back(current_E);
			vec_I.push_back(current_I);
			vec_R.push_back(current_R);
			vec_time.push_back(tfinal);
		}
		
		// Update the exposure times if the number of infectious individuals has changed
		if((minPos_combined == 1) || (minPos_combined == 2)){
			for(int i = 0; i < N; ++i){
				if(state_vec[i] == 0){
					if(current_I == 0){
						exposure_times[i] = tfinal+1;
					} else {
						double dt = rexp(1, beta_int*current_I/(N-current_R) + beta_ext)[0];
						exposure_times[i] = current_time + dt;
					}
				}
			}
		}
	}
	return List::create(vec_time, vec_S, vec_E, vec_I, vec_R, vec_Ind);
}