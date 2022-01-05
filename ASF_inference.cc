#include <iostream>
#include <sstream>
#include <math.h>
#include <map>
#include <algorithm>
#include <vector>
#include <iterator>
#include <fstream>
#include <iomanip>
#include <string>
#include <random>

using namespace std;

#include "struct.hh"

enum Model{ MULTI_INIT, EXT_INF};

Model model = EXT_INF;						 // Chosen EXT_INF here which switches beta_1 on.  If MULTI_INIT is chosen with index_max=1 then this would give a single index case and no external transmission.
	
const auto phi_L = 1.0;                      // Inverse temperature multipling the infection process
const auto phi_G = phi_L;                    // Inverse temperature multipling the transition process

const auto nherd = 1u;                       // The number of herds

const auto alpha = 0.5;                      // Determines an exponentially distributed prior on index cases - irrelavant for ASF data as we will fix to have only one index case
const auto index_max = 1;                    // The maximum limit on the number of index cases - the above line only matters if this is > 1

const auto nupd = 40;			             // The multi_proposals on each iteration

const auto nsamp = 2000000u;                  // The number of MCMC samples
const auto nburnin = nsamp/10;                // The number step for burnin 
unsigned int s;                              // The sample number

const int NO_INF = -10000;                   // The log likelihood of getting exposed when no infected

bool simulated = false;				// Use simulated or observed removal times
//const vector <int> herd_sizes = { 1614,1949,1753,1833,1320,600,600,600,2145};	// Vector option to allow multiple herds to be entered at once if parameters are pooled across herds
const vector <int> herd_sizes = { 1614};										// Herd size
const vector <double> pct_seen = { 0.15};										// Proportion of the herd which dies before cull (only relevant if data is simulated)

// Herd times can be entered as a 2 dimensional array in case running multiple herds at once.  The code as shown below would run for herd 1 only.
std::vector<std::vector<double>> remtimes = {{ 1.0,4.0,7.0,8.0,12.0,13.0,15.0,15.0,15.0,15.0,16.0,16.0,16.0,16.0,16.0,16.0,17.0,17.0,17.0,17.0,17.0,18.0,18.0,18.0,18.0,18.0,18.0,19.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,22.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,24.0,24.0,24.0,24.0,24.0,24.0,24.0,24.0,24.0,24.0,24.0,24.0,24.0}};//,
										  //{ 2.0,3.0,4.0,5.0,5.0,6.0,6.0,6.0,7.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,9.0,9.0,9.0,9.0,9.0,9.0,10.0,10.0,10.0,10.0,11.0,11.0,11.0,11.0,12.0,12.0,12.0,12.0,12.0,12.0,12.0,12.0,12.0,12.0,12.0,13.0,13.0,13.0,14.0,14.0,14.0,14.0,14.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,16.0,16.0,17.0,17.0,17.0},
										  //{ 5.0,6.0,7.0,8.0,8.0,9.0,9.0,10.0,10.0,10.0,10.0,10.0,11.0,11.0,11.0,11.0,12.0,12.0,12.0,12.0,12.0,12.0,12.0,12.0,12.0,13.0,14.0,14.0,14.0,15.0,15.0,15.0,16.0,16.0,16.0,16.0,16.0,17.0,17.0,17.0,17.0,17.0,17.0},
										  //{ 1.0,2.0,2.0,3.0,3.0,3.0,3.0,3.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,5.0,5.0,5.0,5.0,5.0,6.0,6.0,6.0,6.0,6.0,6.0,6.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,9.0,9.0,9.0,9.0,10.0,10.0,10.0,10.0,10.0,11.0,11.0,11.0,11.0,12.0,12.0,12.0,12.0,12.0,12.0,12.0,12.0,12.0,12.0,13.0,13.0,14.0,14.0,14.0,14.0,14.0,14.0,14.0,14.0,14.0,14.0,14.0,14.0,14.0,14.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,16.0,16.0,16.0},
										  //{ 3.0,3.0,3.0,3.0,4.0,4.0,4.0,4.0,4.0,4.0,5.0,5.0,5.0,5.0,6.0,6.0,6.0,6.0,7.0,7.0,7.0,7.0,7.0,7.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,9.0,9.0,9.0,9.0,9.0,9.0,10.0,11.0,11.0,11.0,12.0,12.0,12.0,12.0,12.0,13.0,13.0,13.0,13.0,13.0,13.0,13.0,13.0,13.0,14.0,14.0,14.0,14.0,14.0,14.0,14.0,14.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0},
										  //{ 5.0,6.0,7.0,7.0,8.0,8.0,9.0,12.0,13.0,14.0,14.0,14.0,14.0,14.0,14.0,14.0,14.0},
										  //{ 1.0,1.0,3.0,5.0,5.0,6.0,7.0,7.0,7.0,7.0,8.0,8.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,11.0,11.0,11.0,11.0,11.0,11.0,11.0,11.0,12.0,12.0,12.0,12.0,12.0,12.0,13.0,13.0,13.0,13.0,13.0},
										  //{ 6.0,6.0,8.0,8.0,8.0,9.0,9.0,9.0,10.0,10.0,10.0,10.0,11.0,11.0,11.0,13.0},
										  //{ 2.0,2.0,2.0,2.0,2.0,2.0,2.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,4.0,4.0,4.0,5.0,5.0,7.0,7.0,7.0,7.0,7.0,8.0,10.0,10.0,10.0,10.0,10.0,10.0,12.0,13.0,13.0,13.0,13.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,17.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,19.0,19.0,19.0,19.0,19.0,19.0,19.0,19.0,19.0,19.0,19.0,19.0,19.0,19.0,19.0,19.0,19.0,19.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,21.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,22.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,24.0,24.0,24.0,24.0,24.0,24.0,24.0,24.0,24.0,24.0,24.0,24.0}};

vector <Herd> herd(nherd); // Stores information 
vector <Herd> herd_ppc(nherd);	

vector <Param> param;      // A list of all the model parameters

int mu_L;                  // The parameter giving the mean of the latent period
int sh_L;                  // The parameter giving the standard deviation of the latent period
int mu_I;                  // The parameter giving the mean of the infectious period
int sh_I;                  // The parameter giving the standard deviation of the infectious period

double Pr;                 // The log of the prior

double time_total;
double time_likelihood = 0;

bool is_ppc = false;		// Switched between true and false for generating posterior predictive checks rather than proposing new times.

int slurm_id;				// Parameters below are for running jobs in parallel on HPC using slurm

std::string filename;
std::string filename2;
std::string filename3;

ofstream trace;
ofstream ppc;
ofstream herd_data;

default_random_engine generator;

bool EV_ord (Event ev1, Event ev2){ return (ev1.t < ev2.t); };  

int add_param(string name, PriorType prior_type, double prior_val1, double prior_val2, double upper_bound);
void mcmc_initialise();
void param_prior_sample();
void param_prior_init();
double prior();
void randomwalk_proposal();
void MBP_proposal(MBPType type);
double choose(int N, int n);
void mcmc_diagnostics();
void check(int num);
		
void trace_init();
void herd_init();
void herd_write();
void ppc_init();
void trace_plot();
void ppc_plot();

TransParam set_trans_param();
double gamma_sample(const double mu, const double sh);
double trunc_gamma_sample(const double mu, const double sh, const string nm);
double beta_sample(const double a, const double b);
double gamma_probability(const double x, const double mu, const double sh);
double trunc_gamma_probability(const double x, const double mu, const double sh, const string nm);
double beta_probability(const double x, const double a, const double b);
double exp_sample(const double rate);
double exp_probability(const double x, const double mu);
double normal_sample(const double mu, const double sd);
double normal_probability(const double x, const double mean, const double var);
double uniform_sample(const double lwr, const double upr);
void emsg(const string& msg);
double ran();
void choose_init();

int main(int argc, char** argv)
{
	if(argc != 2) emsg("Must provide a seed");
	
	// Initialise all the model parameters
	
	std::stringstream ss;					// Define filename based on job number
	std::stringstream ss2;
	std::stringstream ss3;
	
	ss << "trace_h1_" << argv[1] << ".txt";			// Open output files
	filename = ss.str();
	ss2 << "herd_data_h1_" << argv[1] << ".txt";
	filename2 = ss2.str();
	ss3 << "ppc_h1_" << argv[1] << ".txt";
	filename3 = ss3.str();
	trace.open(filename);
	herd_data.open(filename2);
	ppc.open(filename3);
	
	slurm_id = std::stoi(argv[1]);			// Set the seed so that results can be replicated and each chain is different
	
	mu_L = add_param("mu_L",TRUNC_GAMMA_PRIOR,6.25,10.0,UNSET);	// Set up some parameters
	sh_L = add_param("sh_L",TRUNC_GAMMA_PRIOR,19.39,5.0,UNSET);
	
	mu_I = add_param("mu_I",TRUNC_GAMMA_PRIOR,9.12,10.0,UNSET);
	sh_I = add_param("sh_I",TRUNC_GAMMA_PRIOR,22.2,5.0,UNSET);
	
	generator.seed(5);
	
	for(auto h = 0u; h < nherd; h++){		// Set up each herd
		herd[h].initialise(h);
		is_ppc = true;
		herd_ppc[h].initialise(h);
		is_ppc = false;
	}
	
	herd_init();							// Sets up an output file for data on removal times in each herd			
	
	if(simulated){                                                    // Simulates some data and writes to a file
		param_prior_init();            								  // Sets the parameter values
		for(auto &he : herd) he.simulate(int(pct_seen[he.index]*herd_sizes[he.index]));
		//for(auto &he : herd) he.simulate(int(1*herd_sizes[he.index]));
		herd_write();
	} 
	else {                                                            // Loads some data
		for(auto &he : herd) he.load_data();
	}
	
	generator.seed(slurm_id+10);
	
	param_prior_sample();											  // Samples values for each parameter from the priors  
	for(auto &he : herd){        									  // Reconstructs other events from recovery times
		he.sample_exp_inf();                                          // Sets exposed and infection times for observed
	
		he.add_unobs();                                               // Adds plausible distribution of unobserved
		he.set_T0();                                                  // Sets initial infection time
	}
	
	mcmc_initialise();												  // For outputing information about the chains
	
	for(auto &he : herd){											  // Outputs some info on the initial state straight away
		auto ninf = 0u, nobs = 0u;
		for(const auto &in : he.ind){
			if(in.inf == true) ninf++;
			if(in.obs == true) nobs++;
		}
			
		cout << "Simulation of herd " << he.index << "     ";
		cout << "  Infected: " << ninf << "      ";
		cout << "  Observed: " << nobs << "      ";
		cout << "  Cull time: " << he.T_cull << endl;
	
		if(true){
			for(const auto &in : he.ind){
				if(in.inf == true){
					cout << "  Exposed: " << in.Exp_t << " ";
					cout << "Infected: " << in.Inf_t << " ";
					cout << "Recovered: " << in.Rec_t << " ";
					cout << "LatentP: " << in.latentP << " ";
					cout << "InfectiousP: " << in.infectiousP << " ";
					if(in.obs == true) cout << "Observed";
					cout << endl;
				}
			}
		}
		cout << endl;
	}

	// Perform mcmc

	trace_init();
	ppc_init();
	
	time_total = -clock();
	for(s = 0; s < nsamp; s++){
		if(s%1000 == 0) cout << "Sample " << s << endl;
		
		if(s%50000 == 0){
			for(auto &he : herd_ppc) he.simulate(int(0.15*herd_sizes[he.index]));
			ppc_plot();
		}
		
		if(s%100 == 0) trace_plot();
	
		for(auto &he : herd){
			if(model == EXT_INF) he.beta1_proposal();              // Proposals to change the external infection rate

			he.R0_proposal();                   // Proposals to change transmission rate
			
			he.multi_proposal();                // Proposals which change exposure / infection times
	
			he.swap_proposal();                 // Proposals which change which individuals are indexes
			
			he.T0_proposal();                   // Proposals which change initial infection time
			
			he.joint_T0_event_proposal();       // Proposals which change T0 and other events
		}
		
		randomwalk_proposal();                                   // Proposals to change parameters for gamma distributed transitions
		
		check(0);                                                // Checks algorithm is performing correctly
	}
	time_total += clock();
	
	if(true){
		for(auto he : herd){
			for(auto i = 0u; i < inf_sampler_bin; i++){
				cout << he.inf_sample[i] << ",";
			}
			cout << "\n";
		}
	}
	
	mcmc_diagnostics();
}	


/// Initialisises mcmc
void mcmc_initialise()
{
	choose_init();
	
	Pr = prior();
	for(auto &he : herd){
		he.Li = he.likelihood();
		
		for(auto &in : he.ind) in.Li_gamma = in.likelihood_gamma();
		
		he.Pr_index = he.prior_index();
	
		he.inf_sample.resize(inf_sampler_bin);
		for(auto &bin : he.inf_sample) bin = 100;
		
		he.nprop = 10;
			
		he.nmulti_propose = 0; he.nmulti_accept = 0; 
		he.nT0_propose = 0; he.nT0_accept = 0;
		
		he.T0_joint_jump = 1;
		he.nT0_joint_propose = 0; he.nT0_joint_accept = 0;
	}
}

/// Initialises a herd with a transmission rate and number of individuals 
void Herd::initialise(int index_)
{
	index = index_;
	
	if(!is_ppc){	// There are duplicate herds for simulating the PPCs.
		if(model == EXT_INF) beta1_param = add_param("beta1_"+to_string(index),UNIFORM_PRIOR,-15.0,0.0,UNSET);  // Lower limit of -15 because lower cases -inf values and numerical issues
		else beta1_param = UNSET;
	
		R0_param = add_param("R0_"+to_string(index),TRUNC_GAMMA_PRIOR,7.4,2.0,UNSET);
	} else {
		beta1_param = param.size()-2;
		R0_param = param.size()-1;
	}
}

/// Simulates from a herd until a certain number of individuals recovers (or epidemic dies out)
void Herd::simulate(const int nrecover)
{
	ind.resize(herd_sizes[index]);

	for(auto &in : ind){
		in.inf = false;
		in.obs = false;
	}
	
	auto beta1 = 0.0; if(model == EXT_INF) beta1 = exp(param[beta1_param].value);
	auto beta2 = param[R0_param].value/param[mu_I].value;

	auto tp = set_trans_param();

	// For simplicity, assume the epidemic starts at time T0
	T0 = 100.0;
	
	switch(model){
		case EXT_INF:
			ind[0].inf = true;
			ind[0].sample_sim(T0,tp);
			break;
			
		case MULTI_INIT:
			auto N = ind.size();                      // Samples the number of index cases
			vector <double> prob_sum(N);
			prob_sum[0] = 0;
			auto sum = 0.0;
			for(auto n = 1u; n < N; n++){
				if(n <= index_max) sum += exp(-alpha*n);
				prob_sum[n] = sum;
			}
			
			auto z = ran()*sum;
			auto n = 0u; while(n < N && z > prob_sum[n]) n++;
			
			cout << "Herd " << index << " number of index cases: " << n << endl;
			
			for(auto j = 0u; j < n; j++){
				ind[j].inf = true;
				ind[j].sample_sim(T0,tp);
			}
			break;
	}
	
	auto t = T0;
	
	vector <int> sus_ind;     // A list of susceptible individuals
	vector <int> inf_ind;     // A list of infected individuals  
		
	for(auto i = 0u; i < ind.size(); i++){
		if(ind[i].inf == false) sus_ind.push_back(i);
		else inf_ind.push_back(i);
	}
		
	auto I = 0u;                       // The number of infectious individuals
	auto N = ind.size();               // The number of alive individuals

	auto N_min = N - nrecover;         // Simulates until a specified number of individuals have recovered (this sets the culling time)
	while(N > N_min){
		auto tnext = LARGE;              // This works out when the next event is
		int tnexti = UNSET;
		for(auto i : inf_ind){
			if(ind[i].inf == true){ 
				auto tt = ind[i].Inf_t;
				if(tt > t && tt < tnext){ tnext = tt; tnexti = i;}
			
				tt = ind[i].Rec_t;
				if(tt > t && tt < tnext){ tnext = tt; tnexti = i;}
			}
		}

		if(I == 0 && tnext == LARGE) break;
		
		auto rate = sus_ind.size()*(beta1 + I*beta2/N);
		
		auto tnew = t + exp_sample(rate);
		
		if(tnew < tnext){                 // Expose a new individual
			t = tnew;
			
			auto j = int(ran()*sus_ind.size());
			auto i = sus_ind[j];
			sus_ind.erase(sus_ind.begin()+j);
			inf_ind.push_back(i);
			
			ind[i].inf = true;
			
			
			ind[i].sample_sim(t,tp);
		}
		else{
			t = tnext;
			
			if(tnext == ind[tnexti].Inf_t){ // An individuals becomes infectious
				I++;	
			}
			else{                           // An individual recovers
				if(tnext != ind[tnexti].Rec_t) emsg("Problem 1");
				I--; N--;
			}
		}
	}
	
	T_cull = t+0.00001;
	for(auto &in : ind){           // Sets individuals which recover before T_cull to be observed
		if(in.inf == true && in.Rec_t <= T_cull) in.obs = true;
		else in.obs = false;
	}	
}


/// This loads up the data into the herd
void Herd::load_data()
{
	ind.resize(herd_sizes[index]);
	
	for(auto &i : ind){
		i.inf = false;
		i.obs = false;
	}
	
	cout << remtimes[index].size() << " size\n";
	
	T_cull = -LARGE;
	for(auto ind_index = 0u; ind_index < remtimes[index].size(); ind_index++){
		auto Rt = remtimes[index][ind_index]+115.0;	  // Adding 115 here so that timings are broadly comparable with simulated data
		
		ind[ind_index].Rec_t = Rt;
		ind[ind_index].inf = true;
		ind[ind_index].obs = true;
		
		if(Rt > T_cull) T_cull = Rt; 
	}
	
	T_cull += TINY;                                 // Cull happens just after the last recovery
}

/// Adds exposure and infection times to the known removal times
void Herd::sample_exp_inf()
{
	auto tp = set_trans_param();
	
	for(auto &in : ind){
		if(in.obs == true){
			in.sample_backwards(tp);
			if(in.Exp_t < 0.0) emsg("Problem");
		}
	}
}

/// This whole chunk below simulates an outbreak in a herd of the same size and uses this to propose unobserved infections
void Herd::add_unobs()
{
	auto nobs = 0u;
	for(auto &in : ind){                     // First we turn off any infected, unobserved individuals (fro simulation)
		if(in.obs == false){
			in.inf = false;
		}
		else nobs++;
	}
	
	Herd herd_sim;                           // Simulates the same number of recoveries 
	herd_sim.index = index;
	herd_sim.beta1_param = beta1_param;
	herd_sim.R0_param = R0_param;
	
	herd_sim.simulate(nobs);
	
	auto time_shift = T_cull - herd_sim.T_cull;  // Shifts herd and herd)sim so they have the same cull time
	for(auto &ind_sim : herd_sim.ind){
		if(ind_sim.inf == true){
			ind_sim.Exp_t += time_shift;
			ind_sim.Inf_t += time_shift;
			ind_sim.Rec_t += time_shift;
		}
	}

	for(const auto &ind_sim : herd_sim.ind){// Copies over infected unobserved individuals
		if(ind_sim.inf == true && ind_sim.obs == false){
			for(auto &in : ind){
				if(in.obs == false && in.inf == false){
					in = ind_sim;
					if(in.Exp_t < 0) emsg("Problem out of range");
					if(in.Rec_t <= T_cull) emsg("unobs problem"); 
					break;
				}
			}
		}
	}
}

/// Sets the time of the initial infection(s)
void Herd::set_T0()
{
	T0 = LARGE;
	for(const auto &in : ind){
		if(in.inf == true){
			if(in.Exp_t < T0) T0 = in.Exp_t;
		}
	}
}

 
/// Performs a random walk MH proposal for beta1
void Herd::beta1_proposal()
{
	auto &par = param[beta1_param];

	auto param_store = par.value;
	
	par.value += normal_sample(0,par.jump);
	
	auto L_prop = likelihood();                 // Calculate the new likelihood and prior
	auto Pr_prop = prior();
	
	auto al = exp(phi_L*(L_prop  - Li) + Pr_prop - Pr);  // Calculates the Metropolis-Hastings acceptance probability

	par.npropose++;
	if(ran() < al){
		par.naccept++;
		Li = L_prop;
		Pr = Pr_prop;
		if(s < nburnin) par.jump *= 1.001;          // The dynamically adapts the jumping size
	}
	else{
		par.value = param_store;
		if(s < nburnin) par.jump *= 0.9995;
	}
}

/// Performs a random walk MH proposal for R0 parameter
void Herd::R0_proposal()
{
	auto &par = param[R0_param];

	auto param_store = par.value;
	
	par.value += normal_sample(0,par.jump);
	
	auto L_prop = likelihood();                 // Calculate the new likelihood and prior
	auto Pr_prop = prior();
	
	if(par.value < 0) L_prop = -LARGE;
	
	auto al = exp(phi_L*(L_prop - Li) + Pr_prop - Pr);  // Calculates the Metropolis-Hastings acceptance probability

	par.npropose++;
	if(ran() < al){
		par.naccept++;
		Li = L_prop;
		Pr = Pr_prop;
		if(s < nburnin) par.jump *= 1.001;          // The dynamically adapts the jumping size
	}
	else{
		par.value = param_store;
		if(s < nburnin) par.jump *= 0.9995;
	}
}

/// This performs random walk proposals for the four parameters associates with latent and infectious times
void randomwalk_proposal()
{
	for(auto loop = 0u; loop < 4; loop++){
		int p;
		switch(loop%4){
			case 0: p = mu_L; break;
			case 1: p = mu_I; break;
			case 2: p = sh_L; break;
			case 3: p = sh_I; break;
		}

		auto &par = param[p];

		auto herd_store = herd;
		auto param_store = par.value;
		
		par.value += normal_sample(0,par.jump);
		
		auto dLi = 0.0, dLi_gamma = 0.0;
		for(auto &he : herd){
			if(p == mu_I){
				auto Li_prop = he.likelihood();
				dLi += Li_prop-he.Li;
				he.Li = Li_prop;
			}
			
			for(auto &in : he.ind){
				auto Li_gamma_prop = in.likelihood_gamma();
				dLi_gamma += Li_gamma_prop - in.Li_gamma;
				in.Li_gamma = Li_gamma_prop;
			}
		}
		
		auto Pr_prop = prior();
			
		auto al = exp(phi_L*dLi + phi_G*dLi_gamma + Pr_prop - Pr); 
	
		par.npropose++;
		if(ran() < al){
			par.naccept++;
			Pr = Pr_prop;
			if(s < nburnin) par.jump *= 1.001;                // The dynamically adapts the jumping size
		}
		else{
			par.value = param_store;
			herd = herd_store;
			if(s < nburnin) par.jump *= 0.9995;
		}
	}
}



/// Performs proposals which change T0 
void Herd::T0_proposal()
{		
	auto tp = set_trans_param();

	vector <int> list_index;               

	auto time_max = LARGE;                                   // Finds the maximum time for T0
	for(auto i = 0u; i < ind.size(); i++){
		const auto &in = ind[i]; 
		if(in.inf == true){
			if(in.Exp_t == T0){
				list_index.push_back(i);
				if(in.Inf_t < time_max) time_max = in.Inf_t;
			}
			else{
				if(in.Exp_t < time_max) time_max = in.Exp_t;
			}
		}
	}

	auto var = tp.mu_L_sample*tp.mu_L_sample/(list_index.size()*tp.sh_L_sample);

	auto T0_store = T0;

	T0 += normal_sample(0,sqrt(var)); 

	if(T0 < 0 || T0 >= time_max){
		T0 = T0_store;
	}
	else{
		auto beta1 = 0.0; if(model == EXT_INF) beta1 = exp(param[beta1_param].value);
		auto dLi = beta1*(T0-T0_store)*(ind.size()-list_index.size()); // DE edit 05 Nov to make this affect everything other than index individuals 
		
		auto dLi_gamma = 0.0;
		for(auto i : list_index){
			auto &in = ind[i];
			in.Exp_t = T0;	
			in.latentP = in.Inf_t-in.Exp_t;
			
			in.Li_gamma_store = in.Li_gamma;
			in.Li_gamma = in.likelihood_gamma(); 
			dLi_gamma += in.Li_gamma - in.Li_gamma_store;
		}
		
		
		auto al = exp(phi_L*dLi + phi_G*dLi_gamma);  // Calculates the Metropolis-Hastings acceptance probability

		nT0_propose++;
		if(ran() < al){
			Li += dLi;
			nT0_accept++;
		}
		else{
			T0 = T0_store;
			for(auto i : list_index){
				auto &in = ind[i];
				in.Exp_t = T0;
				in.latentP = in.Inf_t-in.Exp_t;
				in.Li_gamma = in.Li_gamma_store;
			}
		}
	}
}


/// The performs proposals which swap individuals between being index / non-index cases
void Herd::swap_proposal()
{
	// First we make a list 
	vector <Event> event;
	auto num_index = 0u;
	for(const auto &in : ind){
		if(in.inf == true){
			if(in.Exp_t == T0) num_index++;
					
			if(in.Inf_t <= T_cull){
				Event ev; ev.type = INFECT; ev.t = in.Inf_t;
				event.push_back(ev);
			}
			
			if(in.Rec_t <= T_cull){
				Event ev; ev.type = RECOVERY; ev.t = in.Rec_t;
				event.push_back(ev);
			}
		}
	}
	
	sort(event.begin(),event.end(),EV_ord);
	Event ev; ev.type = TERMINATE; ev.t = T_cull;
	event.push_back(ev);
	
	auto nevent = event.size();

	auto I = 0;             // The number of infected
	auto S = ind.size();    // The number of susceptible
	auto N = S;             // The number alive
	vector <int> Nstore(nevent);
	vector <int> Istore(nevent);
	
	for(auto e = 0u; e < nevent; e++){
		Nstore[e] = N; Istore[e] = I;
		switch(event[e].type){
		case EXPOSE: break;
		case INFECT: I++; break;
		case RECOVERY: I--; N--; break;
		case TERMINATE: break;
		}
	}
	
	auto tp = set_trans_param();
	
	auto beta1 = 0.0; if(model == EXT_INF) beta1 = exp(param[beta1_param].value);
	auto beta2 = param[R0_param].value/param[mu_I].value;
	
	auto Exp_t_max = T0 + 2*tp.mu_L_value;     // This sets the limit for the exposure time
	
	for(auto loop = 0u; loop < 3; loop++){
		for(auto &in : ind){
			if(in.inf == true){
				if(in.Exp_t > T0){    // Add an index case
					if(in.Exp_t < Exp_t_max && num_index < index_max){
						auto t = in.Exp_t;
					
						auto dLi = 0.0;
					
						auto e = 0u; 
						auto tt = T0;
						while(e < nevent && t > event[e].t){
							auto tnew = event[e].t;
							auto rate = beta1 + beta2*Istore[e]/Nstore[e];
							dLi += (tnew-tt)*rate; tt = tnew;
							e++;
						}
						if(e == nevent) emsg("prob");
						auto rate = beta1 + beta2*Istore[e]/Nstore[e];
						dLi += (t-tt)*rate - choose(ind.size(),num_index+1) + choose(ind.size(),num_index);
						if(rate > 0) dLi -= log(rate); else dLi -= NO_INF;
						
						auto dPr = -alpha;
						
						auto probif = 0;
						auto probfi = gamma_probability(in.latentP,tp.mu_L_sample,tp.sh_L_sample);
						
						auto latP_st = in.latentP;
						in.latentP = in.Inf_t-T0;
						auto Li_gamma_prop = in.likelihood_gamma();
						
						auto al = exp(phi_L*dLi + phi_G*(Li_gamma_prop - in.Li_gamma) + dPr + probfi - probif);
						if(ran() < al){
							Li += dLi;
							in.Li_gamma = Li_gamma_prop;
							Pr_index += dPr;
							in.Exp_t = T0;
							num_index++;
						}
						else{
							in.latentP = latP_st;
						}
					}
				}
				else{                 // Remove an index
					if(num_index > 1){
						auto latP = gamma_sample(tp.mu_L_sample,tp.sh_L_sample);
						auto t = in.Inf_t-latP;
						if(t > T0 && t < T_cull && t < Exp_t_max){
							auto dLi = 0.0;
					
							auto e = 0u; 
							auto tt = T0;
							while(e < nevent && t > event[e].t){
								auto tnew = event[e].t;
								auto rate = beta1 + beta2*Istore[e]/Nstore[e];
								dLi -= (tnew-tt)*rate; tt = tnew;
								e++;
							}
							if(e == nevent) emsg("prob");
							auto rate = beta1 + beta2*Istore[e]/Nstore[e];
							dLi += -(t-tt)*rate   - choose(ind.size(),num_index-1) + choose(ind.size(),num_index);
							if(rate > 0) dLi += log(rate); else dLi += NO_INF;
							
							
							auto dPr = alpha;

							auto latP_st = in.latentP;
							in.latentP = latP;
							auto Li_gamma_prop = in.likelihood_gamma();
							
							auto probif = gamma_probability(latP,tp.mu_L_sample,tp.sh_L_sample);
							auto probfi = 0;
							
							auto al = exp(phi_L*dLi + phi_G*(Li_gamma_prop-in.Li_gamma) + dPr + probfi - probif);
							if(ran() < al){
								Li += dLi;
								in.Li_gamma = Li_gamma_prop;
								Pr_index += dPr;
								
								in.Exp_t = t;		
								num_index--;
							}
							else{
								in.latentP = latP_st;
							}
						}	
					}
				}
			}
		}
	}
}

/// Makes multiple changes to events (the initial infection time is not changed)
void Herd::multi_proposal()
{
	const double prob_add_rem = 0.4;                   // This gives the probability of doing an add/rem proposal
	const double prob_obs = 0.2;                       // Prob of changing events on observed individual
	const double prob_unobs = 1-prob_add_rem-prob_obs; // Prob of changing events on unobserved individual
	
	if(prob_unobs < 0) emsg("Samping error");
	
	// This generates a sampler for adding new exposure times
	vector <double> prob(inf_sampler_bin), prob_sum(inf_sampler_bin);

	auto sum = 0.0; for(auto i = 0u; i < inf_sampler_bin; i++) sum += inf_sample[i];
	
	auto sum2 = 0.0;
	for(auto i = 0u; i < inf_sampler_bin; i++){
		prob[i] = inf_sample[i]/sum;
		sum2 += prob[i];
		prob_sum[i] = sum2;
	}

	auto tp = set_trans_param();
	
	for(auto loop = 0u; loop < nupd; loop++){
		auto probif = 0.0, probfi = 0.0;
		
		// Makes lists of observed, infected (unobseved) and susceptible 
		vector <int> list_obs, list_inf, list_sus;
		auto num_index = 0u;
		for(auto i = 0u; i < ind.size(); i++){
			if(ind[i].obs == false){
				if(ind[i].inf == false) list_sus.push_back(i);
				else list_inf.push_back(i);
			}
			else list_obs.push_back(i);
			
			if(ind[i].inf == true && ind[i].Exp_t == T0) num_index++;
		}
		
		auto nlist_inf = list_inf.size();
		auto nlist_sus = list_sus.size();
		auto nlist_obs = list_obs.size();
	
		vector <Individual> ind_store = ind;
	
		auto jmax = (unsigned int)(nprop+0.5); if(jmax == 0) jmax = 1;

		auto dLi_gamma = 0.0;
		for(auto j = 0u; j < jmax; j++){
			
			PropType type;
			
			auto z = ran();                                // Randomly select what change to make
			if(z < prob_add_rem){                          // Adds or removes an individual
				if(ran() < 0.5) type = ADD_UNOBS;
				else type = REM_UNOBS;
			}
			else{
				z -= prob_add_rem;
				if(z < prob_obs) type = RESAMPLE_OBS;        // Resamples times for observed
				else type = RESAMPLE_UNOBS;                  // Resamples  times for unobserved
			}
		
			switch(type){
				case ADD_UNOBS:                              // Adds a new infected individual
					if(nlist_sus > 0){
						auto sus_sel = int(ran()*nlist_sus);     // Randomly selects a susceptible individual
						auto i = list_sus[sus_sel];         
						
						auto z = ran();                          // Selects which bin to sample from
						auto bin = 0u; while(bin < inf_sampler_bin && z > prob_sum[bin]) bin++;
						if(bin == inf_sampler_bin) emsg("Problem 2");
						
						auto t_exp = (bin+ran())*T_cull/inf_sampler_bin; // Samples exposure time 
						if(t_exp > T0){
							auto &in = ind[i];
							in.sample(t_exp,tp);
							
							if(in.Rec_t > T_cull){               // Only accept if recovery time after cull
								in.inf = true;                     // Works out proposal probability (and reverse)
								probif += in.sample_probability(tp);
								probif += log((1.0/nlist_sus)*prob[bin]*(inf_sampler_bin/T_cull));
								probfi += log((1.0/(nlist_inf+1)));
								
								list_inf.push_back(i); nlist_inf++;	   // Updates the lists
								list_sus.erase(list_sus.begin()+sus_sel); nlist_sus--;	
								
								auto Li_gamma_prop = in.likelihood_gamma();
								dLi_gamma += Li_gamma_prop - in.Li_gamma;
								in.Li_gamma = Li_gamma_prop;
							}	
						}
					}
					break;
					
				case REM_UNOBS:                              // Removes an existing infected
					if(nlist_inf > 0){
						auto inf_sel = int(ran()*nlist_inf);     // Randomly selects a infected individual
						auto i = list_inf[inf_sel];       
						
						auto t_exp = ind[i].Exp_t;
						if(t_exp > T0){                           // Only accept if not an index case
							auto &in = ind[i];
							in.inf = false;
							
							auto bin = int(inf_sampler_bin*t_exp/T_cull);
							
							probif += log((1.0/nlist_inf));
							probfi += log((1.0/(nlist_sus+1))*prob[bin]*(inf_sampler_bin/T_cull));
							probfi += in.sample_probability(tp);
								
							list_sus.push_back(i); nlist_sus++;	
							list_inf.erase(list_inf.begin()+inf_sel); nlist_inf--;	
							
							auto Li_gamma_prop = in.likelihood_gamma();
							dLi_gamma += Li_gamma_prop - in.Li_gamma;
							in.Li_gamma = Li_gamma_prop;
						}
					}
					break;
				
				case RESAMPLE_OBS:                          // Makes a change to an observed infected individual
					if(nlist_obs > 0){
						auto &in = ind[list_obs[int(ran()*nlist_obs)]];
						if(in.Exp_t == T0){                     // Resamples infection time for index case
							auto mu1 = in.Exp_t + tp.mu_L_sample;
							auto var1 = tp.mu_L_sample*tp.mu_L_sample/tp.sh_L_sample;
							
							auto mu2 = in.Rec_t - tp.mu_I_sample;
							auto var2 = tp.mu_I_sample*tp.mu_I_sample/tp.sh_I_sample;
							
							auto mu = (mu1*var2 + mu2*var1)/(var1+var2);
							auto var = var1*var2/(var1+var2);
						
							auto Inf_t = normal_sample(mu,sqrt(var));
							if(Inf_t > in.Exp_t && Inf_t < in.Rec_t){
								probfi += normal_probability(in.Inf_t,mu,var);
								probif += normal_probability(Inf_t,mu,var);
								in.Inf_t = Inf_t;
								in.latentP = in.Inf_t - in.Exp_t;
								in.infectiousP = in.Rec_t - in.Inf_t;	
								auto Li_gamma_prop = in.likelihood_gamma();
								dLi_gamma += Li_gamma_prop - in.Li_gamma;
								in.Li_gamma = Li_gamma_prop;								
							}
						}
						else{
							Individual ind_store = in;
							auto dprobfi = in.sample_probability(tp);
							in.sample_backwards(tp);
							if(in.Exp_t <= T0) in = ind_store;		
							else{
								probfi += dprobfi;
								probif += in.sample_probability(tp); 
								auto Li_gamma_prop = in.likelihood_gamma();
								dLi_gamma += Li_gamma_prop - in.Li_gamma;
								in.Li_gamma = Li_gamma_prop;
							}
						}
					}
					break;
					
				case RESAMPLE_UNOBS:                        // Makes a change to an unobserved infected individual
					if(nlist_inf > 0){
						auto &in = ind[list_inf[int(ran()*nlist_inf)]];
						auto dprobfi = in.sample_probability(tp);
							
						Individual ind_store = in;
						if(ran() < 0.5){
							if(in.Exp_t != T0){	
								in.sample_backwards(tp);
								if(in.Exp_t >= T_cull || in.Exp_t <= T0) in = ind_store; 
								else{
									probfi += dprobfi;
									probif += in.sample_probability(tp); 
									auto Li_gamma_prop = in.likelihood_gamma();
									dLi_gamma += Li_gamma_prop - in.Li_gamma;
									in.Li_gamma = Li_gamma_prop;
								}
							}
						}
						else{
							in.sample(in.Exp_t,tp);
							if(in.Rec_t < T_cull) in = ind_store;
							else{
								probfi += dprobfi;
								probif += in.sample_probability(tp); 
								auto Li_gamma_prop = in.likelihood_gamma();
								dLi_gamma += Li_gamma_prop - in.Li_gamma;
								in.Li_gamma = Li_gamma_prop;
							}
						}
					}
					break;
			}
		}

		auto Li_prop = likelihood();  
		
		// Calculates the Metropolis-Hastings acceptance probability
		auto al = exp(phi_L*(Li_prop - Li) + phi_G*dLi_gamma + probfi - probif);  

		nmulti_propose++;
		if(ran() < al){
			nmulti_accept++;
			Li = Li_prop;
			if(s < nburnin){ nprop *= 1.001; if(nprop > 100) nprop = 100;}
		}
		else{
			ind = ind_store;
			if(s < nburnin){ nprop *= 0.999; if(nprop < 1) nprop = 1;}
		}
	}	
	
	if(s < nburnin){                                    // This dynamically adapts the exposure time sampler
		for(auto i = 0u; i < ind.size(); i++){
			auto t_exp = ind[i].Exp_t;
			if(ind[i].inf == true && ind[i].obs == false && t_exp > 0){
				auto bin = int(inf_sampler_bin*t_exp/T_cull);
				if(bin < 0 || bin >= inf_sampler_bin) emsg("Problem 3");
				inf_sample[bin]++;
			}
		}
	}
}

/// Sets the parameter values
void param_prior_init()
{
	for(auto &par : param){
		if(par.name == "mu_I"){
			par.value = 9.0;
		} else if (par.name == "mu_L"){
			par.value = 7.0;
		} else if (par.name == "sh_I"){
			par.value = 22.0;
		} else if (par.name == "sh_L"){
			par.value = 16.0;
		} else if (par.name == "beta1_0"){
			par.value = log(0.01/herd_sizes[0]);
		} else if (par.name == "R0_0"){
			par.value = 8.0;
		}
		cout << "Sampled value for parameter " << par.name << ": " << par.value << "\n";			
	}
}

/// Samples the parameters from the model prior
void param_prior_sample()
{
	for(auto &par : param){
		switch(par.prior_type){
		case UNIFORM_PRIOR:
			par.value = par.prior_val1 + ran()*(par.prior_val2-par.prior_val1);
			break;
			
		case BETA_PRIOR:
			par.value = beta_sample(par.prior_val1,par.prior_val2)*par.upper_bound;
			break;
			
		case GAMMA_PRIOR: 
			par.value = gamma_sample(par.prior_val1,par.prior_val2);
			break;
			
		case TRUNC_GAMMA_PRIOR:
			par.value = trunc_gamma_sample(par.prior_val1,par.prior_val2,par.name);
			break;
			
		case EXP_PRIOR:
			par.value = 1.0/exp_sample(par.prior_val1);
			break;
		}
		cout << "Sampled value for parameter " << par.name << ": " << par.value << "\n";			
	}
}

/// Calculates the prior
double prior()
{
	auto Pr = 0.0;
	for(auto &par : param){
		switch(par.prior_type){
			case UNIFORM_PRIOR:
				if(par.value < par.prior_val1 || par.value > par.prior_val2) Pr -= LARGE; 
				else Pr += log(1.0/(par.prior_val2-par.prior_val1));
				break;
				
			case GAMMA_PRIOR:
				Pr += gamma_probability(par.value, par.prior_val1, par.prior_val2); 
				break;
				
			case TRUNC_GAMMA_PRIOR:
				Pr += trunc_gamma_probability(par.value, par.prior_val1, par.prior_val2, par.name);
				break;
				
			case EXP_PRIOR: 
				Pr += exp_probability(par.value, par.prior_val1); 
				break;	
				
			case BETA_PRIOR: 
				if(par.value > par.upper_bound) Pr -= LARGE;
				else Pr += beta_probability(par.value, par.prior_val1, par.prior_val2)/par.upper_bound;
				break;
		}
	}
	
	return Pr;
}
 

/// Calculate the log likelihood associated with exposure events
double Herd::likelihood()
{
	if(map_out_prior == true) return 0;
	time_likelihood -= clock();
	
	// Gererates a sorted list of event for the entire herd
	
	vector <Event> event;
	for(const auto &in : ind){
		if(in.inf == true){
			Event ev; ev.type = EXPOSE; ev.t = in.Exp_t;
			event.push_back(ev);
			
			if(in.Inf_t <= T_cull){
				Event ev; ev.type = INFECT; ev.t = in.Inf_t;
				event.push_back(ev);
			}
			
			if(in.Rec_t <= T_cull){
				Event ev; ev.type = RECOVERY; ev.t = in.Rec_t;
				event.push_back(ev);
			}
		}
	}
	
	Event ev; ev.type = TERMINATE; ev.t = T_cull;
	event.push_back(ev);
	
	sort(event.begin(),event.end(),EV_ord);

	if(T0 != event[0].t){ cout << T0 << " " << event[0].t << " j\n"; emsg(" Initial event times do not agree");}
	
	// Calculates the log likelihood
	
	auto beta1 = 0.0; if(model == EXT_INF) beta1 = exp(param[beta1_param].value);
	auto beta2 = param[R0_param].value/param[mu_I].value;
		
	auto L = 0.0;           // The log of the likelihood
	
	auto t = 0.0;           // The start time
	
	auto I = 0;             // The number of infected
	auto S = ind.size();    // The number of susceptible
	auto N = S;             // The number alive
	
	auto num_index = 0u;  
	for(auto i = 0u; i < event.size(); i++){
		if(N == 0) break;
		
		auto tt = event[i].t;

		auto rate = 0.0;
		if(i != 0) rate = beta1 + beta2*I/N;
			
		L -= (tt-t)*rate*S;
			
		if(tt == T0) num_index++;
		else{
			if(event[i].type == EXPOSE){
				if(rate == 0.0) L += NO_INF; // Make an exposure with no infected individuals very unlikely
				else L += log(rate);
			}					
		}
		
		switch(event[i].type){
			case EXPOSE: S--; break;
			case INFECT: I++; break;
			case RECOVERY: I--; N--;  if(I < 0) emsg("Negative");  if(N < 0) emsg("Negative"); break;
			case TERMINATE: break;
		}
		t = tt;
	}

	L -= choose(ind.size(),num_index);   // This implements the likelihood of distribution in the index cases at T0
	
	time_likelihood += clock();
	return L;
}

/// Calculate the likelihood for the 
double Individual::likelihood_gamma()
{
	if(inf == false) return 0;

	return gamma_probability(latentP,param[mu_L].value,param[sh_L].value) +
       	 gamma_probability(infectiousP,param[mu_I].value,param[sh_I].value);
}

/// Calculates the prior for the number of index cases
double Herd::prior_index()
{
	auto num_index = 0u;                   
	for(const auto &in : ind){
		if(in.Exp_t == T0) num_index++;
	}
	return -alpha*num_index;
}

/// Sets transition parameters
TransParam set_trans_param()
{
	TransParam tp;
	
	tp.mu_L_value = param[mu_L].value;
	tp.sh_L_value = param[sh_L].value;
	tp.mu_I_value = param[mu_I].value;
	tp.sh_I_value = param[sh_I].value;
	
	tp.sh_L_sample = 1+phi_G*(tp.sh_L_value-1);
	tp.mu_L_sample = tp.mu_L_value*tp.sh_L_sample/(phi_G*tp.sh_L_value);
	tp.sh_I_sample = 1+phi_G*(tp.sh_I_value-1);
	tp.mu_I_sample = tp.mu_I_value*tp.sh_I_sample/(phi_G*tp.sh_I_value);

	return tp;
}



/// Based on a specified exposure time this samples the individual's time line
void Individual::sample(const double expt, const TransParam &tp)
{
	latentP = gamma_sample(tp.mu_L_sample, tp.sh_L_sample);
	infectiousP = gamma_sample(tp.mu_I_sample,tp.sh_I_sample);
	
	Exp_t = expt; 
	Inf_t = Exp_t + latentP;
	Rec_t = Inf_t + infectiousP;
}

/// Very slightly modified version of the above which allows me to pass different values
void Individual::sample_sim(const double expt, const TransParam &tp)
{
	latentP = gamma_sample(tp.mu_L_value, tp.sh_L_value);
	infectiousP = gamma_sample(tp.mu_I_value,tp.sh_I_value);
	
	Exp_t = expt; 
	Inf_t = Exp_t + latentP;
	Rec_t = Inf_t + infectiousP;
}

/// Gets the probability of sampling a certain set of latent and infectious periods
double Individual::sample_probability(const TransParam &tp)
{
	return gamma_probability(latentP,tp.mu_L_sample,tp.sh_L_sample) +
       	 gamma_probability(infectiousP,tp.mu_I_sample,tp.sh_I_sample);
}

/// Based on a recovery time the exposure time and infecttous times are sampled
void Individual::sample_backwards(const TransParam &tp)
{
	latentP = gamma_sample(tp.mu_L_sample, tp.sh_L_sample);
	infectiousP = gamma_sample(tp.mu_I_sample,tp.sh_I_sample);

	Inf_t = Rec_t - infectiousP;
	Exp_t = Inf_t - latentP;
}


/// Adds a parameter to the list of parameters
int add_param(string name, PriorType prior_type, double prior_val1, double prior_val2, double upper_bound)
{
	auto p = param.size();

	double mean = UNSET;
	switch(prior_type){
		case UNIFORM_PRIOR: mean = (prior_val1+prior_val2)/2; break;
		case BETA_PRIOR: mean = upper_bound*prior_val1/(prior_val1+prior_val2); break;
		case GAMMA_PRIOR: mean = prior_val1; break;
		case TRUNC_GAMMA_PRIOR: mean = prior_val1; break;
		case EXP_PRIOR: mean = prior_val1; break;
	}
		
	Param par;
	par.name = name;
	par.value = UNSET;
	par.prior_type = prior_type; par.prior_val1 = prior_val1; par.prior_val2 = prior_val2; 
	par.upper_bound = upper_bound; 
	if(mean >= 0){
		par.jump = mean/10; 
	} else {
		par.jump = -mean/50;
	}
	par.npropose = 0; par.naccept = 0;
	par.jump_MBP = mean/50;
	par.npropose_MBP = 0; par.naccept_MBP = 0; par.nfail_MBP = 0;
	param.push_back(par);
	
	return p;
}

/// Initialises trace plots
void trace_init()
{
	trace << "State"; 
	for(auto &pa : param) trace << "\t" << pa.name; 
	for(auto &he : herd) trace << "\tNinf Herd" << he.index;
	for(auto &he : herd) trace << "\tTinit Herd" << he.index;
	for(auto &he : herd) trace << "\tInit.infected" << he.index;
	for(auto &he : herd) trace << "\tExt.infected" << he.index;	
	for(auto &he : herd) trace << "\tLi Herd" << he.index;
	for(auto &he : herd) trace << "\tPr index Herd" << he.index;
	trace << "\tPrior";
	trace << endl;	
}

/// Initialises ppc plots
void ppc_init()
{
	ppc << "State";
	ppc << "\tHerd";
	ppc << "\tExposure";
	ppc << "\tInfection";
	ppc << "\tRemoval";
	ppc << "\tObserved";
	ppc << endl;	
}


/// Plots useful quantities
void trace_plot()
{
	trace << s;
	for(auto &pa : param) trace << "\t" << pa.value;
	
	for(auto &he : herd){
		auto ninf = 0u; for(const auto &in: he.ind){ if(in.inf == true) ninf++;}
		trace << "\t" << ninf;
	}
	
	for(auto &he : herd) trace << "\t" << he.T0;
	
	for(auto &he : herd){
		auto num = 0u;  for(const auto &in: he.ind){ if(in.inf == true && in.Exp_t == he.T0) num++;}
		trace << "\t" << num;
	}
	
	for(auto &he : herd){
		auto num = 0u;
		auto min_inf = LARGE;
		for(const auto &in: he.ind){
			if(in.inf == true && in.Inf_t < min_inf) min_inf = in.Inf_t;
		}
		for(const auto &in: he.ind){
			if(in.inf == true && in.Exp_t < min_inf) num++;
		}
		trace << "\t" << num;
	}
	
	for(auto &he : herd) trace << "\t" << he.Li;
	
	for(auto &he : herd) trace << "\t" << he.Pr_index;
		
	trace << "\t" << Pr;
	
	trace << endl;
}

/// Plots useful quantities
void ppc_plot()
{
	for(auto &he : herd_ppc){
		for(const auto &in: he.ind){
			if(in.inf == true){
				ppc << s << "\t" << he.index << "\t" << in.Exp_t << "\t" << in.Inf_t << "\t" << in.Rec_t << "\t" << in.obs;
				ppc << endl;
			}
		}
	}
}

/// Plots useful quantities
void herd_write()
{
	for(auto &he : herd){
		for(const auto &in: he.ind){
			if(in.inf == true & in.obs == true){
				herd_data << in.Rec_t << "\t" << he.index;
				herd_data << endl;
			}
		}
	}
}

/// Initialises trace plots
void herd_init()
{
	herd_data << "Removal Time";
	herd_data << "\tHerd";
	herd_data << endl;	
}

/// Outputs information about how well the chain is mixing
void mcmc_diagnostics()
{
	for(auto &pa : param){
		cout << pa.name << " ";
		cout << int(100*pa.naccept/(pa.npropose+0.001)) << "% random walk acceptace probability   ";
		cout << "Jump size: " << pa.jump << "    "; 
		if(pa.npropose_MBP){
			cout << int(100*pa.naccept_MBP/(pa.npropose_MBP+0.001)) << "% MBP   ";
			cout << int(100*pa.nfail_MBP/(pa.npropose_MBP+0.001)) << "% fail MBP   ";
			cout << "MBP Jump size: " << pa.jump_MBP << "    "; 
		}
		
		cout << endl;
	}
	
	for(auto &he : herd){
		cout << "Herd " << he.index << ": ";
		cout << int(100*he.nmulti_accept/(he.nmulti_propose+0.001)) << "% multi AP       ";
		cout << int(100*he.nT0_accept/(he.nT0_propose+0.001)) << "% T0 AP        ";
		cout << int(100*he.nT0_joint_accept/(he.nT0_joint_propose+0.001)) << "% T0 joint AP   size:" << he.T0_joint_jump << "      ";
		cout << he.nprop << " Number of proposals          ";
		cout << endl;
	}
	
	cout << time_total/(60*CLOCKS_PER_SEC) << " Time in minutes" << endl;
	cout << int(100*time_likelihood/(time_total)) << " Percentage in likelihood" << endl;
}

/// Draws a sample from the uniform distribution
double uniform_sample(const double lwr, const double upr)
{
	uniform_real_distribution<double> distribution(lwr,upr);
	return distribution(generator);
}

/// Draws a sample from the gamma distribution with mean and shape
double gamma_sample(const double mu, const double sh)
{
	gamma_distribution<double> distribution(sh,mu/sh);
	return distribution(generator);
}

/// Draws a sample from the gamma distribution with mean and shape.  The name argument allows for truncation in different places dependent on the parameter
double trunc_gamma_sample(const double mu, const double sh, const string nm)
{
	gamma_distribution<double> distribution(sh,mu/sh);
	auto valid = false;
	auto value = 0.0;
	while(valid == false){
		value = distribution(generator);
		if(nm != "sh_L" && nm != "sh_I" && nm != "mu_L" && nm != "mu_I" && value < 20.0) valid = true;
		if((nm == "sh_L" || nm == "sh_I") && value > 5.0) valid = true;
		if(nm == "mu_L" && value < 10.0) valid = true;
		if(nm == "mu_I" && value < 15.0) valid = true;
	}
	return value;
}

/// Draws a sample from the beta distribution 
double beta_sample(const double a, const double b)
{
	auto X1 = gamma_sample(a,a);
	auto X2 = gamma_sample(b,b);
	return X1/(X1+X2);
}

/// The log of the probability from the gamma distribution
double gamma_probability(const double x, const double mu, const double sh)
{
	auto b = sh/mu;
	return (sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh);
}

/// The log of the probability from the gamma distribution.  The nm parameter specifies how to deal with the truncation for each parameter
double trunc_gamma_probability(const double x, const double mu, const double sh, const string nm)
{
	auto b = sh/mu;
	auto prob = -LARGE;
	if(x > 5.0 && nm == "sh_I") prob = ((sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh))/0.9940157;
	if(x > 5.0 && nm == "sh_L") prob = ((sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh))/0.9896794;
	if(x < 15.0 && nm == "mu_L") prob = ((sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh))/0.9995746;
	if(x < 15.0 && nm == "mu_I") prob = ((sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh))/0.9653511;	
	if(x < 30.0 && nm != "sh_L" && nm != "sh_I" && nm != "mu_L" && nm != "mu_I" && nm != "T0") prob = ((sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh))/0.9972577;
	if(x < 30.0 && nm == "T0") prob = ((sh-1)*log(x) - b*x + sh*log(b) - lgamma(sh))/0.9658502;
	
	return prob;
}

/// The log of the probability from the beta distribution
double beta_probability(const double x, const double a, const double b)
{
  return (a-1)*log(x) + (b-1)*log(1-x) - lgamma(a) - lgamma(b) + lgamma(a+b);
}

/// Samples from the exponential distribution using a rate
double exp_sample(const double rate)
{
	exponential_distribution<double> distribution(rate);
	return distribution(generator);
}

/// The log of the probability from the exponential distribution
double exp_probability(const double x, const double mu)
{
	return -log(mu) - x/mu;
}

/// Draws a normally distributed number with mean mu and standard deviation sd
double normal_sample(const double mu, const double sd)
{
	normal_distribution<double> distribution(mu,sd);
	return distribution(generator);
}

/// The log of the probability from the normal distribution
double normal_probability(const double x, const double mean, const double var)
{
  if(var <= 0) emsg("Negative");
  return -0.5*log(2*M_PI*var) - (x-mean)*(x-mean)/(2*var);
}

/// Draws a random number between 0 and 1
double ran()
{
	return double(0.999999999*rand())/RAND_MAX;
}


vector <double> logsum;

/// Initialises calculation of choose function
void choose_init()                 
{
	auto Nmax = 0u;
	for(const auto &he : herd){
		if(he.ind.size() > Nmax) Nmax = he.ind.size();
	}
	Nmax++;
	logsum.resize(Nmax);
	
	auto sum = 0.0;
	for(auto n = 1u; n < Nmax; n++){
		sum += log(n);
		logsum[n] = sum;
	}
}

/// Calculates N choose n
double choose(int N, int n)
{
	if(N < 1 || n < 1 || n > N) emsg("choose problem"); 
	return logsum[N] - logsum[N-n] - logsum[n];
}
	
/// Displays an error message
void emsg(const string& msg)
{
	cout << msg << endl;
	exit (EXIT_FAILURE);
}

/// Performs proposals which changes T0 and alters other event time 
void Herd::joint_T0_event_proposal()
{	
	auto T0_store = T0;
	
	auto d_T0 = normal_sample(0,T0_joint_jump);
	
	T0 += d_T0;
		
	auto T_end = T_cull;                                  // This is the maximum time over which scaling is performed

	auto ind_store = ind;

	auto illegal = false;

	auto dLi_gamma = 0.0;

	auto fac = 0.0;  // This is used to account for changes in phase space under the proposal
	for(auto &in : ind){
		if(in.inf == true){
			auto t = in.Rec_t;
			if(t > T_end) t = T_end;
			auto T = t-T0_store;
			
			if(in.Exp_t == T0_store){
				in.Exp_t = T0;
				if(T0 > t){ illegal = true; break;}
			}
			else{
				auto frac = (t-in.Exp_t)/T;
				if(frac < 0) illegal = true;
				fac += log(1-(d_T0/T));
				in.Exp_t += d_T0*frac; 
				if(in.Exp_t > T_cull){ illegal = true; break;}
			}
					
			auto frac = (t-in.Inf_t)/T;
			if(frac > 0){
				fac += log(1-(d_T0/T));
				in.Inf_t += d_T0*frac; 
			}
			
			in.latentP = in.Inf_t - in.Exp_t;
			in.infectiousP = in.Rec_t - in.Inf_t;
		
			auto Li_gamma_prop = in.likelihood_gamma();
			dLi_gamma += Li_gamma_prop - in.Li_gamma;
			in.Li_gamma = Li_gamma_prop;
			
			if(in.latentP <= 0 || in.infectiousP <= 0){ illegal = true; break;}
		}
	}
		
	double al, Li_prop;
	if(illegal == true) al = 0;
	else{	
		Li_prop = likelihood();   
		
		al = exp(phi_L*(Li_prop - Li) + phi_G*dLi_gamma + fac); 
	}

	nT0_joint_propose++;
	if(ran() < al){
		nT0_joint_accept++;
		Li = Li_prop;
		if(s < nburnin) T0_joint_jump *= 1.001;                // The dynamically adapts the jumping size
	}
	else{
		T0 = T0_store;
		ind = ind_store;
		if(s < nburnin) T0_joint_jump *= 0.9995;
	}
}


/// This checks that all quantities are correctly specified
void check(int num)
{
	double dd;
	
	for(auto &he : herd){
		for(const auto &in : he.ind){                                         // Checks timings are correctly specified
			if(in.latentP < 0) emsg("Problem latentP");
			if(in.infectiousP < 0) emsg("Problem infectiousP");
			
			dd = in.latentP -	(in.Inf_t - in.Exp_t);
			if(dd*dd > TINY) emsg("Problem latentP2");
		
			dd = in.infectiousP -	(in.Rec_t - in.Inf_t);
			if(dd*dd > TINY)	emsg("Problem infectiousP2");
		}
		
		auto T = LARGE;                                                   // CHecks T0 correctly specified
		for(const auto &in : he.ind){
			if(in.inf == true){
				if(in.Exp_t < T) T = in.Exp_t;
				
				if(in.Exp_t > he.T_cull) emsg("Infection must be before T_cull");
				
				if(in.obs == true){
					if(in.Rec_t > he.T_cull) emsg("Observed must recover before T_cull");
				}
				else{
					if(in.Rec_t < he.T_cull) emsg("Unobserved must recover after T_cull");
				}
			}
		}
		if(T != he.T0) emsg ("T0 wrong");
		
		dd = he.likelihood() - he.Li;	
		if(dd*dd > TINY){
			emsg("Problem likelihood"+to_string(num));
		}
		
		for(auto &in : he.ind){     
			dd = in.likelihood_gamma() - in.Li_gamma;
			if(dd*dd > TINY){
				emsg("Problem likelihood_gamma");
			}
		}	
		
		dd = he.prior_index() - he.Pr_index;	
		if(dd*dd > TINY) emsg("Problem Pr_index");
	}
	
	dd = prior() - Pr;	
	if(dd*dd > TINY) emsg("Problem Pr");
}