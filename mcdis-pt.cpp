// #include "mpi.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <complex>
#include <math.h>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <fstream>
#include <functional>

#include <omp.h>

static inline int is_odd(int x) { return x & 1; }

#include "Boundaries.h"
#include "params.h"

using namespace std;

int main(int argc, char **argv) {
	cout << endl;
	clock_t t_total, t; // to measure runtime of simulation
	t_total = clock(); // start the clock t_total

	/* Read in parameters from terminal */
	if (argc == 17) {
    	dim = atof(argv[1]); // dimensionality of the system
    	L = atof(argv[2]); // number of sites in x and y direction
    	Lz = atof(argv[3]); // number of sites in z direction
    
    	J33 = atof(argv[4]); // AF between spin-1 satellites
    	J34 = atof(argv[5]); // FM between core and satellites
    	J44 = atof(argv[6]); // FM between two different core sites
    
    	T_min = atof(argv[7]); // temperature of MC simulation
    	T_max = atof(argv[8]); // temperature of MC simulation
    	T_steps = atof(argv[9]); // temperature of MC simulation

    	disorder_filling = atof(argv[10]); // site filling considered
    	disorder_steps = atof(argv[11]); // number of different disorder realizations considered
    
    	mc_bins = atof(argv[12]); // number of MC steps (per disorder realization)
    	mc_binsize = atof(argv[13]); // number of MC steps in each bin 
    	thermalization_bins = atof(argv[14]); 
    
    	seed_disorder = atof(argv[15]); // initial seed of disorder RNG (determining which sites are filled)
    	seed_mc = atof(argv[16]); // initial seed of MC RNG

    	output_s_q = false; 
		output_observables = true;
		output_data_fraction = true;
    	print_all = false;
		snapshot_output = false; 
		parallel_tempering_step = true; // whether or not to include a parallel tempering step
	}
	else {
		cout << "Default parameters" << endl;
		dim = 3;
		L = 8;
		Lz = 8;
	
		J33 = 0.2;
		J34 = -1.;
		J44 = -1.;

		T_min = 0.5;
		T_max = 6.;
		T_steps = 2;

		disorder_filling = 0.2;
		disorder_steps = 50;

		mc_bins = 1;
		mc_binsize = 1;
		thermalization_bins = 1; 
	
		seed_disorder = 1; // increases by one for each disorder realization (final seed = disorder_steps)
		seed_mc = 1; // 

		output_s_q = false; 
		output_observables = true;
		output_data_fraction = true;
		print_all = true;
		snapshot_output = true;
		parallel_tempering_step = true; 
	}

	N2D = L*L; // total number of sites in 2D system
  	N3D = L*L*Lz; // total number of sites in 3D system
	string outputFilename;
  	// Define number of neighbors (here for 3D cubic lattice)
  	if ( (dim == 2) || (Lz == 1) ) {
  		num_neighbors_direct = 4;
  		num_neighbors_face = 4; 
		num_neighbors_body = 4; // use this as being equal to num_neighbors_face 
	  	num_neighbors_direct_nn = 4; 

	  	N = N2D;
	  	Lz = 1; // set Lz = 1 in case of dim=2 to get correct output values
	}
	else {
  		num_neighbors_direct = 6; 
  		num_neighbors_face = 12; 
  		num_neighbors_body = 8; 
  		num_neighbors_direct_nn = 6; 

	  	N = N3D;
	}

	// initialize arrays for neighboring site labels
  	neighbors_direct = new int[num_neighbors_direct*N]; // initialize array for direct neighbors
  	neighbors_face = new int[num_neighbors_face*N];
  	neighbors_body = new int[num_neighbors_body*N];
  	neighbors_direct_nn = new int[num_neighbors_direct_nn*N];

	// Initialize boundary class to initialize all neighbor table
	Boundaries boundaries(dim, L, Lz, neighbors_direct, neighbors_face, neighbors_body, neighbors_direct_nn);

	// Total number of MCS
  	mc_steps = (mc_bins + thermalization_bins)*mc_binsize; // total number of MCS

  	// Setting up table of temperatures (for parallel tempering)
  	temperatures = new double[T_steps];

  	for (int idx = 0; idx < T_steps; idx++) {
		temperatures[idx] = 0.;
	}
	temperatures[0] = T_min;
	temperatures[T_steps - 1] = T_max;
	for (int idx = 1; idx < T_steps - 1; idx++) {
		temperatures[idx] = T_min*pow(T_max/T_min, static_cast<double>(idx)/static_cast<double>(T_steps - 1));
	}

/*
	for (int idx = 0; idx < T_steps; idx++) {
		cout << "temperatures[" << idx << "] = " << temperatures[idx] << endl;
	}
*/
	int *pt_TtoE; // list where entry i gives ensemble at temperature temperatures[i]
	int *pt_EtoT; // list where entry i gives temperature of ensemble i
	swap_even_pairs = false;
	//cout << "swap_even_pairs = " << static_cast<int>(swap_even_pairs) << endl;

	pt_TtoE = new int[2*T_steps]; // factor of 2 because of 2 replicas per T
	pt_EtoT = new int[2*T_steps];

	for (int idx = 0; idx < 2*T_steps; idx++) {
		pt_TtoE[idx] = idx % T_steps;
		pt_EtoT[idx] = idx % T_steps;
		//cout << "pt_TtoE[" << idx << "] = " << pt_TtoE[idx] << endl;
		//cout << "pt_EtoT[" << idx << "] = " << pt_EtoT[idx] << endl;
	}

	/**********************************************************************************
	*** Set up table of spin flip probabilities used for Metropolis algorithm *********
	**********************************************************************************/
	// Construct flip tables: one for each temperature: prob_1 is for core sites S^{4}, prob_2 is for satellite sites S^{3}
    prob_1 = new double[T_steps*2*(2*num_neighbors_direct + 1) * (2*num_neighbors_direct + 1)]; // S^(4)
    prob_2 = new double[T_steps*2*(2*num_neighbors_direct + 1) * (2*num_neighbors_direct + 1)]; // S^(3)
	
	// Initialize 
	for (int idx = 0; idx < T_steps*2*(2*num_neighbors_direct + 1) * (2*num_neighbors_direct + 1); idx++) {
		prob_1[idx] = 0.; // flip prob. for spin S^(4) with S=1/2: core spin
		prob_2[idx] = 0.; // flip prob. for spin S^(3) with S=1: satellite spin
	}

	for (int T_idx = 0; T_idx < T_steps; T_idx++) {
		for (int bdx = -num_neighbors_direct; bdx <= num_neighbors_direct; bdx++) { // bdx = total spin of neighboring S^(4) spins (counting up=+1 and down=-1). Factor of 1/2 is accounted for below in equation
			for (int adx = -num_neighbors_direct; adx <= num_neighbors_direct; adx++) { // adx = total spin of neighboring S^(3) spins (counting up=+1 and down=-1).
				for (int idx = 0; idx <= 1; idx++) { // idx = direction of central spin: idx=0=spin down, idx=1=spin up
				// according to Eq.(36) of Sandvik's short MC notes. Note |S^(4)|=1, |S^(3)|=1 here. For old convention of spin sizes see oldcode.cpp
					prob_1[idx + (adx + num_neighbors_direct)*2 + (bdx + num_neighbors_direct)*2*(2*num_neighbors_direct + 1) + T_idx*2*(2*num_neighbors_direct + 1) * (2*num_neighbors_direct + 1)] = exp( (2.*static_cast<double>(2*idx - 1)*(J34*static_cast<double>(adx) + J44*static_cast<double>(bdx)) )/temperatures[T_idx] ); 
					prob_2[idx + (adx + num_neighbors_direct)*2 + (bdx + num_neighbors_direct)*2*(2*num_neighbors_direct + 1) + T_idx*2*(2*num_neighbors_direct + 1) * (2*num_neighbors_direct + 1)] = exp( (2.*static_cast<double>(2*idx - 1)*(J33*static_cast<double>(adx) + J34*static_cast<double>(bdx)) )/temperatures[T_idx] );
				}
			}
		}
	}
	/*
	// Output results 
	outputFilename = "Data-ProbabilityTable.dat";
	outputFile.open(outputFilename);
	outputFile.setf(ios::scientific);
	// Output probability tables
	for (int idx = 0; idx < T_steps*2*(2*num_neighbors_direct + 1) * (2*num_neighbors_direct + 1); idx++) {
		outputFile << prob_1[idx] << " " << prob_2[idx] << endl; 
	}
	*/

	// Construct spin array (integer)
	spins = new int[T_steps*2*N]; // 2 replicas per temperature
	order = new int[T_steps*N]; // need only half as many order as the two replicas have the same disorder realization
	spins_energy = new double[T_steps*2]; // array that stores energy of spin ensemble


	vector<mt19937> mt_rand_mc_vector;

	/**************************************************
     *** Create arrays used during MC averaging *******
     *************************************************/
  	magnetization1 = new double[2*T_steps*mc_bins]; // stores magnetization1 for each MC bin at each temperature. Is for spin 1 at each site. Factor of two for the replicas.
  	magnetization2 = new double[2*T_steps*mc_bins]; // stores magnetization2 for each MC bin. Is for spin 3/2 at core Co^{4+} sites. Factor of two for the replicas.
  	magnetizationAF = new double[2*T_steps*mc_bins]; // stores AF magnetization at Q = (pi,pi,pi). Factor of two for the replicas.
  	magnetization_replica = new double[T_steps*mc_bins]; // stores m=1/N\sum_{i} m^{1}_i m^{2}_i correlations between two replicas
  	magnetization1_squared = new double[2*T_steps*mc_bins]; // stores magnetization1 for each MC bin. Is for spin 1 at each site.
  	magnetization2_squared = new double[2*T_steps*mc_bins]; // stores magnetization2 for each MC bin. Is for spin 3/2 at core Co^{4+} sites.
  	magnetizationAF_squared = new double[2*T_steps*mc_bins]; // stores AF magnetization^2
  	magnetization_replica_squared = new double[T_steps*mc_bins];
  	magnetization1_four = new double[2*T_steps*mc_bins]; // stores magnetization1^4 for each MC bin. Needed for Binder ratio of FM order.
  	magnetizationAF_four = new double[2*T_steps*mc_bins]; // stores magnetizationAF^4 for each MC bin. Needed for Binder ratio of AF order.
  	magnetization_replica_four = new double[T_steps*mc_bins];

  	energy = new double[2*T_steps*mc_bins]; // stores average energy for each MC bin. Factor of two for the replicas.
  	energy_squared = new double[2*T_steps*mc_bins]; // stores averag energy^2 for each MC bin
  
  	s_q = new double[T_steps*mc_bins*L*Lz*(L/2 + 1)];
  	xi_a = new double[T_steps*mc_bins];
  	xi_b = new double[T_steps*mc_bins];

	/*****************************************************
	******	MC averages for each disorder realization ****
	*****************************************************/
	magnetization1_av_dis = new double[2*T_steps*disorder_steps]; // stores average magnetization1 for each disorder realization. Is for spin 1 at each site
	magnetization2_av_dis = new double[2*T_steps*disorder_steps]; // is for spin 3/2 at core Co^{4+} sites.
	magnetizationAF_av_dis = new double[2*T_steps*disorder_steps]; 
	magnetization_replica_av_dis = new double[T_steps*disorder_steps];

	sigma_magnetization1_av_dis = new double[2*T_steps*disorder_steps]; // stores standard deviation sigma of magnetization1 for each disorder realization. Is for spin 1 at each site
	sigma_magnetization2_av_dis = new double[2*T_steps*disorder_steps];
	sigma_magnetizationAF_av_dis = new double[2*T_steps*disorder_steps];
	sigma_magnetization_replica_av_dis = new double[T_steps*disorder_steps];

	magnetization1_squared_av_dis = new double[2*T_steps*disorder_steps]; // stores average magnetization1 for each disorder realization. Is for spin 1 at each site
	magnetization2_squared_av_dis = new double[2*T_steps*disorder_steps]; // is for spin 3/2 at core Co^{4+} sites.
	magnetizationAF_squared_av_dis = new double[2*T_steps*disorder_steps]; 
	magnetization_replica_squared_av_dis = new double[T_steps*disorder_steps];
	magnetization1_four_av_dis = new double[2*T_steps*disorder_steps]; 
	magnetizationAF_four_av_dis = new double[2*T_steps*disorder_steps]; 
	magnetization_replica_four_av_dis = new double[T_steps*disorder_steps];

	sigma_magnetization1_squared_av_dis = new double[2*T_steps*disorder_steps]; // stores standard deviation sigma of magnetization1_squared for each disorder realization. Is for spin 1 at each site
	sigma_magnetization2_squared_av_dis = new double[2*T_steps*disorder_steps]; 
	sigma_magnetizationAF_squared_av_dis = new double[2*T_steps*disorder_steps]; 
	sigma_magnetization_replica_squared_av_dis = new double[T_steps*disorder_steps];
	sigma_magnetization1_four_av_dis = new double[2*T_steps*disorder_steps]; 
	sigma_magnetizationAF_four_av_dis = new double[2*T_steps*disorder_steps]; 
	sigma_magnetization_replica_four_av_dis = new double[T_steps*disorder_steps];

	chi_magnetization1_av_dis = new double[2*T_steps*disorder_steps]; 
	sigma_chi_magnetization1_av_dis = new double[2*T_steps*disorder_steps]; 
	chi_magnetizationAF_av_dis = new double[2*T_steps*disorder_steps]; 
	sigma_chi_magnetizationAF_av_dis = new double[2*T_steps*disorder_steps]; 
	chi_magnetization_replica_av_dis = new double[T_steps*disorder_steps]; 
	sigma_chi_magnetization_replica_av_dis = new double[T_steps*disorder_steps]; 

	binder_ratio_av_dis = new double[2*T_steps*disorder_steps]; 
	sigma_binder_ratio_av_dis = new double[2*T_steps*disorder_steps]; 
	binder_ratioAF_av_dis = new double[2*T_steps*disorder_steps]; 
	sigma_binder_ratioAF_av_dis = new double[2*T_steps*disorder_steps]; 
	binder_ratio_replica_av_dis = new double[T_steps*disorder_steps]; 
	sigma_binder_ratio_replica_av_dis = new double[T_steps*disorder_steps];

	energy_av_dis = new double[2*T_steps*disorder_steps]; 
	energy_squared_av_dis = new double[2*T_steps*disorder_steps]; 

	specific_heat_av_dis = new double[2*T_steps*disorder_steps]; 
	sigma_specific_heat_av_dis = new double[2*T_steps*disorder_steps]; 

	//spins_q_av_dis = new double[disorder_steps*N];
	if (output_s_q == true) {
		s_q_av_dis = new double[T_steps*disorder_steps*(L*Lz*(L/2 + 1))];
		cout << "T_steps*disorder_steps*(L*Lz*(L/2 + 1)) = " << T_steps*disorder_steps*(L*Lz*(L/2 + 1));
	}
	xi_a_av_dis = new double[T_steps*disorder_steps];
	xi_b_av_dis = new double[T_steps*disorder_steps];
	sigma_xi_a_av_dis = new double[T_steps*disorder_steps];
	sigma_xi_b_av_dis = new double[T_steps*disorder_steps];

	fraction_filled = new double[disorder_steps];
	fraction_bonds_J33 = new double[disorder_steps];
	fraction_bonds_J34 = new double[disorder_steps];
	fraction_bonds_J44 = new double[disorder_steps];

	/***************************************************
	*** Initialize av_dis variables ********************
	***************************************************/
	for (int idx = 0; idx < 2*T_steps*disorder_steps; idx++) {
		magnetization1_av_dis[idx] = 0.;
		magnetization2_av_dis[idx] = 0.;
		magnetizationAF_av_dis[idx] = 0.;
		sigma_magnetization1_av_dis[idx] = 0.;
		sigma_magnetization2_av_dis[idx] = 0.;
		sigma_magnetizationAF_av_dis[idx] = 0.;
		magnetization1_squared_av_dis[idx] = 0.;
		magnetization2_squared_av_dis[idx] = 0.;
		magnetization1_four_av_dis[idx] = 0.;
		magnetizationAF_four_av_dis[idx] = 0.;
		sigma_magnetization1_squared_av_dis[idx] = 0.;
		sigma_magnetization2_squared_av_dis[idx] = 0.;
		sigma_magnetizationAF_squared_av_dis[idx] = 0.;
		sigma_magnetization1_four_av_dis[idx] = 0.;
		sigma_magnetizationAF_four_av_dis[idx] = 0.;
		chi_magnetization1_av_dis[idx] = 0.;
		sigma_chi_magnetization1_av_dis[idx] = 0.;
		chi_magnetizationAF_av_dis[idx] = 0.;
		sigma_chi_magnetizationAF_av_dis[idx] = 0.;
		binder_ratio_av_dis[idx] = 0.;
		sigma_binder_ratio_av_dis[idx] = 0.;
		binder_ratioAF_av_dis[idx] = 0.;
		sigma_binder_ratioAF_av_dis[idx] = 0.;
		energy_av_dis[idx] = 0.;
		energy_squared_av_dis[idx] = 0.;
		specific_heat_av_dis[idx] = 0.;
		sigma_specific_heat_av_dis[idx] = 0.;
	}

	for (int idx = 0; idx < T_steps*disorder_steps; idx++) {
		magnetization_replica_av_dis[idx] = 0.;
		sigma_magnetization_replica_av_dis[idx] = 0.;
		magnetization_replica_squared_av_dis[idx] = 0.;
		magnetization_replica_four_av_dis[idx] = 0.;
		sigma_magnetization_replica_squared_av_dis[idx] = 0.;
		sigma_magnetization_replica_four_av_dis[idx] = 0.;
		chi_magnetization_replica_av_dis[idx] = 0.;
		sigma_chi_magnetization_replica_av_dis[idx] = 0.;
		binder_ratio_replica_av_dis[idx] = 0.;
		sigma_binder_ratio_replica_av_dis[idx] = 0.;
		if (output_s_q == true) {
			for (int jdx = 0; jdx < L*Lz*(L/2+1); jdx++) {
				s_q_av_dis[jdx + idx*(L*Lz*(L/2+1))] = 0.;
			}
		}	
		xi_a_av_dis[idx] = 0.;
		xi_b_av_dis[idx] = 0.;
		sigma_xi_a_av_dis[idx] = 0.;
		sigma_xi_b_av_dis[idx] = 0.;
}

for (int idx = 0; idx < disorder_steps; idx++) {
	fraction_filled[idx] = 0.;
	fraction_bonds_J33[idx] = 0.;
	fraction_bonds_J34[idx] = 0.;
	fraction_bonds_J44[idx] = 0.;
}

	/**************************************************************************
	***** Final results averaged over MC and disorder ****************
	**************************************************************************/
	magnetization1_av = new double[T_steps];
	magnetization2_av = new double[T_steps];
	magnetizationAF_av = new double[T_steps];
	magnetization_replica_av = new double[T_steps];

	sigma_magnetization1_av = new double[T_steps];
	sigma_magnetization2_av = new double[T_steps];
	sigma_magnetizationAF_av = new double[T_steps];
	sigma_magnetization_replica_av = new double[T_steps];

	magnetization1_squared_av = new double[T_steps];
	magnetization2_squared_av = new double[T_steps];
	magnetizationAF_squared_av = new double[T_steps];
	magnetization_replica_squared_av = new double[T_steps];
	magnetization1_four_av = new double[T_steps];
	magnetizationAF_four_av = new double[T_steps];
	magnetization_replica_four_av = new double[T_steps];

	sigma_magnetization1_squared_av = new double[T_steps];
	sigma_magnetization2_squared_av = new double[T_steps];
	sigma_magnetizationAF_squared_av = new double[T_steps];
	sigma_magnetization_replica_squared_av = new double[T_steps];
	sigma_magnetization1_four_av = new double[T_steps];
	sigma_magnetizationAF_four_av = new double[T_steps];
	sigma_magnetization_replica_four_av = new double[T_steps];

	chi_magnetization1_av = new double[T_steps]; // magnetic susceptibility
	sigma_chi_magnetization1_av = new double[T_steps];
	chi_magnetizationAF_av = new double[T_steps]; // magnetic susceptibility
	sigma_chi_magnetizationAF_av = new double[T_steps];
	chi_magnetization_replica_av = new double[T_steps]; // magnetic susceptibility
	sigma_chi_magnetization_replica_av = new double[T_steps];

	binder_ratio_av = new double[T_steps]; // Binder ratio for FM order
	sigma_binder_ratio_av = new double[T_steps];
	binder_ratioAF_av = new double[T_steps];// Binder ratio for AF order
	sigma_binder_ratioAF_av = new double[T_steps];
	binder_ratio_replica_av = new double[T_steps]; // Binder ratio of replica order
	sigma_binder_ratio_replica_av = new double[T_steps];

	energy_av = new double[T_steps];
	energy_squared_av = new double[T_steps];

	specific_heat_av = new double[T_steps]; // specific heat
	sigma_specific_heat_av = new double[T_steps];

	s_q_av = new double[T_steps*L*Lz*(L/2+1)]; 
	for (int idx = 0; idx < T_steps*L*Lz*(L/2+1); idx++) {
		s_q_av[idx] = 0.;
	}
	xi_a_av = new double[T_steps];
	xi_b_av = new double[T_steps];

	sigma_xi_a_av = new double[T_steps];
	sigma_xi_b_av = new double[T_steps];

	for (int idx = 0; idx < T_steps; idx++) {
		magnetization1_av[idx] = 0.;
		magnetization2_av[idx] = 0.;
		magnetizationAF_av[idx] = 0.;
		magnetization_replica_av[idx] = 0.;

		sigma_magnetization1_av[idx] = 0.;
		sigma_magnetization2_av[idx] = 0.;
		sigma_magnetizationAF_av[idx] = 0.;
		sigma_magnetization_replica_av[idx] = 0.;

		magnetization1_squared_av[idx] = 0.;
		magnetization2_squared_av[idx] = 0.;
		magnetizationAF_squared_av[idx] = 0.;
		magnetization_replica_squared_av[idx] = 0.;
		magnetization1_four_av[idx] = 0.;
		magnetizationAF_four_av[idx] = 0.;
		magnetization_replica_four_av[idx] = 0.;

		sigma_magnetization1_squared_av[idx] = 0.;
		sigma_magnetization2_squared_av[idx] = 0.;
		sigma_magnetizationAF_squared_av[idx] = 0.;
		sigma_magnetization_replica_squared_av[idx] = 0.;
		sigma_magnetization1_four_av[idx] = 0.;
		sigma_magnetizationAF_four_av[idx] = 0.;
		sigma_magnetization_replica_four_av[idx] = 0.;

		chi_magnetization1_av[idx] = 0.; // magnetic susceptibility
		sigma_chi_magnetization1_av[idx] = 0.;
		chi_magnetizationAF_av[idx] = 0.; // magnetic susceptibility
		sigma_chi_magnetizationAF_av[idx] = 0.;
		chi_magnetization_replica_av[idx] = 0.; // magnetic susceptibility
		sigma_chi_magnetization_replica_av[idx] = 0.;

		binder_ratio_av[idx] = 0.; // Binder ratio for FM order
		sigma_binder_ratio_av[idx] = 0.;
		binder_ratioAF_av[idx] = 0.;// Binder ratio for AF order
		sigma_binder_ratioAF_av[idx] = 0.;
		binder_ratio_replica_av[idx] = 0.; // Binder ratio of replica order
		sigma_binder_ratio_replica_av[idx] = 0.;

		energy_av[idx] = 0.;
		energy_squared_av[idx] = 0.;

		specific_heat_av[idx] = 0.; // specific heat
		sigma_specific_heat_av[idx] = 0.;

		xi_a_av[idx] = 0.;
		xi_b_av[idx] = 0.;

		sigma_xi_a_av[idx] = 0.;
		sigma_xi_b_av[idx] = 0.;
	}

	/***** End initialize FINAL RESULT VARIABLES *****/

	// Initialize RNGs
  	// Disorder RNG: Set up RNG "mt_rand_disorder" by first using common lc generator that procudes a seed for the mersenne-twister rng
  	std::minstd_rand0 lc_generator(seed_disorder);
	std::uint_least32_t seed_data[std::mt19937::state_size];
	std::generate_n(seed_data, std::mt19937::state_size, std::ref(lc_generator));
	std::seed_seq q(std::begin(seed_data), std::end(seed_data));
	std::mt19937 mt_rand_disorder(q);

	// Set up a RNG for each MC simulation (2*T_steps)
	for (int idx = 0; idx < 2*T_steps; idx++) {
		// Monte Carlo: Set up RNG "mt_rand_mc" by first using common lc generator that procudes a seed for the mersenne-twister rng
		std::minstd_rand0 lc_generator_mc(seed_mc + idx ); // seed different for each temperature and replica
		std::uint_least32_t seed_data_mc[std::mt19937::state_size];
		std::generate_n(seed_data_mc, std::mt19937::state_size, std::ref(lc_generator_mc));
		std::seed_seq q_mc(std::begin(seed_data_mc), std::end(seed_data_mc));
		mt_rand_mc_vector.push_back(std::mt19937(q_mc));
	}

	std::uniform_real_distribution<double> dis(0.0,1.0); // random distribution of doubles in interval [0, 1]

	cout << endl;
	cout << "dim = " << dim << ", L = " << L << ", Lz = " << Lz << ", J33 = " << J33 << ", J34 = " << J34 << ", J44 = " << J44 << ", T_min = " << T_min << ", T_max = " << T_max << ", T_steps = " << T_steps << ", disorder_filling = " << disorder_filling << ", thermalization_bins = " << thermalization_bins << ", mc_bins = " << mc_bins << ", mc_binsize = " << mc_binsize << ", mc_steps = " << mc_steps << ", disorder_steps = " << disorder_steps << ", disorder_seed = " << seed_disorder << ", mc_seed = " << seed_mc << endl;

	string Lz_string = to_string(Lz);
	string L_string = to_string(L);
	string mc_steps_string = to_string(mc_steps);
	string J33_string = to_string(J33);
	string J34_string = to_string(J34);
	string J44_string = to_string(J44);
	string T_min_string = to_string(T_min);
	string T_max_string = to_string(T_max);
	string T_steps_string = to_string(T_steps);
	string disorder_filling_string = to_string(disorder_filling);

/*************************************************************************
*** MC SIMULATION FOR FIXED DISORDER REALIZATION (Metropolis + PT) *******
*************************************************************************/
	for (int idx_disorder = 0; idx_disorder < disorder_steps; idx_disorder++) {
		// Fill lattice randomly
		// perform MC simulation of this configuration
		// average and compute errors
		// output av_dis variables (MC results for this fixed disorder realization)

		int s1, s2; // variables needed to set up random order of spin filling (permute array order[i])

		// Initialize variables (magnetization per bin, etc)
		for (unsigned long int idx = 0; idx < 2*T_steps*mc_bins; idx++) {
			magnetization1[idx] = 0.;
			magnetization2[idx] = 0.;
			magnetizationAF[idx] = 0.;
			magnetization1_squared[idx] = 0.;
			magnetization2_squared[idx] = 0.; 
			magnetizationAF_squared[idx] = 0.; 
			magnetization1_four[idx] = 0.;
			magnetizationAF_four[idx] = 0.;
			energy[idx] = 0.;
			energy_squared[idx] = 0.;
		}

		// Initialize variables (magnetization per bin, etc)
		for (unsigned long int idx = 0; idx < T_steps*mc_bins; idx++) {
			magnetization_replica[idx] = 0.;	
			magnetization_replica_squared[idx] = 0.; 
			magnetization_replica_four[idx] = 0.;
			for (int jdx = 0; jdx < L*Lz*(L/2+1); jdx++) {
				s_q[jdx + idx*L*Lz*(L/2+1)] = 0.;
			}
			xi_a[idx] = 0.;
			xi_b[idx] = 0.;
		} // for initialize all arrays containing mc bin results
		acceptance_number_Metropolis = 0;
  		acceptance_number_pt = 0;

		/********************************************************************************************************************
		***** FILL LATTICE RANDOMLY. Different T have same disorder realization.                                          ***
		***** Determine the order in which the sites will be filled with Co^{4+}. Set up permuted list of sites: order[i] ***
		********************************************************************************************************************/
		int N_max = round(disorder_filling*N); // number of core spins that are put in the system
		//cout << "N_max = " << N_max << endl;

		for (int idx = 0; idx < T_steps*2*N; idx++) {
			spins[idx] = 0; // initialize all spins being empty
		}
	
		for (int idx = 0; idx < N; idx++) {
			order[idx] = idx;
		}
		
		// randomly permute the array order[i]
   		for (int i = 0; i < N; i++) {
   			int j = i + (N-i)*dis(mt_rand_disorder); // use Merseene Twister RNG to generate uniform distribution in interval [0, 1]
	   		int temp = order[i];
   			order[i] = order[j];
   			order[j] = temp;
   		}

	for (int i = 0; i < N_max; i++) { // loop over all sites of the lattice 
    	s1 = order[i]; // s1 are labels of the next site to be filled in the lattice
    	// site labels are: {0=empty, 1=+1/2=Co^(4+) up, -1=-1/2=Co^(4+) down, 2=+1=Co^(3+) up, -2=-1=Co^(3+) down
    	
		// Do not parallelize this loop (using #pragma as it gives wrong results! 
    	// Fill site s1 with spin up or down in ensembles at all T and replicas (2*T_steps)
    	for (int T_idx = 0; T_idx < T_steps; T_idx++) {
    		for (int replica_idx = 0; replica_idx <= 1; replica_idx++) {
    		// create core Co^(4+_ spin up and spin down with equal probabilities
    			if (dis(mt_rand_disorder) < 0.5 ){
	    			spins[T_idx*2*N + replica_idx*N + s1] = 1; // create spin up at site s1
	    		}
	    		else {
	    			spins[T_idx*2*N + replica_idx*N + s1] = -1; // create spin down at site s1 
				}
			// Update all empty direct nn sites of \pm 2 (since they are part of a polaron). Randomly initialize satellite spin state.
				for (int jdx = 0; jdx < num_neighbors_direct; jdx++) {
    				s2 = neighbors_direct[num_neighbors_direct*s1 + jdx]; // direct neighbor of site s1 (one of num_neighbors_direct = 6)
    				if (spins[T_idx*2*N + replica_idx*N + s2] == 0) {
						if (dis(mt_rand_disorder) < 0.5 ) {
    						spins[T_idx*2*N + replica_idx*N + s2] = 2; // put direct neighbors into \pm 2 state (Co^{3+} with same direction as polaron core site \pm 1)	
						}
						else {
							spins[T_idx*2*N + replica_idx*N + s2] = -2;
						}
					} // if spins[s2] == 0 
   				} // for jdx over direct neighbors
   			} // for replica_idx 
   		} // for T_idx
   	} // for i (runs over filled sites)

	   // output initial spin state (for debugging)
	/*for (int replica_idx = 0; replica_idx <= 1; replica_idx++) {
		cout << "Replica " << replica_idx << ": " << endl;
	   	for (int idx = 0; idx < N; idx++) {
			cout << "spins[" << idx << "] = " << spins[replica_idx*N + idx] << endl;
	   	}
		cout << endl;
	}
	*/
	
	// Count number of filled sites (= number of spins) for this disorder realization. All 2*T_steps ensembles have the same disorder realization. 
	int counter_filled = 0;
	vector<int> filled_sites(N,0);
	for (int idx = 0; idx < N; idx++) {
		// generate list of filled sites (these contain a spin and can be flipped)
		if (spins[idx] != 0) {
			filled_sites[counter_filled] = idx; 
			counter_filled++;
		}
	}
	if (print_all == true) {
		cout << "Total number of filled sites = " << counter_filled << ". Fraction of filled sites = " << static_cast<double>(counter_filled)/static_cast<double>(N) << endl;
	}

	fraction_filled[idx_disorder] = static_cast<double>(counter_filled)/static_cast<double>(N);

	// Measure histogram distribution of bond interactions {J33, J34, J44}
	int number_bonds_J33 = 0;
	int number_bonds_J34 = 0;
	int number_bonds_J44 = 0;

	for (int idx = 0; idx < counter_filled; idx++) {
		int site = filled_sites[idx];
		for (int neighbor_idx = 0; neighbor_idx < num_neighbors_direct; neighbor_idx++) {
			int neighbor_site = neighbors_direct[num_neighbors_direct*site + neighbor_idx];
			if (abs(spins[site]) == 1) { // central spin is Co^{4+}
				if (abs(spins[neighbor_site]) == 1) {
					number_bonds_J44++;
				}
				else if (abs(spins[neighbor_site]) == 2) {
					number_bonds_J34++;
				}
				else {
					// do nothing as neighbor site is empty
				}
				}
			else if (abs(spins[site]) == 2) { // central spin is Co^{3+}
				if (abs(spins[neighbor_site]) == 1) {
					number_bonds_J34++;
				}
				else if (abs(spins[neighbor_site]) == 2) {
					number_bonds_J33++;
				}
				else {
					// do nothing as neighbor site is empty
				}	
			}
			else {
				// do nothing as core site is empty (should not occur anyway as I am summing over filled sites only)
			}
		} // for neighbor_idx
	} // for idx (over filled sites)

	// compute fraction from number of bonds
	fraction_bonds_J33[idx_disorder] = static_cast<double>(number_bonds_J33)/static_cast<double>((num_neighbors_direct*N));
	fraction_bonds_J34[idx_disorder] = static_cast<double>(number_bonds_J34)/static_cast<double>((num_neighbors_direct*N));
	fraction_bonds_J44[idx_disorder] = static_cast<double>(number_bonds_J44)/static_cast<double>((num_neighbors_direct*N));

	if (print_all == true) {
		cout << "Number of bonds J33 = " << number_bonds_J33 << ", Number of bonds J34 = " << number_bonds_J34 << ", Number of bonds J44 = " << number_bonds_J44 << endl;
		cout << "Fraction J33 = " << fraction_bonds_J33[idx_disorder] << ", Fraction J34 = " << fraction_bonds_J34[idx_disorder] << ", Fraction J44 = " << fraction_bonds_J44[idx_disorder] << endl;
		// cout << "Fraction Av J33 = " << fraction_bonds_J33_av << ", Fraction Av J34 = " << fraction_bonds_J34_av << ", Fraction Av J44 = " << fraction_bonds_J44_av << ", Fraction Bonds Av All = " << fraction_bonds_J33_av + fraction_bonds_J34_av + fraction_bonds_J44_av << endl;
	}
	// Output snapshot of spin state
	// Output spin state configuration
   if (snapshot_output == true) {
		outputFilename = "Data-SpinState-L_" + L_string + "-Lz_" + Lz_string + "-disorder_filling_" + disorder_filling_string + "-Before.dat";
		outputFile.open(outputFilename);
		outputFile.setf(ios::scientific);

		for (int idx = 0; idx < N; idx++) {
			outputFile << idx << " " << spins[idx] << endl;
		}
		outputFile.close();   
	}
	if ( idx_disorder == 0 ) {
		t = clock();
	}

	// Perform Metropolis step in parallel for all 2*T_steps ensembles

	// Measure energy each ensemble: spins_energy[], 2*T_steps ensembles
	#pragma omp parallel for 
	for (int idx_ensemble = 0; idx_ensemble < 2*T_steps; idx_ensemble++) {
		spins_energy[idx_ensemble] = 0.; // initialize to zero
		double en = 0.; 
		for (int i = 0; i < counter_filled; i++) {
			int site = filled_sites[i];
			int spin1 = spins[idx_ensemble * N + site];
			if (abs(spin1) == 1) { // spin[i] is core spin Co^(4+)
				for (int j = 0; j < num_neighbors_direct; j++) { // loop over nearest-neighbors of site "site"
					int neighbor_site = neighbors_direct[num_neighbors_direct*site + j];
					if (abs(spins[idx_ensemble * N + neighbor_site]) == 1) { // neighbor site is Co^(4+) spin 1/2
						en += J44*static_cast<double>(spin1*spins[idx_ensemble * N + neighbor_site]);
					}
					else if (abs(spins[idx_ensemble * N + neighbor_site]) == 2) { 
						en += J34*static_cast<double>(spin1*spins[idx_ensemble * N + neighbor_site])/2.; // divide by 2 to account for fact that we take all spins to be of size 1 (but store Co^(3+) spins as \pm 2)
					}
					else {
						// do nothing as neighboring site is empty, then spins[n] = 0
					}
				} // for j over neighboring sites
			} // if spin = core spin at given site 
			else if (abs(spin1) == 2) { // spin[i] is satellite spin Co^(3+)
				for (int j = 0; j < num_neighbors_direct; j++) { // loop over nearest-neighbors of site "site"
					int neighbor_site = neighbors_direct[num_neighbors_direct*site + j];
					if (abs(spins[idx_ensemble * N + neighbor_site]) == 1) { // neighbor site is Co^(4+) spin 1/2
						en += J34*static_cast<double>(spin1*spins[idx_ensemble * N + neighbor_site])/2.; // divide by 2 as spin1 is \pm 2
					}
					else if (abs(spins[idx_ensemble * N + neighbor_site]) == 2) { // this includes if neighboring site is empty, then spins[n] = 0
						en += J33*static_cast<double>(spin1*spins[idx_ensemble * N + neighbor_site])/4.; // divide by 4 as both spins are Co^(3+) and thus stored as \pm 2
					}
					else {
						// do nothing as neighboring site is empty
					}
				} // for j over neighboring sites
			}
			else { 
				// do nothing if no spin at site i. This should not occur as we only loop over filled sites. 
			}
  		} // for i over all filled sites

  		en = en/2.; // we have summed over all bonds twice above, so we must divide by 2 here
		spins_energy[idx_ensemble] = en; // store energy of ensemble idx_ensemble into array spins_energy
	} // for idx_ensemble

	// Print energy of each ensemble
	/*for (int idx = 0; idx < 2*T_steps; idx++) {
		cout << "Energy of ensemble[" << idx << "] = " << spins_energy[idx] << endl;
	}*/

//	cout << endl << "Thermalization started." << endl << endl;
	// Thermalization
	for (unsigned long int idx_mc_bin = 0; idx_mc_bin < thermalization_bins; idx_mc_bin++) {
		for (unsigned long int idx_mcs = 0; idx_mcs < mc_binsize; idx_mcs++) {
			// Metropolis step for each ensemble
			#pragma omp parallel for 
			for (int idx_ensemble = 0; idx_ensemble < 2*T_steps; idx_ensemble++) {
				for (int idx = 0; idx < N; idx++) { // loop over all N lattice sites within each MC step
					int j = dis(mt_rand_mc_vector[idx_ensemble])*counter_filled; // randomly select site which has non-zero spin (N_max = disorder_filling*N such sites exist): j \in [0, N_max-1]
					int site = filled_sites[j]; // site index of this site is determined by the order of which sites are (randomly) filled
					int spin1 = spins[idx_ensemble*N + site]; // spin at selected site j
					int a = 0; // sum over spin 1 neighbors \in {-6, -5, ..., 6}
					int b = 0; // sum over spin 1/2 neighbors \in {-6, -5, ..., 6}

					double flip_prob = 0.;
					// Determine flip probability by counting number of nearest-neighbor Co^{4+} and Co^{3+} spins
					if (abs(spin1) == 1) { // if 4+ site -> apply prob_1 table of flip probabilities
						int spin_idx = (spin1 + 1)/2; // = {0, 1}
						//cout << "spin_idx = " << spin_idx << endl;
						for (int jdx = 0; jdx < num_neighbors_direct; jdx++) { // sum over six nearest neighbors
							int neighbor_site = neighbors_direct[num_neighbors_direct*site + jdx];
							//cout << "spin_neighbor[" << neighbor_site << "] = " << spin[neighbor_site] << endl;
							if (abs(spins[idx_ensemble*N + neighbor_site]) == 1) { // neighbor site is spin 1/2
								b += spins[idx_ensemble*N + neighbor_site];
							}
							else if (abs(spins[idx_ensemble*N + neighbor_site]) == 2) { // else: neighboring site is spin 1 (spin[n] = \pm 2) or empty (spin[n] = 0)
								a += spins[idx_ensemble*N + neighbor_site]/2; // a is the sum of S=\tilde{S}/2 \in {-6, -5, ..., 6}
							}
							else {
								// do nothing as neighboring site is empty
							}
						} // for jdx over neighbors
						// Co^{4+} flip probability from table prob_1
						// cout << "spin1_initial = " << spin1 <<  ", b = " << b << ", a = " << a << ", spin_idx = " << spin_idx << endl;
						// determine temperature of this ensemble
						int T_index = pt_EtoT[idx_ensemble];
						flip_prob = prob_1[spin_idx + (a + num_neighbors_direct)*2 + (b + num_neighbors_direct)*2*(2*num_neighbors_direct + 1) + T_index*2*(2*num_neighbors_direct + 1) * (2*num_neighbors_direct + 1)];	
					} // if core site = 4+ (spin 1/2)
					else if (abs(spin1) == 2) { // if 3+ site -> apply prob_2 table of flip probabilities
						int spin_idx = (spin1/2 + 1)/2; // = {0, 1}
						for (int jdx = 0; jdx < num_neighbors_direct; jdx++) { // sum over six nearest neighbors
							int neighbor_site = neighbors_direct[num_neighbors_direct*site + jdx];
							if (abs(spins[idx_ensemble*N + neighbor_site]) == 1) { // neighbor site is spin 1/2
								b += spins[idx_ensemble*N + neighbor_site];
							}
							else if (abs(spins[idx_ensemble*N + neighbor_site]) == 2) { // this includes if neighboring site is empty, then spin[n] = 0
								a += spins[idx_ensemble*N + neighbor_site]/2; // a is the sum of S=\tilde{S}/2 \in {-6, -5, ..., 6}
							}
							else {
								// do nothing as neighboring site is empty
							}
						} // for jdx over neighbors
						// cout << "b = " << b << ", a = " << a << ", spin_idx = " << spin_idx << endl;
						int T_index = pt_EtoT[idx_ensemble];
						flip_prob = prob_2[spin_idx + (a + num_neighbors_direct)*2 + (b + num_neighbors_direct)*2*(2*num_neighbors_direct + 1)+ T_index*2*(2*num_neighbors_direct + 1) * (2*num_neighbors_direct + 1)];
					} // if core site = 3+ (spin 1)
					double prand = dis(mt_rand_mc_vector[idx_ensemble]); // draw random number to compare with flip_prob
					if (prand < flip_prob) {
						// flip spin
						spins[idx_ensemble*N + site] = - spin1;
						// update energy due to spin flip. Here, we assume all spins are spin 1
						if (abs(spin1) == 1) {
							spins_energy[idx_ensemble] += -2.*spin1*(J44*b + J34*a); // Delta E = -2*S_i * (local_field)
						}
						else if (abs(spin1) == 2) {
							spins_energy[idx_ensemble] += -spin1*(J34*b + J33*a); // Delta E = -2*S_i * (local_field). Factor of 2 is already in spin1 
						}
						if (pt_EtoT[idx_ensemble] == 0) {
							acceptance_number_Metropolis++; // Metropolis acceptance rate at lowest T
						} // increase acceptance prob by one
					} // if prand < flip_prob, the flip spin
				} // for idx < N: one Metropolis step
			} // omp for idx_ensemble
	
			// Print energy of each ensemble
			/*for (int idx = 0; idx < 2*T_steps; idx++) {
				cout << "Energy of ensemble[" << idx << "] = " << spins_energy[idx] << endl;
			}*/
			if (parallel_tempering_step == true) {
			// Parallel Tempering step 
			#pragma omp parallel for
			for (int idx_replica = 0; idx_replica <= 1; idx_replica++) {
				#pragma omp parallel for
				for (int idx_T = 0; idx_T < T_steps/2 - static_cast<int>(swap_even_pairs); idx_T++) {
					double delta_en = spins_energy[idx_replica*T_steps + pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs)]] - spins_energy[idx_replica*T_steps + pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs) + 1]];
					//cout << "Pairs of temperatures = (" << 2*idx_T + static_cast<int>(swap_even_pairs) << " = " << temperatures[2*idx_T + static_cast<int>(swap_even_pairs) ] << ", " << 2*idx_T + static_cast<int>(swap_even_pairs) + 1 << " = " << temperatures[2*idx_T + static_cast<int>(swap_even_pairs) + 1] << "). Attempt to exchange temperatures." << endl;
					//cout << "pt_TtoE[" << idx_replica*T_steps + 2*idx_T << "] = " << pt_TtoE[idx_replica*T_steps + 2*idx_T] << ", pt_TtoE[" << idx_replica*T_steps + 2*idx_T + 1 << "] = " << pt_TtoE[idx_replica*T_steps +2*idx_T + 1] << endl;
					//cout << "spins_energy[" << idx_replica*T_steps + pt_TtoE[idx_replica*T_steps +2*idx_T] << "] = " << spins_energy[idx_replica*T_steps + pt_TtoE[idx_replica*T_steps +2*idx_T]] << ", spins_energy[" << idx_replica*T_steps + pt_TtoE[idx_replica*T_steps +2*idx_T+1] << "] = " << spins_energy[idx_replica*T_steps + pt_TtoE[idx_replica*T_steps +2*idx_T+1]] << endl;
					//cout << "delta_en of this pair = " << delta_en << endl;
					if (delta_en > 0.) {
						// exchange for sure
						//cout << "Swap temperatures as positive delta_en = " << delta_en << endl;
						int a = pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs)];
						int b = pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs) + 1];
						pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs)] = b;
						pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs) + 1] = a;
						pt_EtoT[idx_replica*T_steps + a] = 2*idx_T + static_cast<int>(swap_even_pairs) + 1;
						pt_EtoT[idx_replica*T_steps + b] = 2*idx_T + static_cast<int>(swap_even_pairs);
						acceptance_number_pt++;
					}
					else {
						double exchange_prob = exp((1./temperatures[2*idx_T + static_cast<int>(swap_even_pairs)] - 1./temperatures[2*idx_T + static_cast<int>(swap_even_pairs) + 1])*delta_en);
						//cout << "exchange_prob = " << exchange_prob << endl;
						double prand = dis(mt_rand_mc_vector[idx_replica*T_steps + pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs)]]); // draw random number to compare with flip_prob
						if (prand < exchange_prob) { // exchange
							//cout << "Swap temperatures even though delta_en is negative as prand =  " << prand << " < exchange_prob = " << exchange_prob << endl;
							int a = pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs)];
							int b = pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs) + 1];
							pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs)] = b;
							pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs) + 1] = a;
							pt_EtoT[idx_replica*T_steps + a] = 2*idx_T + static_cast<int>(swap_even_pairs) + 1;
							pt_EtoT[idx_replica*T_steps + b] = 2*idx_T + static_cast<int>(swap_even_pairs);
							acceptance_number_pt++;
						}
					}
				} // for idx_T
			} // for idx_replica: end of Parallel Tempering step
		} // if Parallel Tempering == true			
		swap_even_pairs = !swap_even_pairs;
		//cout << "swap_even_pairs = " << swap_even_pairs << endl;
		} // for idx_mcs
	} // for idx_mc_bin (end of thermalization)

/*
	for (int idx = 0; idx < 2*T_steps; idx++) {
		cout << "pt_TtoE[" << idx << "] = " << pt_TtoE[idx] << ", pt_EtoT[" << idx << "] = " << pt_EtoT[idx] << endl;
	}
*/

//	cout << endl << "Thermalization finished. Measurements begin. " << endl << endl;
/***************************************************************************************
 * ***** MC simulation including measurements in each MC bin ***************************
***************************************************************************************/
	/*for (int idx = 0; idx < 2*N; idx++) {
		cout << "spins[" << idx << "] = " << spins[idx] << endl;
	}
	*/
	for (unsigned long int idx_mc_bin = 0; idx_mc_bin < mc_bins; idx_mc_bin++) {
		for (unsigned long int idx_mcs = 0; idx_mcs < mc_binsize; idx_mcs++) {
			// Metropolis step for each ensemble
			#pragma omp parallel for 
			for (int idx_ensemble = 0; idx_ensemble < 2*T_steps; idx_ensemble++) {
				//cout << endl;
				//cout << "idx_ensemble = " << idx_ensemble << endl;
				//cout << endl;
				//cout << "I am thread number " << omp_get_thread_num() << " working on idx_ensemble = " << idx_ensemble << endl;
				int T_index = pt_EtoT[idx_ensemble];
				for (int idx = 0; idx < N; idx++) { // loop over all N lattice sites within each MC step
					int j = dis(mt_rand_mc_vector[idx_ensemble])*counter_filled; // randomly select site which has non-zero spin (N_max = disorder_filling*N such sites exist): j \in [0, N_max-1]
					//omp_get_thread_num()
					int site = filled_sites[j]; // site index of this site is determined by the order of which sites are (randomly) filled
					int spin1 = spins[idx_ensemble*N + site]; // spin at selected site j
					/*if (idx < 10) {
						cout << "j = " << j << ", site = " << site << ", spin1 = " << spin1 << endl;
					}*/
					int a = 0; // sum over spin 1 neighbors \in {-6, -5, ..., 6}
					int b = 0; // sum over spin 1/2 neighbors \in {-6, -5, ..., 6}

					double flip_prob = 0.;
					// Determine flip probability by counting number of nearest-neighbor Co^{4+} and Co^{3+} spins
					if (abs(spin1) == 1) { // if 4+ site -> apply prob_1 table of flip probabilities
						int spin_idx = (spin1 + 1)/2; // = {0, 1}
						/*if (idx < 10) {
							cout << "spin_idx at site 4+ = " << spin_idx << endl;
						}*/
						for (int jdx = 0; jdx < num_neighbors_direct; jdx++) { // sum over six nearest neighbors
							int neighbor_site = neighbors_direct[num_neighbors_direct*site + jdx];
							/*if (idx < 10) {
								cout << "spin_neighbor[" << neighbor_site << "] = " << spins[idx_ensemble*N + neighbor_site] << endl;	
							}*/
							if (abs(spins[idx_ensemble*N + neighbor_site]) == 1) { // neighbor site is spin 1/2
								b += spins[idx_ensemble*N + neighbor_site];
							}
							else if (abs(spins[idx_ensemble*N + neighbor_site]) == 2) { // else: neighboring site is spin 1 (spin[n] = \pm 2) or empty (spin[n] = 0)
								a += spins[idx_ensemble*N + neighbor_site]/2; // a is the sum of S=\tilde{S}/2 \in {-6, -5, ..., 6}
							}
						} // for jdx over neighbors
						// Co^{4+} flip probability from table prob_1
						/*if (idx < 10) {
							cout << "spin1_initial = " << spin1 <<  ", b = " << b << ", a = " << a << ", spin_idx = " << spin_idx << endl;
						}*/
						// determine temperature of this ensemble
						flip_prob = prob_1[spin_idx + (a + num_neighbors_direct)*2 + (b + num_neighbors_direct)*2*(2*num_neighbors_direct + 1) + T_index*2*(2*num_neighbors_direct + 1) * (2*num_neighbors_direct + 1)];	
					} // if core site = 4+ (spin 1/2)
					else if (abs(spin1) == 2) { // if 3+ site -> apply prob_2 table of flip probabilities
						int spin_idx = (spin1/2 + 1)/2; // = {0, 1}
						/*if (idx < 10) {
							cout << "spin_idx at site 3+ = " << spin_idx << endl;
						}*/
						for (int jdx = 0; jdx < num_neighbors_direct; jdx++) { // sum over six nearest neighbors
							int neighbor_site = neighbors_direct[num_neighbors_direct*site + jdx];
							/*if (idx < 10) {
								cout << "spin_neighbor[" << neighbor_site << "] = " << spins[idx_ensemble*N + neighbor_site] << endl;	
							}*/
							if (abs(spins[idx_ensemble*N + neighbor_site]) == 1) { // neighbor site is spin 1/2
								b += spins[idx_ensemble*N + neighbor_site];
							}
							else if (abs(spins[idx_ensemble*N + neighbor_site]) == 2) { // this includes if neighboring site is empty, then spin[n] = 0
								a += spins[idx_ensemble*N + neighbor_site]/2; // a is the sum of S=\tilde{S}/2 \in {-6, -5, ..., 6}
							}
						} // for jdx over neighbors
						/*if (idx < 10) {
							cout << "spin1_initial = " << spin1 <<  ", b = " << b << ", a = " << a << ", spin_idx = " << spin_idx << endl;
						}*/
						flip_prob = prob_2[spin_idx + (a + num_neighbors_direct)*2 + (b + num_neighbors_direct)*2*(2*num_neighbors_direct + 1)+ T_index*2*(2*num_neighbors_direct + 1) * (2*num_neighbors_direct + 1)];
						/*if (idx < 10) {
							cout << "flip_prob = prob_2 = " << flip_prob << endl;
						}*/
					} // if core site = 3+ (spin 1)
					double prand = dis(mt_rand_mc_vector[idx_ensemble]); // draw random number to compare with flip_prob
					if (prand < flip_prob) {
						// flip spin
						spins[idx_ensemble*N + site] = - spin1;
						/*if (idx < 10) {
							cout << "prand = " << prand << ", flip_prob = " << flip_prob << ". Thus spin is flipped and new spin entry = " << spins[idx_ensemble*N + site] << endl;
						}*/
						// update energy due to spin flip. Here, we assume all spins are spin 1
						if (abs(spin1) == 1) {
							spins_energy[idx_ensemble] += -2.*spin1*(J44*b + J34*a); // Delta E = -2*S_i * (local_field)
						}
						else if (abs(spin1) == 2) {
							spins_energy[idx_ensemble] += -spin1*(J34*b + J33*a); // Delta E = -2*S_i * (local_field). Factor of 2 is already in spin1 
						}
						if (pt_EtoT[idx_ensemble] == 0) {
							acceptance_number_Metropolis++; // Metropolis acceptance rate at lowest T
						} // increase acceptance prob by one
					} // if prand < flip_prob, the flip spin
					else {
						// do not flip spin, so do nothing. 
					}
				} // for idx < N: one Metropolis step
			} // omp for idx_ensemble (Metropolis step for each ensemble)

		if (parallel_tempering_step == true) {
			// Parallel Tempering step 
			#pragma omp parallel for 
			for (int idx_replica = 0; idx_replica <= 1; idx_replica++) {
				#pragma omp parallel for
				for (int idx_T = 0; idx_T < T_steps/2 - static_cast<int>(swap_even_pairs); idx_T++) {
					double delta_en = spins_energy[idx_replica*T_steps + pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs)]] - spins_energy[idx_replica*T_steps + pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs) + 1]];
					//cout << "Pairs of temperatures = (" << 2*idx_T + static_cast<int>(swap_even_pairs) << " = " << temperatures[2*idx_T + static_cast<int>(swap_even_pairs) ] << ", " << 2*idx_T + static_cast<int>(swap_even_pairs) + 1 << " = " << temperatures[2*idx_T + static_cast<int>(swap_even_pairs) + 1] << "). Attempt to exchange temperatures." << endl;
					//cout << "pt_TtoE[" << idx_replica*T_steps + 2*idx_T << "] = " << pt_TtoE[idx_replica*T_steps + 2*idx_T] << ", pt_TtoE[" << idx_replica*T_steps + 2*idx_T + 1 << "] = " << pt_TtoE[idx_replica*T_steps +2*idx_T + 1] << endl;
					//cout << "spins_energy[" << idx_replica*T_steps + pt_TtoE[idx_replica*T_steps +2*idx_T] << "] = " << spins_energy[idx_replica*T_steps + pt_TtoE[idx_replica*T_steps +2*idx_T]] << ", spins_energy[" << idx_replica*T_steps + pt_TtoE[idx_replica*T_steps +2*idx_T+1] << "] = " << spins_energy[idx_replica*T_steps + pt_TtoE[idx_replica*T_steps +2*idx_T+1]] << endl;
					//cout << "delta_en of this pair = " << delta_en << endl;
					if (delta_en > 0.) {
						// exchange for sure
						//cout << "Swap temperatures as positive delta_en = " << delta_en << endl;
						int a = pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs)];
						int b = pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs) + 1];
						pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs)] = b;
						pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs) + 1] = a;
						pt_EtoT[idx_replica*T_steps + a] = 2*idx_T + static_cast<int>(swap_even_pairs) + 1;
						pt_EtoT[idx_replica*T_steps + b] = 2*idx_T + static_cast<int>(swap_even_pairs);
						acceptance_number_pt++;
					}
					else {
						double exchange_prob = exp((1./temperatures[2*idx_T + static_cast<int>(swap_even_pairs)] - 1./temperatures[2*idx_T + static_cast<int>(swap_even_pairs) + 1])*delta_en);
						//cout << "exchange_prob = " << exchange_prob << endl;
						double prand = dis(mt_rand_mc_vector[idx_replica*T_steps + pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs)]]); // draw random number to compare with flip_prob
						if (prand < exchange_prob) { // exchange
							//cout << "Swap temperatures even though delta_en is negative as prand =  " << prand << " < exchange_prob = " << exchange_prob << endl;
							int a = pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs)];
							int b = pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs) + 1];
							pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs)] = b;
							pt_TtoE[idx_replica*T_steps + 2*idx_T + static_cast<int>(swap_even_pairs) + 1] = a;
							pt_EtoT[idx_replica*T_steps + a] = 2*idx_T + static_cast<int>(swap_even_pairs) + 1;
							pt_EtoT[idx_replica*T_steps + b] = 2*idx_T + static_cast<int>(swap_even_pairs);
							acceptance_number_pt++;
						}
					}
				} // for idx_T
			} // for idx_replica: end of parallel tempering step
		} // end if parallel tempering step

/*			
for (int idx = 0; idx < 2*T_steps; idx++) {
				cout << "pt_TtoE[" << idx << "] = " << pt_TtoE[idx] << ", pt_EtoT[" << idx << "] = " << pt_EtoT[idx] << endl;
			}
*/
			swap_even_pairs = !swap_even_pairs;
			//cout << "swap even pairs = " << static_cast<int>(swap_even_pairs) << endl;
		} // for idx_mcs (within one bin)

		/*******************************************************************************************
		***** Measurement of observables at the end of each mc_bin  ********************************
		*******************************************************************************************/
		// Magnetization
		#pragma omp parallel for
		for (int idx_T = 0; idx_T < T_steps; idx_T++) {
			#pragma omp parallel for
			for (int idx_replica = 0; idx_replica <= 1; idx_replica++) {
				int mag1 = 0; // all spins 1
				int mag2 = 0; // spin 3/2 at core Co^(4+) sites and spin 1 at satellite Co^(3+) sites
				int magAF = 0; // AF magnetization for all spins 1
				int magReplica = 0; // Equal to q in literature. Replica magnetization for all spins 1. 
				double *spins_data; // pointer to spins_data that will be used to compute Fourier transform
				double *speq; // pointer to Nyquist component (used in real-space Fourier transform program rlft3 of NR)
				spins_data = new double[N]; // use for measuring spins_q. Fill with values of array spins below and then compute Fourier transform. 
				speq = new double[2*L*Lz]; // will store Nyquist critical frequency of third=z component (see p. 634 of NR book) 
				for (int idx = 0; idx < N; idx++) {
					spins_data[idx] = 0.;
				}

				for (int idx = 0; idx < 2*L*Lz; idx++) {
					speq[idx] = 0.;
				}

				// determine ensemble that sits at temperature idx_T
				int ensemble = pt_TtoE[idx_replica*T_steps + idx_T];
				int ensemble_replica = ensemble; // initialize (will be modified in next line). Only used for computation of replica correlations, i.e., when idx_replica == 0
				if (idx_replica == 0) {
					ensemble_replica = pt_TtoE[T_steps + idx_T]; // determine ensemble in second set of replicas (idx_replica = 1). Needed to compute replica cross correlations. Here, idx_replica == 1 in the brackets. This is the ensemble at temperature idx_T in the second set of simulated systems (T_steps, ..., T_steps + idx_T, ..., 2*T_steps - 1) 
				}

				//cout << "ensemble that sits at temperature T[" << idx_T << "] = " << temperatures[idx_T] << " is ensemble = " << ensemble << endl;

				// Measure mag1, mag2, magAF, magReplica
				for (int idx = 0; idx < N; idx++) {
					int odd = is_odd(idx%L) + is_odd((idx/L)%L) + is_odd((idx/(L*L))%Lz); // check if site is odd. Then gets a minus sign for AF
					//cout << "idx = " << idx << ", odd = " << odd << endl;
					if (abs(spins[(idx_replica*T_steps + ensemble)*N + idx]) == 1) { // core site
						mag1 += 2*spins[(idx_replica*T_steps + ensemble)*N + idx]; // put spin 1 (divide by 2 in the end when computing magnetization)
						mag2 += 2*spins[(idx_replica*T_steps + ensemble)*N + idx]; // put spin 1 (divide by 2 in the end when computing magnetization)
						if (idx_replica == 0) {
							magReplica += 4*spins[ensemble*N + idx]*spins[(T_steps + ensemble_replica)*N + idx]; // take into account that different ensembles in the two replica sets can have different temperatures. 

						}
						if (odd == 0 || odd == 2) {
							magAF += 2*spins[(idx_replica*T_steps + ensemble)*N + idx];	// cos(pi*x)*cos(pi*y)*cos(pi*z) = +1
						} 
						else {
							magAF -= 2*spins[(idx_replica*T_steps + ensemble)*N + idx]; // cos(pi*x)*cos(pi*y)*cos(pi*z) = -1
						}
						spins_data[idx] = static_cast<double>(spins[(idx_replica*T_steps + ensemble)*N + idx]);
					}
					else { // satellite sites and empty sites
						mag1 += spins[(idx_replica*T_steps + ensemble)*N + idx]; // put spin 1 (divide by 2 in the end when computing magnetization)
						mag2 += spins[(idx_replica*T_steps + ensemble)*N + idx]; // put spin 1 (divide by 2 in the end when computing magnetization)
						if (idx_replica == 0) {
							magReplica += spins[ensemble*N + idx]*spins[(T_steps + ensemble_replica)*N + idx];
						}
						if (odd == 0 || odd == 2) {
							magAF += spins[(idx_replica*T_steps + ensemble)*N + idx]; // cos(pi*x)*cos(pi*y)*cos(pi*z) = +1
						} 
						else {
							magAF -= spins[(idx_replica*T_steps + ensemble)*N + idx]; // cos(pi*x)*cos(pi*y)*cos(pi*z) = -1
						}
						spins_data[idx] = static_cast<double>(spins[(idx_replica*T_steps + ensemble)*N + idx]/2); // spins_data contains \pm 1 (all spins are equal. Should only affect normalization of s_q.)
					} // if spin == 1 or else 
				} // for idx < N (running over all lattice sites)

				// Magnetization and Magnetization^2
				magnetization1[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] = static_cast<double>(mag1)/2./static_cast<double>(N); // normalize with respect to lattice sites (including empty ones). First /2. is to consider spin 1, second /2. is since we sum over two replicas. 
				//cout << "mag1 = " << static_cast<double>(mag1)/2./static_cast<double>(N) << endl;
				magnetization2[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] = static_cast<double>(mag2)/2./static_cast<double>(counter_filled); // normalize with respect to total number of spins
				magnetizationAF[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] = static_cast<double>(magAF)/2./static_cast<double>(N);

				magnetization1_squared[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] = magnetization1[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica]*magnetization1[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica];
				magnetization2_squared[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] = magnetization2[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica]*magnetization2[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica];
				magnetizationAF_squared[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] = magnetizationAF[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica]*magnetizationAF[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica];

				magnetization1_four[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] = magnetization1[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica]*magnetization1[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica]*magnetization1[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica]*magnetization1[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica];
				magnetizationAF_four[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] = magnetizationAF[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica]*magnetizationAF[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica]*magnetizationAF[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica]*magnetizationAF[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica];

				energy[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] = static_cast<double>(spins_energy[idx_replica*T_steps + ensemble]); 
				energy_squared[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] = energy[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica]*energy[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica];

				if (idx_replica == 0) { // save replica magnetization functions
					magnetization_replica[idx_T*mc_bins + idx_mc_bin] = static_cast<double>(magReplica)/4./static_cast<double>(N); // divide by four to consider spin 1 everywhere
					magnetization_replica_squared[idx_T*mc_bins + idx_mc_bin] = magnetization_replica[idx_T*mc_bins + idx_mc_bin]*magnetization_replica[idx_T*mc_bins + idx_mc_bin];
					magnetization_replica_four[idx_T*mc_bins + idx_mc_bin] = magnetization_replica[idx_T*mc_bins + idx_mc_bin]*magnetization_replica[idx_T*mc_bins + idx_mc_bin]*magnetization_replica[idx_T*mc_bins + idx_mc_bin]*magnetization_replica[idx_T*mc_bins + idx_mc_bin];
				
					if (output_s_q == true) {
					// Spectral function S(q), (ferromagnetic) correlation length \xi

					// perform 3D FFT on spins_data and store result in speq. This was done using function rlft3 from nr.
					// rlft3(spins_data, speq, 1, Lz, L, L);
	
					for (int zdx = 0; zdx < Lz; zdx++) {
						for (int ydx = 0; ydx < L; ydx++) {
							for (int xdx = 0; xdx < L/2; xdx++) {
								s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1) + xdx + ydx*(L/2+1) + zdx*(L/2+1)*L] =  (spins_data[2*xdx + ydx*L + zdx*L*L]*spins_data[2*xdx + ydx*L + zdx*L*L] + spins_data[2*xdx+1 + ydx*L + zdx*L*L]*spins_data[2*xdx+1 + ydx*L + zdx*L*L])/N/2.; // divide by two due to averaging over two replicas
							}
							s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1) + L/2 + ydx*(L/2+1) + zdx*(L/2+1)*L] = (speq[2*ydx + zdx*2*L]*speq[2*ydx + zdx*2*L] + speq[1 + 2*ydx + zdx*2*L]*speq[1 + 2*ydx + zdx*2*L])/N/2.; // divide by two due to averaging over two replicas
						}
					}

					if ( (s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1)]/s_q[idx_mc_bin*L*Lz*(L/2+1) + 1] - 1. > 0.) && (abs(s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1) + 1]) > 1.e-10 ) ) {
						xi_a[idx_T*mc_bins + idx_mc_bin] = 1/2./M_PI*sqrt(s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1)]/s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1) + 1] - 1.);
					}
					if ( ( (s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1) + 1]/s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1) + 2] - 1.)/(4. - s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1) + 1]/s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1) + 2]) > 0.) && (abs(s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1) + 2]) > 1.e-10 ) && ( abs( s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1)] ) > 1.e-10 ) ) {
						xi_b[idx_T*mc_bins + idx_mc_bin] = 1/2./M_PI*sqrt((s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1) + 1]/s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1) + 2] - 1.)/(4. - s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1) + 1]/s_q[idx_T*mc_bins*L*Lz*(L/2+1) + idx_mc_bin*L*Lz*(L/2+1) + 2]));
					}
					}
				} // if idx_replica == 0 

				delete[] spins_data;
				delete[] speq;				
			} // for idx_replica
		} // for idx_T
	} // for idx_mc_bin

	// Output final spin state after MC simulation 
	// Output snapshot of spin state
// Output spin state configuration
   if (snapshot_output == true) {
   outputFilename = "Data-SpinState-L_" + L_string + "-Lz_" + Lz_string + "-disorder_filling_" + disorder_filling_string + "-After.dat";
   outputFile.open(outputFilename);
   outputFile.setf(ios::scientific);

   for (int idx = 0; idx < N; idx++) {
   	outputFile << idx << " " << spins[idx] << endl;
   }
   outputFile.close();   
}


	// Output data for all bins for given idx_disorder into file
	/*
	outputFilename = "Data-Observables-AllBins-L_" + L_string + "-Lz_" + Lz_string + "-disorder_filling_" + disorder_filling_string + "-J33_" + J33_string + "-J34_" + J34_string + "-J44_" + J44_string + "-TMin_" + T_min_string + "-TMax_" + T_max_string + "-TSteps_" + T_steps_string + ".dat";
	outputFile.open(outputFilename, std::ios_base::app);

	outputFile.setf(ios::scientific);

	for (int idx_T = 0; idx_T < T_steps; idx_T++) {
		for (int idx_disorder = 0; idx_disorder < disorder_steps; idx_disorder++) {
			for (unsigned long int idx_mc_bin = 0; idx_mc_bin < mc_bins; idx_mc_bin++) {
				for (int idx_replica = 0; idx_replica <= 1; idx_replica++) {
					outputFile << magnetization1[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] << " " << magnetization2[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] << " " << magnetizationAF[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] << " " << energy[idx_T*2*mc_bins + 2*idx_mc_bin + idx_replica] << " ";
					if (idx_replica == 0) { 
						outputFile << magnetization_replica[idx_T*mc_bins + idx_mc_bin] << " " << endl;
					}
					else {
						outputFile << "0." << " " << endl;
					}
				}
			}
		}
	} // for idx_T < T_steps
 	outputFile.close(); 
	*/

	// Error analysis
	// Variables needed for jackknife method of computing error sigma. They will be averages using all bins but one. 
	vector<double> mag1(T_steps*2*mc_bins, 0.);
	vector<double> mag2(T_steps*2*mc_bins, 0.);
	vector<double> magAF(T_steps*2*mc_bins, 0.);
	vector<double> magReplica(T_steps*mc_bins, 0.);
	vector<double> mag1_squared(T_steps*2*mc_bins, 0.);
	vector<double> mag2_squared(T_steps*2*mc_bins, 0.);
	vector<double> magAF_squared(T_steps*2*mc_bins, 0.);
	vector<double> magReplica_squared(T_steps*mc_bins, 0.);
	vector<double> mag1_four(T_steps*2*mc_bins, 0.);
	vector<double> magAF_four(T_steps*2*mc_bins, 0.);
	vector<double> magReplica_four(T_steps*mc_bins, 0.);
	vector<double> chi1(T_steps*2*mc_bins, 0.);
	vector<double> chiAF(T_steps*2*mc_bins, 0.);
	vector<double> chiReplica(T_steps*mc_bins, 0.);
	vector<double> bind_ratio(T_steps*2*mc_bins, 0.);
	vector<double> bind_ratioAF(T_steps*2*mc_bins, 0.);
	vector<double> bind_ratioReplica(T_steps*mc_bins, 0.);

	vector<double> en(T_steps*2*mc_bins, 0.);
	vector<double> en_squared(T_steps*2*mc_bins, 0.);
	vector<double> cV(T_steps*2*mc_bins, 0.);

	vector<double> corr_a(T_steps*mc_bins, 0.);
	vector<double> corr_b(T_steps*mc_bins, 0.);

	// Store results in variables magnetization1_av_dis[T_steps*disorder_steps]
	#pragma omp parallel for
	for (int idx_T = 0; idx_T < T_steps; idx_T++) {
		#pragma omp parallel for
		for (int idx_replica = 0; idx_replica <= 1; idx_replica++) {
			for (unsigned long int idx = 0; idx < mc_bins; idx++) {
				magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += abs(magnetization1[idx_T*2*mc_bins + 2*idx + idx_replica])/static_cast<double>(mc_bins);
				magnetization2_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += abs(magnetization2[idx_T*2*mc_bins + 2*idx + idx_replica])/static_cast<double>(mc_bins);
				magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += abs(magnetizationAF[idx_T*2*mc_bins + 2*idx + idx_replica])/static_cast<double>(mc_bins);

				magnetization1_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += magnetization1_squared[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins);
				magnetization2_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += magnetization2_squared[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins);
				magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += magnetizationAF_squared[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins);
			
				magnetization1_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += magnetization1_four[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins);
				magnetizationAF_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += magnetizationAF_four[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins);
		
				energy_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += energy[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins);
				energy_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += energy_squared[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins);	
			} // for idx < mc_bins
		
			// Magnetic susceptibility from <m^2> and <m>
			chi_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = static_cast<double>(N)/temperatures[idx_T]*(magnetization1_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] - magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]*magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
			// AF magnetic susceptibility from <m_AF^2> and <m_AF>
			chi_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = static_cast<double>(N)/temperatures[idx_T]*(magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] - magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]*magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
				
			binder_ratio_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = magnetization1_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]/(magnetization1_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]*magnetization1_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
			if (abs(magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]) > 1.e-5 ) {
				binder_ratioAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = magnetizationAF_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]/(magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]*magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
			} 
			else 
			{ 
				binder_ratioAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = magnetizationAF_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]*1.e10;
			}
			// Specific heat from <E^2> and <E>	
			specific_heat_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = 1./static_cast<double>(N)/temperatures[idx_T]/temperatures[idx_T]*(energy_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] - energy_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]*energy_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
			} // for idx_replica
		} // for idx_T

	#pragma omp parallel for
	for (int idx_T = 0; idx_T < T_steps; idx_T++) {
		for (unsigned long int idx = 0; idx < mc_bins; idx++) {
			magnetization_replica_av_dis[idx_T*disorder_steps + idx_disorder] += magnetization_replica[idx_T*mc_bins + idx]/static_cast<double>(mc_bins);
			magnetization_replica_squared_av_dis[idx_T*disorder_steps + idx_disorder] += magnetization_replica_squared[idx_T*mc_bins + idx]/static_cast<double>(mc_bins);
			magnetization_replica_four_av_dis[idx_T*disorder_steps + idx_disorder] += magnetization_replica_four[idx_T*mc_bins + idx]/static_cast<double>(mc_bins);
			//cout << "magnetization1[" << idx << "] = " << magnetization1[idx] << ", magnetization1_squared[" << idx << "] = " << magnetization1_squared[idx] << endl;
			//cout << "magnetization1_four[" << idx << "] = " << magnetization1_four[idx] << endl;

			// save each disorder realization of s_q (this is very big for large systems and disorder_steps)
			if (output_s_q == true) {
				for (int jdx = 0; jdx < L*Lz*(L/2+1); jdx++) {
					s_q_av_dis[idx_T*disorder_steps*L*Lz*(L/2+1) + idx_disorder*L*Lz*(L/2+1) + jdx] += s_q[idx_T*mc_bins + idx*L*Lz*(L/2+1) + jdx]/static_cast<double>(mc_bins);
				}
			}
			/*
			for (int jdx = 0; jdx < L*Lz*(L/2+1); jdx++) { 
				s_q_av[idx_T*L*Lz*(L/2+1) + jdx] = static_cast<double>(idx_disorder + 1)*s_q_av[idx_T*L*Lz*(L/2+1) + jdx] + s_q_av_dis_small[idx_T*disorder_steps + idx_disorder*L*Lz*(L/2+1) + jdx]
			}
			*/
			// s_q_av[idx_T*L*Lz*(L/2+1) + jdx] += s_q_av_dis[idx_T*disorder_steps*L*Lz*(L/2+1) + idx*L*Lz*(L/2+1) + jdx]/static_cast<double>(idx_disorder + 1)/2.;

			xi_a_av_dis[idx_T*disorder_steps + idx_disorder] += xi_a[idx_T*mc_bins + idx]/static_cast<double>(mc_bins);
			xi_b_av_dis[idx_T*disorder_steps + idx_disorder] += xi_b[idx_T*mc_bins + idx]/static_cast<double>(mc_bins);
		} // for idx < mc_bins
		// Replica magnetic susceptibility from <m_replica^2> and <m_replica>
		chi_magnetization_replica_av_dis[idx_T*disorder_steps + idx_disorder] = static_cast<double>(N)*(magnetization_replica_squared_av_dis[idx_T*disorder_steps + idx_disorder] - magnetization_replica_av_dis[idx_T*disorder_steps + idx_disorder]*magnetization_replica_av_dis[idx_T*disorder_steps + idx_disorder]);
		binder_ratio_replica_av_dis[idx_T*disorder_steps + idx_disorder] = magnetization_replica_four_av_dis[idx_T*disorder_steps + idx_disorder]/(magnetization_replica_squared_av_dis[idx_T*disorder_steps + idx_disorder]*magnetization_replica_squared_av_dis[idx_T*disorder_steps + idx_disorder]);
	} // for idx_T < T_steps
	
	// JACKNIFE ERROR ANALYSIS
	#pragma omp parallel for
	for (int idx_T = 0; idx_T < T_steps; idx_T++){
		for (unsigned long int idx = 0; idx < mc_bins; idx++) {
			for (int idx_replica = 0; idx_replica <= 1; idx_replica++) {
				for (unsigned long int jdx = 0; jdx < mc_bins; jdx++) {
					mag1[idx_T*2*mc_bins + 2*idx + idx_replica] += abs(magnetization1[idx_T*2*mc_bins + 2*jdx + idx_replica]);
					mag2[idx_T*2*mc_bins + 2*idx + idx_replica] += abs(magnetization2[idx_T*2*mc_bins + 2*jdx + idx_replica]);
					magAF[idx_T*2*mc_bins + 2*idx + idx_replica] += abs(magnetizationAF[idx_T*2*mc_bins + 2*jdx + idx_replica]);
					mag1_squared[idx_T*2*mc_bins + 2*idx + idx_replica] += magnetization1_squared[idx_T*2*mc_bins + 2*jdx + idx_replica];
					mag2_squared[idx_T*2*mc_bins + 2*idx + idx_replica] += magnetization2_squared[idx_T*2*mc_bins + 2*jdx + idx_replica];
					magAF_squared[idx_T*2*mc_bins + 2*idx + idx_replica] += magnetizationAF_squared[idx_T*2*mc_bins + 2*jdx + idx_replica];
					mag1_four[idx_T*2*mc_bins + 2*idx + idx_replica] += magnetization1_four[idx_T*2*mc_bins + 2*jdx + idx_replica];
					magAF_four[idx_T*2*mc_bins + 2*idx + idx_replica] += magnetizationAF_four[idx_T*2*mc_bins + 2*jdx + idx_replica];

					en[idx_T*2*mc_bins + 2*idx + idx_replica] += energy[idx_T*2*mc_bins + 2*jdx + idx_replica];
					en_squared[idx_T*2*mc_bins + 2*idx + idx_replica] += energy_squared[idx_T*2*mc_bins + 2*jdx + idx_replica];
				} // for jdx
				// Magnetization 1 and chi1
				mag1[idx_T*2*mc_bins + 2*idx + idx_replica] -= abs(magnetization1[idx_T*2*mc_bins + 2*idx + idx_replica]);
				mag1[idx_T*2*mc_bins + 2*idx + idx_replica] = mag1[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins - 1);
				mag1_squared[idx_T*2*mc_bins + 2*idx + idx_replica] -= magnetization1_squared[idx_T*2*mc_bins + 2*idx + idx_replica];
				mag1_squared[idx_T*2*mc_bins + 2*idx + idx_replica] = mag1_squared[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins - 1);

				chi1[idx_T*2*mc_bins + 2*idx + idx_replica] = static_cast<double>(N)/temperatures[idx_T]*(mag1_squared[idx_T*2*mc_bins + 2*idx + idx_replica] - mag1[idx_T*2*mc_bins + 2*idx + idx_replica]*mag1[idx_T*2*mc_bins + 2*idx + idx_replica]);

				mag1_four[idx_T*2*mc_bins + 2*idx + idx_replica] -= magnetization1_four[idx_T*2*mc_bins + 2*idx + idx_replica];
				mag1_four[idx_T*2*mc_bins + 2*idx + idx_replica] = mag1_four[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins - 1);

				bind_ratio[idx_T*2*mc_bins + 2*idx + idx_replica] = mag1_four[idx_T*2*mc_bins + 2*idx + idx_replica]/(mag1_squared[idx_T*2*mc_bins + 2*idx + idx_replica]*mag1_squared[idx_T*2*mc_bins + 2*idx + idx_replica]);

				// Magnetization 2
				mag2[idx_T*2*mc_bins + 2*idx + idx_replica] -= abs(magnetization2[idx_T*2*mc_bins + 2*idx + idx_replica]);
				mag2[idx_T*2*mc_bins + 2*idx + idx_replica] = mag2[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins - 1);
				mag2_squared[idx_T*2*mc_bins + 2*idx + idx_replica] -= magnetization2_squared[idx_T*2*mc_bins + 2*idx + idx_replica];
				mag2_squared[idx_T*2*mc_bins + 2*idx + idx_replica] = mag2_squared[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins - 1);

				// Magnetization AF
				magAF[idx_T*2*mc_bins + 2*idx + idx_replica] -= abs(magnetizationAF[idx_T*2*mc_bins + 2*idx + idx_replica]);
				magAF[idx_T*2*mc_bins + 2*idx + idx_replica] = magAF[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins - 1);
				magAF_squared[idx_T*2*mc_bins + 2*idx + idx_replica] -= magnetizationAF_squared[idx_T*2*mc_bins + 2*idx + idx_replica];
				magAF_squared[idx_T*2*mc_bins + 2*idx + idx_replica] = magAF_squared[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins - 1);

				chiAF[idx_T*2*mc_bins + 2*idx + idx_replica] = static_cast<double>(N)/temperatures[idx_T]*(magAF_squared[idx_T*2*mc_bins + 2*idx + idx_replica] - magAF[idx_T*2*mc_bins + 2*idx + idx_replica]*magAF[idx_T*2*mc_bins + 2*idx + idx_replica]);

				magAF_four[idx_T*2*mc_bins + 2*idx + idx_replica] -= magnetizationAF_four[idx_T*2*mc_bins + 2*idx + idx_replica];
				magAF_four[idx_T*2*mc_bins + 2*idx + idx_replica] = magAF_four[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins - 1);
				if (abs(magAF_squared[idx_T*2*mc_bins + 2*idx + idx_replica]) > 1.e-5 ) {
					bind_ratioAF[idx_T*2*mc_bins + 2*idx + idx_replica] = magAF_four[idx_T*2*mc_bins + 2*idx + idx_replica]/(magAF_squared[idx_T*2*mc_bins + 2*idx + idx_replica]*magAF_squared[idx_T*2*mc_bins + 2*idx + idx_replica]);
				}
				else {
					bind_ratioAF[idx_T*2*mc_bins + 2*idx + idx_replica] = 1.e10*magAF_four[idx_T*2*mc_bins + 2*idx + idx_replica];	
				}
				// Energy and specific heat
				en[idx_T*2*mc_bins + 2*idx + idx_replica] -= energy[idx_T*2*mc_bins + 2*idx + idx_replica];
				en[idx_T*2*mc_bins + 2*idx + idx_replica] = en[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins - 1);

				en_squared[idx_T*2*mc_bins + 2*idx + idx_replica] -= energy_squared[idx_T*2*mc_bins + 2*idx + idx_replica];
				en_squared[idx_T*2*mc_bins + 2*idx + idx_replica] = en_squared[idx_T*2*mc_bins + 2*idx + idx_replica]/static_cast<double>(mc_bins - 1);

				cV[idx_T*2*mc_bins + 2*idx + idx_replica] = 1./static_cast<double>(N)/temperatures[idx_T]/temperatures[idx_T]*(en_squared[idx_T*2*mc_bins + 2*idx + idx_replica] - en[idx_T*2*mc_bins + 2*idx + idx_replica]*en[idx_T*2*mc_bins + 2*idx + idx_replica]);
			} // idx_replica	
		} // for idx

		for (unsigned long int idx = 0; idx < mc_bins; idx++) { // quantities that have mc_bins (not 2*mc_bins)
			for (unsigned long int jdx = 0; jdx < mc_bins; jdx++) { 
				corr_a[idx_T*mc_bins + idx] += xi_a[idx_T*mc_bins + jdx];
				corr_b[idx_T*mc_bins + idx] += xi_b[idx_T*mc_bins + jdx];
			
				// Replica magnetization, chi, binder_ratio
				magReplica[idx_T*mc_bins + idx] += magnetization_replica[idx_T*mc_bins + jdx];
				magReplica_squared[idx_T*mc_bins + idx] += magnetization_replica_squared[idx_T*mc_bins + jdx];
				magReplica_four[idx_T*mc_bins + idx] += magnetization_replica_four[idx_T*mc_bins + jdx];
			} // for jdx
			corr_a[idx_T*mc_bins + idx] -= xi_a[idx];
			corr_a[idx_T*mc_bins + idx] = corr_a[idx_T*mc_bins + idx]/static_cast<double>(mc_bins - 1);
			corr_b[idx_T*mc_bins + idx] -= xi_b[idx];
			corr_b[idx_T*mc_bins + idx] = corr_b[idx_T*mc_bins + idx]/static_cast<double>(mc_bins - 1);
			
			magReplica[idx_T*mc_bins + idx] -= magnetization_replica[idx_T*mc_bins + idx];
			magReplica[idx_T*mc_bins + idx] /= static_cast<double>(mc_bins - 1);
			magReplica_squared[idx_T*mc_bins + idx] -= magnetization_replica_squared[idx_T*mc_bins + idx];
			magReplica_squared[idx_T*mc_bins + idx] /= static_cast<double>(mc_bins - 1);

			chiReplica[idx_T*mc_bins + idx] = static_cast<double>(N)*(magReplica_squared[idx_T*mc_bins + idx] - magReplica[idx_T*mc_bins + idx]*magReplica[idx_T*mc_bins + idx]);

			magReplica_four[idx_T*mc_bins + idx] -= magnetization_replica_four[idx_T*mc_bins + idx];
			magReplica_four[idx_T*mc_bins + idx] = magReplica_four[idx_T*mc_bins + idx]/static_cast<double>(mc_bins - 1);
			if (abs(magReplica_squared[idx_T*mc_bins + idx]) > 1.e-5 ) {
				bind_ratioReplica[idx_T*mc_bins + idx] = magReplica_four[idx_T*mc_bins + idx]/(magReplica_squared[idx_T*mc_bins + idx]*magReplica_squared[idx_T*mc_bins + idx]);
			}
			else {
				bind_ratioReplica[idx_T*mc_bins + idx] = magReplica_four[idx_T*mc_bins + idx]*1e10;
			}
			/*
			if (print_all == true) {
				cout << "magReplica_squared[" << idx_T*mc_bins + idx << "] = " << magReplica_squared[idx_T*mc_bins + idx] << endl;
				cout << "magReplica_four[" << idx_T*mc_bins + idx << "] = " << magReplica_four[idx_T*mc_bins + idx] << endl;
				cout << "bind_ratioReplica[" << idx_T*mc_bins + idx << "] = " << bind_ratioReplica[idx_T*mc_bins + idx] << endl;
			}
			*/
		} // for idx

		for (unsigned long int jdx = 0; jdx < mc_bins; jdx++) {
			for (int idx_replica = 0; idx_replica <= 1; idx_replica++) {

				sigma_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += (mag1[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica])*(mag1[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
				sigma_magnetization2_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += (mag2[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetization2_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica])*(mag2[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetization2_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
				sigma_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += (magAF[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica])*(magAF[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);

				sigma_magnetization1_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += (mag1_squared[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetization1_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica])*(mag1_squared[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetization1_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
				sigma_magnetization2_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += (mag2_squared[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetization2_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica])*(mag2_squared[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetization2_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
				sigma_magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += (magAF_squared[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica])*(magAF_squared[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);

				sigma_magnetization1_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += (mag1_four[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetization1_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica])*(mag1_four[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetization1_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
				sigma_magnetizationAF_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += (magAF_four[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetizationAF_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica])*(magAF_four[idx_T*2*mc_bins + 2*jdx + idx_replica] - magnetizationAF_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);

				sigma_chi_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += (chi1[idx_T*2*mc_bins + 2*jdx + idx_replica] - chi_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica])*(chi1[idx_T*2*mc_bins + 2*jdx + idx_replica] - chi_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
				sigma_chi_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += (chiAF[idx_T*2*mc_bins + 2*jdx + idx_replica] - chi_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica])*(chiAF[idx_T*2*mc_bins + 2*jdx + idx_replica] - chi_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);

				sigma_binder_ratio_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += (bind_ratio[idx_T*2*mc_bins + 2*jdx + idx_replica] - binder_ratio_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica])*(bind_ratio[idx_T*2*mc_bins + 2*jdx + idx_replica] - binder_ratio_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
				sigma_binder_ratioAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += (bind_ratioAF[idx_T*2*mc_bins + 2*jdx + idx_replica] - binder_ratioAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica])*(bind_ratioAF[idx_T*2*mc_bins + 2*jdx + idx_replica] - binder_ratioAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);

				sigma_specific_heat_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] += (cV[idx_T*2*mc_bins + 2*jdx + idx_replica] - specific_heat_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica])*(cV[idx_T*2*mc_bins + 2*jdx + idx_replica] - specific_heat_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);

			} // for idx_replica

			sigma_xi_a_av_dis[idx_T*disorder_steps + idx_disorder] += (corr_a[idx_T*mc_bins + jdx] - xi_a_av_dis[idx_T*disorder_steps + idx_disorder])*(corr_a[idx_T*mc_bins + jdx] - xi_a_av_dis[idx_T*disorder_steps + idx_disorder]);
			sigma_xi_b_av_dis[idx_T*disorder_steps + idx_disorder] += (corr_b[idx_T*mc_bins + jdx] - xi_b_av_dis[idx_T*disorder_steps + idx_disorder])*(corr_b[idx_T*mc_bins + jdx] - xi_b_av_dis[idx_T*disorder_steps + idx_disorder]);

			sigma_magnetization_replica_av_dis[idx_T*disorder_steps + idx_disorder] += (magReplica[idx_T*mc_bins + jdx] - magnetization_replica_av_dis[idx_T*disorder_steps + idx_disorder])*(magReplica[idx_T*mc_bins + jdx] - magnetization_replica_av_dis[idx_T*disorder_steps + idx_disorder]);
			sigma_magnetization_replica_squared_av_dis[idx_T*disorder_steps + idx_disorder] += (magReplica_squared[idx_T*mc_bins + jdx] - magnetization_replica_squared_av_dis[idx_T*disorder_steps + idx_disorder])*(magReplica_squared[idx_T*mc_bins + jdx] - magnetization_replica_squared_av_dis[idx_T*disorder_steps + idx_disorder]);
			sigma_chi_magnetization_replica_av_dis[idx_T*disorder_steps + idx_disorder] += (chiReplica[idx_T*mc_bins + jdx] - chi_magnetization_replica_av_dis[idx_T*disorder_steps + idx_disorder])*(chiReplica[idx_T*mc_bins + jdx] - chi_magnetization_replica_av_dis[idx_T*disorder_steps + idx_disorder]);
			sigma_magnetization_replica_four_av_dis[idx_T*disorder_steps + idx_disorder] += (magReplica_four[idx_T*mc_bins + jdx] - magnetization_replica_four_av_dis[idx_T*disorder_steps + idx_disorder])*(magReplica_four[idx_T*mc_bins + jdx] - magnetization_replica_four_av_dis[idx_T*disorder_steps + idx_disorder]);
			sigma_binder_ratio_replica_av_dis[idx_T*disorder_steps + idx_disorder] += (bind_ratioReplica[idx_T*mc_bins + jdx] - binder_ratio_replica_av_dis[idx_T*disorder_steps + idx_disorder])*(bind_ratioReplica[idx_T*mc_bins + jdx] - binder_ratio_replica_av_dis[idx_T*disorder_steps + idx_disorder]);

		} // for jdx

		for (int idx_replica = 0; idx_replica <= 1; idx_replica++) {
			sigma_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = sqrt(sigma_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
			sigma_magnetization2_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]  = sqrt(sigma_magnetization2_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
			sigma_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = sqrt(sigma_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]); 

			sigma_magnetization1_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = sqrt(sigma_magnetization1_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]); 
			sigma_magnetization2_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = sqrt(sigma_magnetization2_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
			sigma_magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = sqrt(sigma_magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);

			sigma_magnetization1_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = sqrt(sigma_magnetization1_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
			sigma_magnetizationAF_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = sqrt(sigma_magnetizationAF_four_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);

			sigma_chi_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = sqrt(sigma_chi_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
			sigma_chi_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = sqrt(sigma_chi_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);

			sigma_binder_ratio_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = sqrt(sigma_binder_ratio_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);
			sigma_binder_ratioAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = sqrt(sigma_binder_ratioAF_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);

			sigma_specific_heat_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica] = sqrt(sigma_specific_heat_av_dis[idx_T*2*disorder_steps + 2*idx_disorder + idx_replica]);

		} // for idx_replica

		sigma_xi_a_av_dis[idx_T*disorder_steps + idx_disorder] = sqrt(sigma_xi_a_av_dis[idx_T*disorder_steps + idx_disorder]);
		sigma_xi_b_av_dis[idx_T*disorder_steps + idx_disorder] = sqrt(sigma_xi_b_av_dis[idx_T*disorder_steps + idx_disorder]);

		sigma_magnetization_replica_av_dis[idx_T*disorder_steps + idx_disorder] = sqrt(sigma_magnetization_replica_av_dis[idx_T*disorder_steps + idx_disorder]);
		sigma_magnetization_replica_squared_av_dis[idx_T*disorder_steps + idx_disorder] = sqrt(sigma_magnetization_replica_squared_av_dis[idx_T*disorder_steps + idx_disorder]);
		sigma_magnetization_replica_four_av_dis[idx_T*disorder_steps + idx_disorder] = sqrt(sigma_magnetization_replica_four_av_dis[idx_T*disorder_steps + idx_disorder]);
		sigma_binder_ratio_replica_av_dis[idx_T*disorder_steps + idx_disorder] = sqrt(sigma_binder_ratio_replica_av_dis[idx_T*disorder_steps + idx_disorder]);

	} // for idx_T	
	
	if (idx_disorder == 0) {
		cout << "Acceptance ratio Metropolis for T_min = " << static_cast<double>(acceptance_number_Metropolis)/static_cast<double>((thermalization_bins + mc_bins)*mc_binsize*N) << endl;
		cout << "Acceptance ratio PT = " << static_cast<double>(acceptance_number_pt)/static_cast<double>((thermalization_bins + mc_bins)*mc_binsize*T_steps) << endl;
		cout << endl;
		t = clock() - t;
		cout << "Time per MC step = " << static_cast<double>(t)/CLOCKS_PER_SEC/((thermalization_bins + mc_bins)*mc_binsize) << " sec = " << static_cast<double>(t)/CLOCKS_PER_SEC/60./((thermalization_bins + mc_bins)*mc_binsize) << " min = " << static_cast<double>(t)/CLOCKS_PER_SEC/60./60./((thermalization_bins + mc_bins)*mc_binsize) << " h" << endl;
		cout << "Time for simulation = " << static_cast<double>(t)/CLOCKS_PER_SEC/((thermalization_bins + mc_bins)*mc_binsize)*(thermalization_bins + mc_bins)*mc_binsize*disorder_steps << " sec = " << static_cast<double>(t)/CLOCKS_PER_SEC/((thermalization_bins + mc_bins)*mc_binsize)*(thermalization_bins + mc_bins)*mc_binsize*disorder_steps/60. << " min = " << static_cast<double>(t)/CLOCKS_PER_SEC/((thermalization_bins + mc_bins)*mc_binsize)*(thermalization_bins + mc_bins)*mc_binsize*disorder_steps/60./60. << " h = " << static_cast<double>(t)/CLOCKS_PER_SEC/((thermalization_bins + mc_bins)*mc_binsize)*(thermalization_bins + mc_bins)*mc_binsize*disorder_steps/60./60./24. << " days." << endl;
	}

	seed_mc += 2*T_steps; // re-initialize RNG for MC simulation for each disorder realization
	if (print_all == true) {
		cout << "idx_disorder = " << idx_disorder << " of disorder_steps = " << disorder_steps << endl;
	}

	// Output results
	if (output_observables == true) {
	if (idx_disorder % 2 == 0) {
		outputFilename = "Data-AveragedObservables-L_" + L_string + "-Lz_" + Lz_string + "-disorder_filling_" + disorder_filling_string + "-J33_" + J33_string + "-J34_" + J34_string + "-J44_" + J44_string + "-TMin_" + T_min_string + "-TMax_" + T_max_string + "-TSteps_" + T_steps_string + ".dat";
	}
	else {
		outputFilename = "Data-AveragedObservables-L_" + L_string + "-Lz_" + Lz_string + "-disorder_filling_" + disorder_filling_string + "-J33_" + J33_string + "-J34_" + J34_string + "-J44_" + J44_string + "-TMin_" + T_min_string + "-TMax_" + T_max_string + "-TSteps_" + T_steps_string + "-2.dat";

	}
	outputFile.open(outputFilename);

	outputFile << "# idx_disorder = " << idx_disorder << " of disorder_steps = " << disorder_steps << endl;
	outputFile << "# L = " << L << endl;
	outputFile << "# Lz = " << Lz << endl;
	outputFile << "# disorder_filling = " << disorder_filling << endl;
	outputFile << "# J33 = " << J33 << endl;
	outputFile << "# J34 = " << J34 << endl;
	outputFile << "# J44 = " << J44 << endl;
	outputFile << "# T_min = " << temperatures[0] << endl;
	outputFile << "# T_max = " << temperatures[T_steps-1] << endl;
	outputFile << "# T_steps = " << T_steps << endl;
	outputFile << "# disorder_steps = " << disorder_steps << endl;
	outputFile << "# mc_bins = " << mc_bins << endl;
	outputFile << "# mc_binsize = " << mc_binsize << endl;
	outputFile << "# thermalization_bins = " << thermalization_bins << endl;
	outputFile << "# seed_disorder = " << seed_disorder << endl;
	outputFile << "# seed_mc = " << seed_mc << endl ;
	outputFile << "# 1: L; 2: Lz; 3: disorder_filling; 4: temperature; 5: magnetization1; 6: sigma_magnetization1; 7: magnetization2; 8: sigma_magnetization2; 9: magnetization1_squared; 10: sigma_magnetization1_squared; 11: magnetization2_squared; 12: sigma_magnetization2_squared; 13: chi_magnetization1; 14: sigma_chi_magnetization1; 15: energy; 16: energy_squared; 17: specific_heat; 18: sigma_specific_heat; 19: xi_a; 20: sigma_xi_a; 21: xi_b; 22: sigma_xi_b; 23: magnetization AF; 24: sigma_magnetizationAF; 25: magnetizationAF_squared; 26: sigma_magnetizationAF_squared; 27: chi_magnetizationAF; 28: sigma_chi_magnetizationAF; 29: magnetization_four; 30: sigma_magnetization_four; 31: binder_ratio; 32: sigma_binder_ratio; 33: magnetizationAF_four; 34: sigma_magnetizationAF_four; 35: binder_ratioAF; 36: sigma_binder_ratioAF; 37:magnetization_replica; 38: sigma_magnetization_replica; 39:magnetization_replica_squared; 40: sigma_magnetization_replica_squared; 41:chi_magnetization_replica; 42: sigma_chi_magnetization_replica; 43:magnetization_replica_four; 44: sigma_magnetization_replica_four; 45:binder_ratio_replica; 46: sigma_binder_ratio_replica;  " << endl;

	outputFile.setf(ios::scientific);

	// Initialize all magnetization_av[idx_T], etc to zero. Average over all idx_disorder + 1) disorder realizations thus far. 
	for (int idx = 0; idx < T_steps; idx++) {
		magnetization1_av[idx] = 0.;
		magnetization2_av[idx] = 0.;
		magnetizationAF_av[idx] = 0.;
		magnetization_replica_av[idx] = 0.;

		sigma_magnetization1_av[idx] = 0.;
		sigma_magnetization2_av[idx] = 0.;
		sigma_magnetizationAF_av[idx] = 0.;
		sigma_magnetization_replica_av[idx] = 0.;

		magnetization1_squared_av[idx] = 0.;
		magnetization2_squared_av[idx] = 0.;
		magnetizationAF_squared_av[idx] = 0.;
		magnetization_replica_squared_av[idx] = 0.;
		magnetization1_four_av[idx] = 0.;
		magnetizationAF_four_av[idx] = 0.;
		magnetization_replica_four_av[idx] = 0.;

		sigma_magnetization1_squared_av[idx] = 0.;
		sigma_magnetization2_squared_av[idx] = 0.;
		sigma_magnetizationAF_squared_av[idx] = 0.;
		sigma_magnetization_replica_squared_av[idx] = 0.;
		sigma_magnetization1_four_av[idx] = 0.;
		sigma_magnetizationAF_four_av[idx] = 0.;
		sigma_magnetization_replica_four_av[idx] = 0.;

		chi_magnetization1_av[idx] = 0.; // magnetic susceptibility
		sigma_chi_magnetization1_av[idx] = 0.;
		chi_magnetizationAF_av[idx] = 0.; // magnetic susceptibility
		sigma_chi_magnetizationAF_av[idx] = 0.;
		chi_magnetization_replica_av[idx] = 0.; // magnetic susceptibility
		sigma_chi_magnetization_replica_av[idx] = 0.;

		binder_ratio_av[idx] = 0.; // Binder ratio for FM order
		sigma_binder_ratio_av[idx] = 0.;
		binder_ratioAF_av[idx] = 0.;// Binder ratio for AF order
		sigma_binder_ratioAF_av[idx] = 0.;
		binder_ratio_replica_av[idx] = 0.; // Binder ratio of replica order
		sigma_binder_ratio_replica_av[idx] = 0.;

		energy_av[idx] = 0.;
		energy_squared_av[idx] = 0.;

		specific_heat_av[idx] = 0.; // specific heat
		sigma_specific_heat_av[idx] = 0.;

		xi_a_av[idx] = 0.;
		xi_b_av[idx] = 0.;

		sigma_xi_a_av[idx] = 0.;
		sigma_xi_b_av[idx] = 0.;
	}
	if (output_s_q == true) {
		for (int idx = 0; idx < T_steps*L*Lz*(L/2+1); idx++) {
			s_q_av[idx] = 0.;
		}
	}
	// Average everything over disorder
	#pragma omp parallel for
	for (int idx_T = 0; idx_T < T_steps; idx_T++) {
		#pragma omp parallel for
		for (int idx = 0; idx <= idx_disorder; idx++){
			for (int idx_replica = 0; idx_replica <= 1; idx_replica++) {
				magnetization1_av[idx_T] += magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				magnetization2_av[idx_T] += magnetization2_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				magnetizationAF_av[idx_T] += magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				magnetization1_squared_av[idx_T] += magnetization1_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				magnetization2_squared_av[idx_T] += magnetization2_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				magnetizationAF_squared_av[idx_T] += magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				magnetization1_four_av[idx_T] += magnetization1_four_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				magnetizationAF_four_av[idx_T] += magnetizationAF_four_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				chi_magnetization1_av[idx_T] += chi_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				chi_magnetizationAF_av[idx_T] += chi_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				binder_ratio_av[idx_T] += binder_ratio_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				binder_ratioAF_av[idx_T] += binder_ratioAF_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;

				sigma_magnetization1_av[idx_T] += sigma_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				sigma_magnetization2_av[idx_T] += sigma_magnetization2_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				sigma_magnetizationAF_av[idx_T] += sigma_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				sigma_magnetization1_squared_av[idx_T] += sigma_magnetization1_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				sigma_magnetization2_squared_av[idx_T] += sigma_magnetization2_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				sigma_magnetizationAF_squared_av[idx_T] += sigma_magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				sigma_magnetization1_four_av[idx_T] += sigma_magnetization1_four_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				sigma_magnetizationAF_four_av[idx_T] += sigma_magnetizationAF_four_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				sigma_chi_magnetization1_av[idx_T] += sigma_chi_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				sigma_chi_magnetizationAF_av[idx_T] += sigma_chi_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				sigma_binder_ratio_av[idx_T] += sigma_binder_ratio_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				sigma_binder_ratioAF_av[idx_T] += sigma_binder_ratioAF_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
			
				energy_av[idx_T] += energy_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				energy_squared_av[idx_T] += energy_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				specific_heat_av[idx_T] += specific_heat_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
				sigma_specific_heat_av[idx_T] += sigma_specific_heat_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;

				if (output_s_q == true) {
					for (int jdx = 0; jdx < L*Lz*(L/2+1); jdx++) {
						s_q_av[idx_T*L*Lz*(L/2+1) + jdx] += s_q_av_dis[idx_T*disorder_steps*L*Lz*(L/2+1) + idx*L*Lz*(L/2+1) + jdx]/static_cast<double>(idx_disorder + 1)/2.;
					}
				}

			xi_a_av[idx_T] += xi_a_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
			xi_b_av[idx_T] += xi_b_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
			sigma_xi_a_av[idx_T] += sigma_xi_a_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
			sigma_xi_b_av[idx_T] += sigma_xi_b_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(idx_disorder + 1)/2.;
		} // for idx_replica
		
		magnetization_replica_av[idx_T] += magnetization_replica_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(idx_disorder + 1);
		magnetization_replica_squared_av[idx_T] += magnetization_replica_squared_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(idx_disorder + 1);
		magnetization_replica_four_av[idx_T] += magnetization_replica_four_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(idx_disorder + 1);
		chi_magnetization_replica_av[idx_T] += chi_magnetization_replica_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(idx_disorder + 1);
		binder_ratio_replica_av[idx_T] += binder_ratio_replica_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(idx_disorder + 1);
		sigma_magnetization_replica_av[idx_T] += sigma_magnetization_replica_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(idx_disorder + 1);
		sigma_magnetization_replica_squared_av[idx_T] += sigma_magnetization_replica_squared_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(idx_disorder + 1);
		sigma_chi_magnetization_replica_av[idx_T] += sigma_chi_magnetization_replica_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(idx_disorder + 1);
		sigma_magnetization_replica_four_av[idx_T] += sigma_magnetization_replica_four_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(idx_disorder + 1);
		sigma_binder_ratio_replica_av[idx_T] += sigma_binder_ratio_replica_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(idx_disorder + 1);

	} // for idx < disorder_steps


outputFile << L << " " << Lz << " " << disorder_filling << " " << temperatures[idx_T] << " " << magnetization1_av[idx_T] << " " << sigma_magnetization1_av[idx_T] << " " << magnetization2_av[idx_T] << " " << sigma_magnetization2_av[idx_T] << " " << magnetization1_squared_av[idx_T] << " " << sigma_magnetization1_squared_av[idx_T] << " " << magnetization2_squared_av[idx_T] << " " << sigma_magnetization2_squared_av[idx_T] << " " << chi_magnetization1_av[idx_T] << " " << sigma_chi_magnetization1_av[idx_T] << " " << energy_av[idx_T] << " " << energy_squared_av[idx_T] << " " << specific_heat_av[idx_T] << " " << sigma_specific_heat_av[idx_T] << " "  << xi_a_av[idx_T] << " " << sigma_xi_a_av[idx_T] << " " << xi_b_av[idx_T] << " " << sigma_xi_b_av[idx_T] << " " << magnetizationAF_av[idx_T] << " " << sigma_magnetizationAF_av[idx_T] << " " << magnetizationAF_squared_av[idx_T] << " " << sigma_magnetizationAF_squared_av[idx_T] << " " << chi_magnetizationAF_av[idx_T] << " " << sigma_chi_magnetizationAF_av[idx_T] << " " << magnetization1_four_av[idx_T] << " " << sigma_magnetization1_four_av[idx_T] << " " << binder_ratio_av[idx_T] << " " << sigma_binder_ratio_av[idx_T] << " " << magnetizationAF_four_av[idx_T] << " " << sigma_magnetizationAF_four_av[idx_T] << " " << binder_ratioAF_av[idx_T] << " " << sigma_binder_ratioAF_av[idx_T] << " " << magnetization_replica_av[idx_T] << " " << sigma_magnetization_replica_av[idx_T] << " " << magnetization_replica_squared_av[idx_T] << " " << sigma_magnetization_replica_squared_av[idx_T] << " " << chi_magnetization_replica_av[idx_T] << " " << sigma_chi_magnetization_replica_av[idx_T] << " " << magnetization_replica_four_av[idx_T] << " " << sigma_magnetization_replica_four_av[idx_T] << " " << binder_ratio_replica_av[idx_T] << " " << sigma_binder_ratio_replica_av[idx_T] << endl;

} // for idx_T < T_steps
 outputFile.close(); 
} // end if output_observables == true

if (output_s_q == true) {
if (idx_disorder % 2 == 0) {
	outputFilename = "Data-Sq-L_" + L_string + "-Lz_" + Lz_string + "-disorder_filling_" + disorder_filling_string + "-J33_" + J33_string + "-J34_" + J34_string + "-J44_" + J44_string + "-TMin_" + T_min_string + "-TMax_" + T_max_string + "-TSteps_" + T_steps_string + ".dat";
}
else 
{
	outputFilename = "Data-Sq-L_" + L_string + "-Lz_" + Lz_string + "-disorder_filling_" + disorder_filling_string + "-J33_" + J33_string + "-J34_" + J34_string + "-J44_" + J44_string + "-TMin_" + T_min_string + "-TMax_" + T_max_string + "-TSteps_" + T_steps_string + "-2.dat";

}
outputFile.open(outputFilename);
outputFile << "# L = " << L << endl;
outputFile << "# Lz = " << Lz << endl;
outputFile << "# disorder_filling = " << disorder_filling << endl;
outputFile << "# J33 = " << J33 << endl;
outputFile << "# J34 = " << J34 << endl;
outputFile << "# J44 = " << J44 << endl;
outputFile << "# T_min = " << temperatures[0] << endl;
outputFile << "# T_max = " << temperatures[T_steps-1] << endl;
outputFile << "# T_steps = " << T_steps << endl;
outputFile << "# disorder_steps = " << disorder_steps << endl;
outputFile << "# mc_bins = " << mc_bins << endl;
outputFile << "# mc_binsize = " << mc_binsize << endl;
outputFile << "# thermalization_bins = " << thermalization_bins << endl;
outputFile << "# seed_disorder = " << seed_disorder << endl;
outputFile << "# seed_mc = " << seed_mc << endl ;
outputFile << "# idx = xdx + ydx*L + zdx*L*Lz with 0 <= idx L*Lz*(L/2+1) covering Fourier momentum space: 0 <= k_x <= pi, -pi < k_y <= pi, -pi < kz <= pi. " << endl;
outputFile << "# 1: s_q(idx); " << endl << endl;
outputFile.setf(ios::scientific);

for (int idx_T = 0; idx_T < T_steps; idx_T++) {
	for (int idx = 0; idx < L*Lz*(L/2+1); idx++) {
		outputFile << s_q_av[idx_T*L*Lz*(L/2+1) + idx] << endl;	
	}
}
outputFile.close();
}

	
} // for idx_disorder

/***********************************************************************
****     ALL FOR LOOPS FINISHED. NOW JUST OUTPUT THE FINAL DATA. *******
***********************************************************************/

if (print_all == true) {
	cout << "idx_disorder loop finished" << endl;
}

if (output_observables == true) {
// Initialize all magnetization_av[idx_T], etc to zero. Average over all idx_disorder + 1) disorder realizations thus far. 
	for (int idx = 0; idx < T_steps; idx++) {
		magnetization1_av[idx] = 0.;
		magnetization2_av[idx] = 0.;
		magnetizationAF_av[idx] = 0.;
		magnetization_replica_av[idx] = 0.;

		sigma_magnetization1_av[idx] = 0.;
		sigma_magnetization2_av[idx] = 0.;
		sigma_magnetizationAF_av[idx] = 0.;
		sigma_magnetization_replica_av[idx] = 0.;

		magnetization1_squared_av[idx] = 0.;
		magnetization2_squared_av[idx] = 0.;
		magnetizationAF_squared_av[idx] = 0.;
		magnetization_replica_squared_av[idx] = 0.;
		magnetization1_four_av[idx] = 0.;
		magnetizationAF_four_av[idx] = 0.;
		magnetization_replica_four_av[idx] = 0.;

		sigma_magnetization1_squared_av[idx] = 0.;
		sigma_magnetization2_squared_av[idx] = 0.;
		sigma_magnetizationAF_squared_av[idx] = 0.;
		sigma_magnetization_replica_squared_av[idx] = 0.;
		sigma_magnetization1_four_av[idx] = 0.;
		sigma_magnetizationAF_four_av[idx] = 0.;
		sigma_magnetization_replica_four_av[idx] = 0.;

		chi_magnetization1_av[idx] = 0.; // magnetic susceptibility
		sigma_chi_magnetization1_av[idx] = 0.;
		chi_magnetizationAF_av[idx] = 0.; // magnetic susceptibility
		sigma_chi_magnetizationAF_av[idx] = 0.;
		chi_magnetization_replica_av[idx] = 0.; // magnetic susceptibility
		sigma_chi_magnetization_replica_av[idx] = 0.;

		binder_ratio_av[idx] = 0.; // Binder ratio for FM order
		sigma_binder_ratio_av[idx] = 0.;
		binder_ratioAF_av[idx] = 0.;// Binder ratio for AF order
		sigma_binder_ratioAF_av[idx] = 0.;
		binder_ratio_replica_av[idx] = 0.; // Binder ratio of replica order
		sigma_binder_ratio_replica_av[idx] = 0.;

		energy_av[idx] = 0.;
		energy_squared_av[idx] = 0.;

		specific_heat_av[idx] = 0.; // specific heat
		sigma_specific_heat_av[idx] = 0.;

		xi_a_av[idx] = 0.;
		xi_b_av[idx] = 0.;

		sigma_xi_a_av[idx] = 0.;
		sigma_xi_b_av[idx] = 0.;
	}
	if (output_s_q == true) {
		for (int idx = 0; idx < T_steps*L*Lz*(L/2+1); idx++) {
			s_q_av[idx] = 0.;
		}
	}

// Output results 
outputFilename = "Data-AveragedObservables-L_" + L_string + "-Lz_" + Lz_string + "-disorder_filling_" + disorder_filling_string + "-J33_" + J33_string + "-J34_" + J34_string + "-J44_" + J44_string + "-TMin_" + T_min_string + "-TMax_" + T_max_string + "-TSteps_" + T_steps_string + ".dat";
outputFile.open(outputFilename);

outputFile << "# L = " << L << endl;
outputFile << "# Lz = " << Lz << endl;
outputFile << "# disorder_filling = " << disorder_filling << endl;
outputFile << "# J33 = " << J33 << endl;
outputFile << "# J34 = " << J34 << endl;
outputFile << "# J44 = " << J44 << endl;
outputFile << "# T_min = " << temperatures[0] << endl;
outputFile << "# T_max = " << temperatures[T_steps-1] << endl;
outputFile << "# T_steps = " << T_steps << endl;
outputFile << "# disorder_steps = " << disorder_steps << endl;
outputFile << "# mc_bins = " << mc_bins << endl;
outputFile << "# mc_binsize = " << mc_binsize << endl;
outputFile << "# thermalization_bins = " << thermalization_bins << endl;
outputFile << "# seed_disorder = " << seed_disorder << endl;
outputFile << "# seed_mc = " << seed_mc << endl ;
outputFile << "# 1: L; 2: Lz; 3: disorder_filling; 4: temperature; 5: magnetization1; 6: sigma_magnetization1; 7: magnetization2; 8: sigma_magnetization2; 9: magnetization1_squared; 10: sigma_magnetization1_squared; 11: magnetization2_squared; 12: sigma_magnetization2_squared; 13: chi_magnetization1; 14: sigma_chi_magnetization1; 15: energy; 16: energy_squared; 17: specific_heat; 18: sigma_specific_heat; 19: xi_a; 20: sigma_xi_a; 21: xi_b; 22: sigma_xi_b; 23: magnetization AF; 24: sigma_magnetizationAF; 25: magnetizationAF_squared; 26: sigma_magnetizationAF_squared; 27: chi_magnetizationAF; 28: sigma_chi_magnetizationAF; 29: magnetization_four; 30: sigma_magnetization_four; 31: binder_ratio; 32: sigma_binder_ratio; 33: magnetizationAF_four; 34: sigma_magnetizationAF_four; 35: binder_ratioAF; 36: sigma_binder_ratioAF; 37:magnetization_replica; 38: sigma_magnetization_replica; 39:magnetization_replica_squared; 40: sigma_magnetization_replica_squared; 41:chi_magnetization_replica; 42: sigma_chi_magnetization_replica; 43:magnetization_replica_four; 44: sigma_magnetization_replica_four; 45:binder_ratio_replica; 46: sigma_binder_ratio_replica;  " << endl;

outputFile.setf(ios::scientific);

// Average everything over disorder
for (int idx_T = 0; idx_T < T_steps; idx_T++) {
	for (int idx = 0; idx < disorder_steps; idx++) {
		for (int idx_replica = 0; idx_replica <= 1; idx_replica++) {
			magnetization1_av[idx_T] += magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			magnetization2_av[idx_T] += magnetization2_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			magnetizationAF_av[idx_T] += magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			magnetization1_squared_av[idx_T] += magnetization1_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			magnetization2_squared_av[idx_T] += magnetization2_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			magnetizationAF_squared_av[idx_T] += magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			magnetization1_four_av[idx_T] += magnetization1_four_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			magnetizationAF_four_av[idx_T] += magnetizationAF_four_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			chi_magnetization1_av[idx_T] += chi_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			chi_magnetizationAF_av[idx_T] += chi_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			binder_ratio_av[idx_T] += binder_ratio_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			binder_ratioAF_av[idx_T] += binder_ratioAF_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			
			sigma_magnetization1_av[idx_T] += sigma_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_magnetization2_av[idx_T] += sigma_magnetization2_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_magnetizationAF_av[idx_T] += sigma_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_magnetization1_squared_av[idx_T] += sigma_magnetization1_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_magnetization2_squared_av[idx_T] += sigma_magnetization2_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_magnetizationAF_squared_av[idx_T] += sigma_magnetizationAF_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_magnetization1_four_av[idx_T] += sigma_magnetization1_four_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_magnetizationAF_four_av[idx_T] += sigma_magnetizationAF_four_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_chi_magnetization1_av[idx_T] += sigma_chi_magnetization1_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_chi_magnetizationAF_av[idx_T] += sigma_chi_magnetizationAF_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_binder_ratio_av[idx_T] += sigma_binder_ratio_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_binder_ratioAF_av[idx_T] += sigma_binder_ratioAF_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			
			energy_av[idx_T] += energy_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			energy_squared_av[idx_T] += energy_squared_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			specific_heat_av[idx_T] += specific_heat_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_specific_heat_av[idx_T] += sigma_specific_heat_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;


			if (output_s_q == true) {
					for (int jdx = 0; jdx < L*Lz*(L/2+1); jdx++) {
						s_q_av[idx_T*L*Lz*(L/2+1) + jdx] += s_q_av_dis[idx_T*disorder_steps*L*Lz*(L/2+1) + idx*L*Lz*(L/2+1) + jdx]/static_cast<double>(disorder_steps)/2.;
					}
				}
				/*
			for (int jdx = 0; jdx < L*Lz*(L/2+1); jdx++) {
				s_q_av[idx_T*L*Lz*(L/2+1) + jdx] += s_q_av_dis[idx_T*disorder_steps*L*Lz*(L/2+1) + idx*L*Lz*(L/2+1) + jdx]/static_cast<double>(disorder_steps)/2.;
			}
			*/

			xi_a_av[idx_T] += xi_a_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			xi_b_av[idx_T] += xi_b_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_xi_a_av[idx_T] += sigma_xi_a_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
			sigma_xi_b_av[idx_T] += sigma_xi_b_av_dis[idx_T*2*disorder_steps + 2*idx + idx_replica]/static_cast<double>(disorder_steps)/2.;
		} // for idx_replica
		
		magnetization_replica_av[idx_T] += magnetization_replica_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(disorder_steps);
		magnetization_replica_squared_av[idx_T] += magnetization_replica_squared_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(disorder_steps);
		magnetization_replica_four_av[idx_T] += magnetization_replica_four_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(disorder_steps);
		chi_magnetization_replica_av[idx_T] += chi_magnetization_replica_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(disorder_steps);
		binder_ratio_replica_av[idx_T] += binder_ratio_replica_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(disorder_steps);
		sigma_magnetization_replica_av[idx_T] += sigma_magnetization_replica_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(disorder_steps);
		sigma_magnetization_replica_squared_av[idx_T] += sigma_magnetization_replica_squared_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(disorder_steps);
		sigma_chi_magnetization_replica_av[idx_T] += sigma_chi_magnetization_replica_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(disorder_steps);
		sigma_magnetization_replica_four_av[idx_T] += sigma_magnetization_replica_four_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(disorder_steps);
		sigma_binder_ratio_replica_av[idx_T] += sigma_binder_ratio_replica_av_dis[idx_T*disorder_steps + idx]/static_cast<double>(disorder_steps);
					
	} // for idx < disorder_steps

outputFile << L << " " << Lz << " " << disorder_filling << " " << temperatures[idx_T] << " " << magnetization1_av[idx_T] << " " << sigma_magnetization1_av[idx_T] << " " << magnetization2_av[idx_T] << " " << sigma_magnetization2_av[idx_T] << " " << magnetization1_squared_av[idx_T] << " " << sigma_magnetization1_squared_av[idx_T] << " " << magnetization2_squared_av[idx_T] << " " << sigma_magnetization2_squared_av[idx_T] << " " << chi_magnetization1_av[idx_T] << " " << sigma_chi_magnetization1_av[idx_T] << " " << energy_av[idx_T] << " " << energy_squared_av[idx_T] << " " << specific_heat_av[idx_T] << " " << sigma_specific_heat_av[idx_T] << " "  << xi_a_av[idx_T] << " " << sigma_xi_a_av[idx_T] << " " << xi_b_av[idx_T] << " " << sigma_xi_b_av[idx_T] << " " << magnetizationAF_av[idx_T] << " " << sigma_magnetizationAF_av[idx_T] << " " << magnetizationAF_squared_av[idx_T] << " " << sigma_magnetizationAF_squared_av[idx_T] << " " << chi_magnetizationAF_av[idx_T] << " " << sigma_chi_magnetizationAF_av[idx_T] << " " << magnetization1_four_av[idx_T] << " " << sigma_magnetization1_four_av[idx_T] << " " << binder_ratio_av[idx_T] << " " << sigma_binder_ratio_av[idx_T] << " " << magnetizationAF_four_av[idx_T] << " " << sigma_magnetizationAF_four_av[idx_T] << " " << binder_ratioAF_av[idx_T] << " " << sigma_binder_ratioAF_av[idx_T] << " " << magnetization_replica_av[idx_T] << " " << sigma_magnetization_replica_av[idx_T] << " " << magnetization_replica_squared_av[idx_T] << " " << sigma_magnetization_replica_squared_av[idx_T] << " " << chi_magnetization_replica_av[idx_T] << " " << sigma_chi_magnetization_replica_av[idx_T] << " " << magnetization_replica_four_av[idx_T] << " " << sigma_magnetization_replica_four_av[idx_T] << " " << binder_ratio_replica_av[idx_T] << " " << sigma_binder_ratio_replica_av[idx_T] << endl;

} // for idx_T < T_steps
 outputFile.close(); 
}

if (output_s_q == true) {
outputFilename = "Data-Sq-L_" + L_string + "-Lz_" + Lz_string + "-disorder_filling_" + disorder_filling_string + "-J33_" + J33_string + "-J34_" + J34_string + "-J44_" + J44_string + "-TMin_" + T_min_string + "-TMax_" + T_max_string + "-TSteps_" + T_steps_string + ".dat";
outputFile.open(outputFilename);
outputFile << "# L = " << L << endl;
outputFile << "# Lz = " << Lz << endl;
outputFile << "# disorder_filling = " << disorder_filling << endl;
outputFile << "# J33 = " << J33 << endl;
outputFile << "# J34 = " << J34 << endl;
outputFile << "# J44 = " << J44 << endl;
outputFile << "# T_min = " << temperatures[0] << endl;
outputFile << "# T_max = " << temperatures[T_steps-1] << endl;
outputFile << "# T_steps = " << T_steps << endl;
outputFile << "# disorder_steps = " << disorder_steps << endl;
outputFile << "# mc_bins = " << mc_bins << endl;
outputFile << "# mc_binsize = " << mc_binsize << endl;
outputFile << "# thermalization_bins = " << thermalization_bins << endl;
outputFile << "# seed_disorder = " << seed_disorder << endl;
outputFile << "# seed_mc = " << seed_mc << endl ;
outputFile << "# idx = xdx + ydx*L + zdx*L*Lz with 0 <= idx L*Lz*(L/2+1) covering Fourier momentum space: 0 <= k_x <= pi, -pi < k_y <= pi, -pi < kz <= pi. " << endl;
outputFile << "# 1: s_q(idx); " << endl << endl;
outputFile.setf(ios::scientific);

for (int idx_T = 0; idx_T < T_steps; idx_T++) {
	for (int idx = 0; idx < L*Lz*(L/2+1); idx++) {
		outputFile << s_q_av[idx_T*L*Lz*(L/2+1) + idx] << endl;	
	}
}
outputFile.close();
}

double fraction_filled_av = 0.;
// disorder averaged probability of bond interactions {J33, J34, J44}
double fraction_bonds_J33_av = 0.;
double fraction_bonds_J34_av = 0.;
double fraction_bonds_J44_av = 0.;

double sigma_fraction_filled = 0.;
double sigma_fraction_bonds_J33 = 0.;
double sigma_fraction_bonds_J34 = 0.;
double sigma_fraction_bonds_J44 = 0.;

// only use one of the two replicas for now
for (int idx = 0; idx < disorder_steps; idx++) {
	fraction_filled_av += fraction_filled[idx]/static_cast<double>(disorder_steps);
	fraction_bonds_J33_av += fraction_bonds_J33[idx]/static_cast<double>(disorder_steps);
	fraction_bonds_J34_av += fraction_bonds_J34[idx]/static_cast<double>(disorder_steps);
	fraction_bonds_J44_av += fraction_bonds_J44[idx]/static_cast<double>(disorder_steps);
}

// Perform statistical analysis of filled sites and bond interactions {J33, J34, J44}
for (int idx = 0; idx < disorder_steps; idx++) {
	sigma_fraction_filled += (fraction_filled[idx] - fraction_filled_av)*(fraction_filled[idx] - fraction_filled_av);
	sigma_fraction_bonds_J33 += (fraction_bonds_J33[idx] - fraction_bonds_J33_av)*(fraction_bonds_J33[idx] - fraction_bonds_J33_av);
	sigma_fraction_bonds_J34 += (fraction_bonds_J34[idx] - fraction_bonds_J34_av)*(fraction_bonds_J34[idx] - fraction_bonds_J34_av);
	sigma_fraction_bonds_J44 += (fraction_bonds_J44[idx] - fraction_bonds_J44_av)*(fraction_bonds_J44[idx] - fraction_bonds_J44_av);
}
sigma_fraction_filled /= disorder_steps*(disorder_steps - 1);
sigma_fraction_bonds_J33 /= disorder_steps*(disorder_steps - 1);
sigma_fraction_bonds_J34 /= disorder_steps*(disorder_steps - 1);
sigma_fraction_bonds_J44 /= disorder_steps*(disorder_steps - 1);

sigma_fraction_filled = sqrt(sigma_fraction_filled);
sigma_fraction_bonds_J33 = sqrt(sigma_fraction_bonds_J33);
sigma_fraction_bonds_J34 = sqrt(sigma_fraction_bonds_J34);
sigma_fraction_bonds_J44 = sqrt(sigma_fraction_bonds_J44);

double fraction_bonds_all = fraction_bonds_J33_av + fraction_bonds_J34_av + fraction_bonds_J44_av;
double sigma_fraction_bonds_all = sigma_fraction_bonds_J33 + sigma_fraction_bonds_J34 + sigma_fraction_bonds_J44;

if (output_data_fraction == true) {
outputFilename = "Data-Fraction.dat";
outputFile.open(outputFilename, ios::out | ios::app);
//outputFile << "# 1: disorder_filling; 2: fraction_filled_av; 3: sigma_fraction_filled_av; 4: fraction_bonds_J33_av; 5: sigma_J33; 6: fraction_bonds_J34_av; 7: sigma_J34; 8: fraction_bonds_J44_av; 9: sigma_J44; 10 fraction_bonds_all; 11: sigma_all "<< endl;
outputFile << disorder_filling << " " << fraction_filled_av << " " << sigma_fraction_filled << " " << fraction_bonds_J33_av << " " << sigma_fraction_bonds_J33 << " " << fraction_bonds_J34_av << " " << sigma_fraction_bonds_J34 << " " << fraction_bonds_J44_av << " " << sigma_fraction_bonds_J44 <<  " " << fraction_bonds_all << " " << sigma_fraction_bonds_all << endl; 

outputFile.close();
}
/*cout << "Replica 1 at low T = " << temper	atures[0] << ": " << endl;
	for (int idx = 0; idx < N; idx++) {
		if (spins[idx] == 1) {cout << " +1 ";}
		else if  (spins[idx] == 2) {cout << " +2 ";}
		else if  (spins[idx] == -1) {cout << " -1 ";}
		else if  (spins[idx] == -2) {cout << " -2 ";}
		else {cout << "  0 ";}
		if ((idx + 1) % L == 0) {cout << endl;}
	}*/


	t_total = clock() - t_total; // time difference between current clock and clock at beginning of main
	cout << "Disorder_filling = " << disorder_filling << ": Actual total runtime = " << static_cast<double>(t_total)/CLOCKS_PER_SEC << " sec = " << static_cast<double>(t_total)/CLOCKS_PER_SEC/60. << " min = " << static_cast<double>(t_total)/CLOCKS_PER_SEC/60./60. << " h" << endl;

	cout << endl;
	return 0;

} // end main