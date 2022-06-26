// Define global parameters of the program
bool output_observables;
bool output_data_fraction;
bool output_s_q; 
bool print_all;
bool snapshot_output;
bool parallel_tempering_step; 

std::ofstream outputFile; // output file stream
std::ifstream inputFile; // output file stream

// Parameters read in when program is started
int dim;
int L, Lz;
double J33, J34, J44;
double T_min, T_max;
int T_steps;
double disorder_filling;
int disorder_steps;
unsigned long int mc_bins, mc_binsize, thermalization_bins, mc_steps;
unsigned long int seed_disorder, seed_mc;

// Temperature table
double *temperatures;

// Number of sites
int N2D, N3D, N;

// Number of neighbors (different in 2D and 3D)
int num_neighbors_direct; // number of nearest neighbor sites
int num_neighbors_face; // number of face diagonal neighbor sites
int num_neighbors_body; // number of body diagonal neighbor sites
int num_neighbors_direct_nn; // number of direct next-nearest neighbor sites
// Table of neighbor sites
int *neighbors_direct;
int *neighbors_face;
int *neighbors_body;
int *neighbors_direct_nn;

// Table of spin flip probabilities
double *prob_1, *prob_2;

// Array of spins 
int *spins;
int *order;
double *spins_energy;
bool swap_even_pairs;

// Variables used in MC runs for a fixed disorder realization
double *magnetization1; 
double *magnetization2;
double *magnetizationAF;
double *magnetization_replica;
double *magnetization1_squared; 
double *magnetization2_squared;
double *magnetizationAF_squared;
double *magnetization_replica_squared;
double *magnetization1_four; 
double *magnetizationAF_four;
double *magnetization_replica_four;
double *energy;
double *energy_squared;
double *specific_heat;
double *s_q;
double *xi_a;
double *xi_b;
unsigned long int acceptance_number_Metropolis;
unsigned long int acceptance_number_pt;

// Variables used as MC averages for FIXED disorder realization
double *magnetization1_av_dis;
double *magnetization2_av_dis;
double *magnetizationAF_av_dis;
double *magnetization_replica_av_dis;

double *magnetization1_squared_av_dis;
double *magnetization2_squared_av_dis;
double *magnetizationAF_squared_av_dis;
double *magnetization_replica_squared_av_dis;

double *magnetization1_four_av_dis;
double *magnetizationAF_four_av_dis;
double *magnetization_replica_four_av_dis;

double *sigma_magnetization1_av_dis;
double *sigma_magnetization2_av_dis;
double *sigma_magnetizationAF_av_dis;
double *sigma_magnetization_replica_av_dis;

double *sigma_magnetization1_squared_av_dis;
double *sigma_magnetization2_squared_av_dis;
double *sigma_magnetizationAF_squared_av_dis;
double *sigma_magnetization_replica_squared_av_dis;

double *sigma_magnetization1_four_av_dis;
double *sigma_magnetizationAF_four_av_dis;
double *sigma_magnetization_replica_four_av_dis;

double *chi_magnetization1_av_dis;
double *sigma_chi_magnetization1_av_dis;
double *chi_magnetizationAF_av_dis;
double *sigma_chi_magnetizationAF_av_dis;
double *chi_magnetization_replica_av_dis;
double *sigma_chi_magnetization_replica_av_dis;

double *binder_ratio_av_dis;
double *sigma_binder_ratio_av_dis;
double *binder_ratioAF_av_dis;
double *sigma_binder_ratioAF_av_dis;
double *binder_ratio_replica_av_dis;
double *sigma_binder_ratio_replica_av_dis;

double *energy_av_dis;;
double *energy_squared_av_dis;
double *specific_heat_av_dis;
double *sigma_specific_heat_av_dis;
double *s_q_av_dis;
double *xi_a_av_dis;
double *xi_b_av_dis;
double *sigma_xi_a_av_dis;
double *sigma_xi_b_av_dis;

// Output variables averaged over disorder. Array in T_steps.
double *magnetization1_av;
double *magnetization2_av;
double *magnetizationAF_av;
double *magnetization_replica_av;

double *sigma_magnetization1_av;
double *sigma_magnetization2_av;
double *sigma_magnetizationAF_av;
double *sigma_magnetization_replica_av;

double *magnetization1_squared_av;
double *magnetization2_squared_av;
double *magnetizationAF_squared_av;
double *magnetization_replica_squared_av;
double *magnetization1_four_av;
double *magnetizationAF_four_av;
double *magnetization_replica_four_av;

double *sigma_magnetization1_squared_av;
double *sigma_magnetization2_squared_av;
double *sigma_magnetizationAF_squared_av;
double *sigma_magnetization_replica_squared_av;
double *sigma_magnetization1_four_av;
double *sigma_magnetizationAF_four_av;
double *sigma_magnetization_replica_four_av;

double *chi_magnetization1_av; // magnetic susceptibility
double *sigma_chi_magnetization1_av;
double *chi_magnetizationAF_av; // magnetic susceptibility
double *sigma_chi_magnetizationAF_av;
double *chi_magnetization_replica_av; // replica susceptibility
double *sigma_chi_magnetization_replica_av;
double *binder_ratio_av; // Binder ratio for FM order
double *sigma_binder_ratio_av;
double *binder_ratioAF_av; // Binder ratio for AF order
double *sigma_binder_ratioAF_av;
double *binder_ratio_replica_av; // Binder ratio of replica order
double *sigma_binder_ratio_replica_av;

double *energy_av;
double *energy_squared_av;

double *specific_heat_av; // specific heat
double *sigma_specific_heat_av;


double *s_q_av;

double *xi_a_av;
double *xi_b_av;
double *sigma_xi_a_av;
double *sigma_xi_b_av;

double *fraction_filled; // measure the fraction of filled sites
double *fraction_bonds_J33; // measure the fraction of J33 bonds
double *fraction_bonds_J34; // J34 bonds
double *fraction_bonds_J44; // J44 bonds

