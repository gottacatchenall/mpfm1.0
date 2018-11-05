#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <assert.h>
#include <cmath>
#include <map>
#include <set>
#include <unordered_map>

class Individual;
class EnvFactor;
class Individual;
class GenomeDict;
class Patch;
class AlleleTracker;
class MigrationTracker;

void init();
void migration();
void selection();
void mating();
void logging();
void census();
void update_progress_bar(int gen);
void check_dispersal();
// init
void read_params_file();
void initialize_rand_generators(int random_seed);
void initialize_env_factors();
void initialize_patches();
void init_dist_matrix();
void initialize_genome_dict();
void initialize_individuals();
void initialize_genomes();

void setup_initial_alleles();
std::vector<std::vector<double>> gen_alleles();
std::vector<std::vector<double>> generate_allele_freq_from_beta(std::vector<std::vector<double>> alleles);
std::vector<std::vector<int>>  get_num_with_each_allele(std::vector<std::vector<double>> props, int patch_size);


void initialize_gene_tracker();
std::vector<std::vector<double>> gen_init_alleles();
std::vector<std::vector<double>> gen_init_allele_freqs(std::vector<std::vector<double>> alleles);
double get_allele_from_dist(std::vector<std::vector<double>>  alleles, std::vector<std::vector<double>> freqs, int locus);
int expected_num_alleles();
int get_total_population_size();

// Logging
void log_patch(Patch* patch_i);
void log_fst(std::vector<double> f_st);
void log_eff_migration(Patch* patch_i);
void log_attempted_migration(Patch* patch_i, Patch* patch_j);
void log_successful_migration(Patch* patch_i, Patch* patch_j);
void log_population(Patch* patch_i);
void log_linkage(int patch_num, int l1, double al1,  int l2, double al2, double D);
void log_allele_freq(int patch_num, int locus, double allele_val, double freq);
void log_env_factors(Patch* patch_i);
double mean_w(Patch* patch_i);
double sd_w(Patch* patch_i, double w_bar);
bool is_file_empty(std::string path);

std::vector<double> get_fst();


// Random Gen
double beta_dist(double alpha, double beta, std::mt19937* gen);
double std_normal(std::mt19937* gen);
double normal(double mu, double sigma, std::mt19937* gen);
double poisson(int lambda, std::mt19937* gen);
int int_uniform(int a, int b, std::mt19937* gen);
double real_uniform(double a, double b, std::mt19937* gen);

extern int generation;
extern std::map<std::string, double> params;
extern std::vector<EnvFactor*>* envFactors;
extern std::vector<Patch*>* patches;
extern std::vector<std::vector<double>> dist_matrix;
extern std::mt19937* ef_generator;
extern std::mt19937* genome_generator;
extern std::mt19937* patch_generator;
extern std::mt19937* main_generator;
extern MigrationTracker* migration_tracker;
extern GenomeDict* genome_dict;
