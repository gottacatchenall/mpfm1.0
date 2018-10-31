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

void init();
void migration();
void selection();
void mating();
void logging();
void census();
void update_progress_bar(int gen);

// init
void read_params_file();
void initialize_rand_generators(int random_seed);
void initialize_env_factors();
void initialize_patches();
void init_dist_matrix();
void initialize_genome_dict();
void initialize_individuals();
void initialize_genomes();
std::vector<std::vector<double>> gen_init_alleles();
std::vector<std::vector<double>> gen_init_allele_freqs(std::vector<std::vector<double>> alleles);
double get_allele_from_dist(std::vector<std::vector<double>>  alleles, std::vector<std::vector<double>> freqs, int locus);
int expected_num_alleles(int n_indiv);
int get_total_population_size();

// Random Gen
double beta_dist(double alpha, double beta, std::mt19937* gen);
double std_normal(std::mt19937* gen);
double poisson(int lambda, std::mt19937* gen);
int int_uniform(int a, int b, std::mt19937* gen);
double real_uniform(double a, double b, std::mt19937* gen);

extern std::map<std::string, double> params;
extern std::vector<EnvFactor*>* envFactors;
extern std::vector<Patch*>* patches;
extern std::vector<std::vector<double>> dist_matrix;
extern std::mt19937* ef_generator;
extern std::mt19937* genome_generator;
extern std::mt19937* patch_generator;
extern std::mt19937* main_generator;
extern GenomeDict* genome_dict;
