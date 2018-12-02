#include "include.h"

std::mt19937* ef_generator;
std::mt19937* genome_generator;
std::mt19937* patch_generator;
std::mt19937* main_generator;

void initialize_rand_generators(int seed){
    ef_generator = new std::mt19937(seed);
    genome_generator = new std::mt19937(seed);
    patch_generator = new std::mt19937(seed);
    main_generator = new std::mt19937(seed);
}

double beta_dist(double alpha, double beta, std::mt19937* gen){
    std::gamma_distribution<double> dis1(alpha, 1);
    std::gamma_distribution<double> dis2(beta, 1);

    double x1 = dis1(*gen);
    double x2 = dis2(*gen);

    double res = double(x1)/double(x1+x2);
    return res;
}


double std_normal(std::mt19937 *gen){
    std::normal_distribution<double> std_norm_dis(0,1);
    return std_norm_dis(*gen);
}

double normal(double mu, double sigma, std::mt19937* gen){
    std::normal_distribution<double> norm_dis(mu,sigma);
    return norm_dis(*gen);
}


double poisson(double lambda, std::mt19937* gen){
    std::poisson_distribution<int> p_dis(lambda);
    return p_dis(*gen);
}

int int_uniform(int a, int b, std::mt19937* gen){
    std::uniform_int_distribution<int> dis(a,b);
    return dis(*gen);
}


double real_uniform(double a, double b, std::mt19937* gen){
    std::uniform_real_distribution<double> dis(a,b);
    return dis(*gen);
}
