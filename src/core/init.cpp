#include "include.h"
#include "Patch.h"
#include "EnvFactor.h"
#include "GenomeDict.h"
#include "Individual.h"

void init(){
    read_params_file();
    initialize_rand_generators((int) params["RANDOM_SEED"]);
    initialize_env_factors();
    initialize_patches();
    init_dist_matrix();
    initialize_genome_dict();
    initialize_individuals();
    initialize_genomes();
}

void read_params_file(){
    std::ifstream infile("params.ini");
    std::string line;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string name;
        double val;

        iss >> name;
        iss >> val;

        params.insert(std::pair<std::string, double>(name, val));
    }
}

void initialize_env_factors(){
    int n_ef = params["NUM_ENV_FACTORS"];
    int size = params["ENV_FACTOR_RESOLUTION"];
    int h_val = params["ENV_FACTOR_H_VALUE"];

    envFactors = new std::vector<EnvFactor*>;

    for (int i = 0; i < n_ef; i++){
        EnvFactor* ef_i = new EnvFactor(size, h_val);
        envFactors->push_back(ef_i);
    }
}

void initialize_patches(){
    int n_patches = int(params["NUM_PATCHES"]);
    double k_mean =  params["PATCH_K_MEAN"];
    double k_sigma = params["PATCH_K_SD"];
    double side_len = params["SIDE_LENGTH"];

    std::normal_distribution<double> k_dis(k_mean, k_sigma);
    std::uniform_real_distribution<double> x_dis(0, side_len);
    std::uniform_real_distribution<double> y_dis(0, side_len);

    patches = new std::vector<Patch*>;

    for (int i = 0; i < n_patches; i++){
        Patch *tmp = new Patch(x_dis(*patch_generator), y_dis(*patch_generator), k_dis(*patch_generator));
        patches->push_back(tmp);
    }
}

void init_dist_matrix(){
    int n_patches = params["NUM_PATCHES"];
    Patch* patch_i;
    Patch* patch_j;
    for (int i = 0; i < n_patches; i++){
        std::vector<double> row;
        for (int j = 0; j < n_patches; j++){
            patch_i = (*patches)[i];
            patch_j = (*patches)[j];
            double x_dist = patch_i->get_x() - patch_j->get_x();
            double y_dist = patch_i->get_y() - patch_j->get_y();
            double d2 = pow((x_dist),2)+pow((y_dist),2);
            row.push_back(sqrt(d2));
        }
        dist_matrix.push_back(row);
    }
}

void initialize_genome_dict(){
    genome_dict = new GenomeDict();
}

void initialize_individuals(){
    double k;
    int n_indiv, i;
    Individual* indiv_i;

    for (Patch* patch_i : *patches){
        k = patch_i->get_K();
        std::uniform_int_distribution<int> n_indiv_dis(0, int(k));
        n_indiv = n_indiv_dis(*patch_generator);

        for (i = 0; i < n_indiv; i++){
            indiv_i = new Individual(patch_i);
            patch_i->add_individual(indiv_i);
        }
    }
}


void initialize_genomes(){
    std::vector<std::vector<double>> alleles = gen_init_alleles();
    std::vector<std::vector<double>> freqs = gen_init_allele_freqs(alleles);
    int n_loci = params["NUM_OF_LOCI"];

    for (Patch* patch_i : *patches){
        for (Individual* indiv_i : patch_i->get_all_individuals()){
            for (int locus = 0; locus < n_loci; locus++){
                double al1 = get_allele_from_dist(alleles, freqs, locus);
                double al2 = get_allele_from_dist(alleles, freqs, locus);
                indiv_i->set_locus(locus, 0, al1);
                indiv_i->set_locus(locus, 1, al2);
            }
        }
    }

}

double get_allele_from_dist(std::vector<std::vector<double>>  alleles, std::vector<std::vector<double>> freqs, int locus){
    double draw, sum;
    double allele_val = -1.0;
    int n_alleles;

    std::uniform_real_distribution<double> dis(0.0,1.0);

    n_alleles = freqs[locus].size();
    sum = 0.0;
    draw = dis(*patch_generator);

    for (int allele_i = 0; allele_i < n_alleles; allele_i++){
        sum += freqs[locus][allele_i];
        if (sum > draw){
            allele_val = alleles[locus][allele_i];
            break;
        }
    }
    return allele_val;
}

std::vector<std::vector<double>> gen_init_alleles(){
    int n_loci = params["NUM_OF_LOCI"];

    std::uniform_real_distribution<double> allele_dis(0.0, 1.0);
    int n_indiv = get_total_population_size();

    std::vector<std::vector<double>> alleles;
    int n_alleles, allele_i;
    double allele;

    for (int locus = 0; locus < n_loci; locus++){
        n_alleles = expected_num_alleles(n_indiv);
        alleles.push_back(std::vector<double>());
        for (allele_i = 0; allele_i < n_alleles; allele_i++){
            allele = allele_dis(*genome_generator);
            alleles[locus].push_back(allele);
        }
    }
    return alleles;
}

int expected_num_alleles(int n_indiv){
    std::uniform_int_distribution<int> num_allele_dis(20, 40);
    return num_allele_dis(*genome_generator);
}

std::vector<std::vector<double>> gen_init_allele_freqs(std::vector<std::vector<double>> alleles){
    int n_loci = params["NUM_OF_LOCI"];
    std::vector<std::vector<double>> freqs;


    for (int locus = 0; locus < n_loci; locus++){
        double rem = 1.0;
        double p, prop;
        std::vector<double> props;

        for (int al = 0; al < alleles[locus].size(); al++){
                p = beta_dist(0.6, 1.7, genome_generator);
                prop = rem * p;
                rem -= prop;
                props.push_back(prop);
        }
        props.push_back(rem);
        freqs.push_back(props);
    }

    return freqs;
}

int get_total_population_size(){
    int ct = 0;
    for (Patch* patch_i : *patches){
        ct += patch_i->get_size();
    }
    return ct;
}
