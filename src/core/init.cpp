#include "include.h"
#include "Patch.h"
#include "EnvFactor.h"
#include "GenomeDict.h"
#include "Individual.h"
#include "AlleleTracker.h"
#include "MigrationTracker.h"

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
        std::vector<std::string> record;
        while (iss){
          std::string s;
          if (!getline(iss, s, ',' )) break;
          record.push_back( s );
        }

        std::string name = record[0];
        double val = atof(record[1].c_str());

        std::cout << "name: " << name << " val: " << val << "\n";

        params.insert(std::pair<std::string, double>(name, val));
    }
}

void initialize_env_factors(){
    int n_ef = params["NUM_ENV_FACTORS"];
    int size = params["SIDE_LENGTH"];
    double h_val = params["ENV_FACTOR_H_VALUE"];

    envFactors = new std::vector<EnvFactor*>;

    for (int i = 0; i < n_ef; i++){
        EnvFactor* ef_i = new EnvFactor(size, h_val);
        envFactors->push_back(ef_i);
    }

    write_efs();
}

void initialize_patches(){
    int n_patches = params["NUM_PATCHES"];
    int n_indiv = params["N_INDIVIDUALS"];
    double side_len = params["SIDE_LENGTH"];
    double alpha = params["PATCH_DECAY"];

    double pow_val = (1.0/(alpha+1.0));

    patches = new std::vector<Patch*>;
    std::vector<double> k_vals;
    double k, y;

    double k_sum = 0;
    for (int i = 0; i < n_patches; i++){
        y = real_uniform(0, 1.0, patch_generator);
        k = pow(y, pow_val);
        k_vals.push_back(k);
        k_sum += k;
    }

    for (int i = 0; i < n_patches; i++){
        k_vals[i] = k_vals[i]/k_sum;
        printf("k_val: %f\n", k_vals[i]);
    }
    printf("\n");
    for (int i = 0; i < n_patches; i++){
        double k = double(n_indiv)*k_vals[i];
        double x = real_uniform(0, side_len-1, patch_generator);
        double y = real_uniform(0, side_len-1, patch_generator);
        printf("K: %.2f\n", k);
        Patch *tmp = new Patch(x,y,k);
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
    migration_tracker = new MigrationTracker();
    genome_dict = new GenomeDict();
}

void initialize_individuals(){
    double k;
    int i;
    int n_indiv;
    Individual* indiv_i;

    for (Patch* patch_i : *patches){
        k = patch_i->get_K();
        for (i = 0; i < k; i++){
            indiv_i = new Individual(patch_i, false);
            patch_i->add_individual(indiv_i);
        }
    }
}


void initialize_genomes(){
    int n_loci = params["NUM_OF_LOCI"];

    std::vector<std::vector<double>> alleles = gen_alleles();
    std::vector<std::vector<double>> props = generate_allele_freq_from_beta(alleles);

    for (Patch* patch_i : *patches){
        int patch_size = patch_i->get_size();
        std::vector<std::vector<int>> num_per_allele = get_num_with_each_allele(props, patch_size);

        for (int locus = 0; locus < n_loci; locus++){
            std::vector<double> allele_map;
            int n_alleles = num_per_allele[locus].size();
            for (int i = 0; i < n_alleles; i++){
                int n = num_per_allele[locus][i];
                for (int j = 0; j < n; j++){
                    allele_map.push_back(alleles[locus][i]);
                }
            }

            random_shuffle(std::begin(allele_map), std::end(allele_map));

            int allele_ct = 0;
            for (Individual* indiv: patch_i->get_all_individuals()){
                indiv->set_locus(locus, 0, allele_map[allele_ct]);
                allele_ct++;
                indiv->set_locus(locus, 1, allele_map[allele_ct]);
                allele_ct++;
            }
        }
    }

    #if __DEBUG__
        for (Patch* patch_i : *patches){
            for (Individual* indiv_i : patch_i->get_all_individuals()){
                for (int locus = 0; locus < n_loci; locus++){
                    double al1 = indiv_i->get_locus(locus, 0);
                    double al2 = indiv_i->get_locus(locus, 1);
                    if (al1 < 0.0 || al1 > 1.0){
                        assert(0 && "why has god abandoned us\n");
                    }
                    if (al2 < 0.0 || al2 > 1.0){
                        assert(0 && "god dammit all to hell\n");
                    }
                }
            }
        }
    #endif
}

int expected_num_alleles(){
    std::uniform_int_distribution<int> num_allele_dis(15, 30);
    return num_allele_dis(*genome_generator);
}


/*  Patch::generate_allele_freq_from_beta(int n_alleles)
        Generates the allelic frequencies in the population for n_alleles alleles from a beta distribution
*/
std::vector<std::vector<double>> generate_allele_freq_from_beta(std::vector<std::vector<double>> alleles){
    int n_loci = params["NUM_OF_LOCI"];


    std::vector<std::vector<double>> props;

    for (int l = 0; l < n_loci; l++){
        props.push_back(std::vector<double>());
        int n_alleles = alleles[l].size();
        double rem = 1.0;
        double p, prop;
        for (int i = 0; i < n_alleles - 1; i++){
            p = beta_dist(0.6, 1.7, main_generator);
            prop = rem * p;
            rem -= prop;
            props[l].push_back(prop);
        }
        props[l].push_back(rem);
    }

    return props;
}

std::vector<std::vector<int>>  get_num_with_each_allele(std::vector<std::vector<double>> props, int patch_size){
    int n_loci = params["NUM_OF_LOCI"];
    int total_n_sites = 2*patch_size;

    std::vector<std::vector<int>> num_per_allele;


    for (int l = 0; l < n_loci; l++){
        num_per_allele.push_back(std::vector<int>());
        int num_als_this_locus = props[l].size();
        int n_left = total_n_sites;
        for (int al_num = 0; al_num < num_als_this_locus; al_num++){
            double prop = props[l][al_num];
            int n_this_al = int(total_n_sites*prop);
            n_left -= n_this_al;
            num_per_allele[l].push_back(n_this_al);
        }
        num_per_allele[l][0] += n_left;
    }
    return num_per_allele;
}

std::vector<std::vector<double>> gen_alleles(){
    std::vector<std::vector<double>> alleles;

    int n_loci = params["NUM_OF_LOCI"];

    for (int l = 0; l < n_loci; l++){
        alleles.push_back(std::vector<double>());

        int n_alleles = expected_num_alleles();

        for (int i = 0; i < n_alleles; i++){
            //alleles[l].push_back(1.0);
            alleles[l].push_back(real_uniform(0, 1.0, main_generator));
        }
    }
    return alleles;
}
