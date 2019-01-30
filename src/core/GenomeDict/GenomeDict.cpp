#include "GenomeDict.h"

GenomeDict::GenomeDict(){
    int n_loci = params["NUM_OF_LOCI"];
    int n_chromo = params["NUM_OF_CHROMOSOMES"];
    int n_ef = params["NUM_ENV_FACTORS"];
    int n_loci_per_ef = params["NUM_LOCI_PER_EF"];

    int* chromo = this->generate_perm_with_uniq_ints(n_chromo, n_loci);

    for (int i = 0; i < n_chromo; i++){
        this->chromo_map.push_back(chromo[i]);
    }


    int n_unique_indecies = n_loci_per_ef*n_ef;
    int* uniq_perm = generate_perm_with_uniq_ints(n_unique_indecies, n_loci);

    int ind = 0;
    double weight;

    this->selection_strengths = new double[n_loci];
    for (int i = 0; i < n_loci; i++){
        this->selection_strengths[i] = 0.0;
    }

    double MEAN_LOCUS_WEIGHT = params["MEAN_LOCUS_WEIGHT"];

    for (int i = 0; i < n_ef; i++){
        std::vector<int> tmp_vec;
        for (int j = 0; j < n_loci_per_ef; j++){
            int val = uniq_perm[ind];
            //weight = real_uniform(0.0, 1.0, genome_generator);
            weight = MEAN_LOCUS_WEIGHT;
            this->selection_strengths[val] = (weight);

            log_locus(val, i, weight, "fitness");
            tmp_vec.push_back(uniq_perm[ind]);
            ind++;
        }

        this->fitness_loci.push_back(tmp_vec);
    }

    for (int i = 0; i < n_loci; i++){
        bool neut = true;
        for (int j = 0; j < n_unique_indecies; j++){
            if (i == uniq_perm[j]){
                neut = false;
            }
        }
        if (neut){
            this->neutral_loci.push_back(i);
            log_locus(i, -1, 0, "neutral");
        }
    }

    assert(this->chromo_map.size() == n_chromo);
    free(chromo);
    free(uniq_perm);
}


// Generates a permutation of length size of integers from 0...n-1
int* GenomeDict::generate_perm_with_uniq_ints(int size, int n){
    int* result = new int[size];
    int r;
    int exists;
    for(int i = 0; i < size; i++){
        do{
            exists = 0;
            r = int_uniform(0,n-1, genome_generator);
            for (int j = 0; i < j; j++){
                if (result[j] == r){
                    exists = 1;
                    break;
                }
            }
        } while(exists);
        result[i] = r;
    }
    return result;
}
