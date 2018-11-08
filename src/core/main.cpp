#include "include.h"
#include "Patch.h"
#include "Individual.h"
#include "AlleleTracker.h"
#include "MigrationTracker.h"

int generation;
MigrationTracker* migration_tracker;
std::map<std::string, double> params;
std::vector<EnvFactor*>* envFactors;
std::vector<Patch*>* patches;
std::vector<std::vector<double>> dist_matrix;
GenomeDict* genome_dict;

int main(int argc, char* argv[]){

    init();

    int num_gen = params["NUM_GENERATIONS"];
    int census_freq = params["CENSUS_FREQ"];

    for (generation = 0; generation <= num_gen; generation++){
        migration();
        selection();
        check_dispersal();
        logging();
        mating();
        if (generation % census_freq == 0 && generation > 0){
            census();
        }
        update_progress_bar(generation);
        migration_tracker->reset_migration_matrix();
    }

    return 0;
}

void migration(){
    for (Patch* patch_i: *patches){
        for (Individual* indiv: patch_i->get_all_individuals() ){
            indiv->migrate();
        }
    }
}

void selection(){
    for (Patch* patch_i: *patches){
        for (Individual* indiv: patch_i->get_all_individuals() ){
            indiv->calc_fitness();
        }
        patch_i->selection();
    }
}

void logging(){
    for (Patch* patch_i: *patches){
        if (patch_i->get_size() > 0){
            log_population(patch_i);
            log_eff_migration(patch_i);
            log_patch(patch_i);
            for (Patch* patch_j: *patches){
                if (patch_j->get_size()){
                    log_attempted_migration(patch_i, patch_j);
                    log_successful_migration(patch_i, patch_j);
                }
            }
        }
    }
    log_fst(get_fst());
}

void mating(){
    int mean_off_per_female = params["AVG_NUM_OFFSPRING_PER_FEMALE"];

    for (Patch* patch_i: *patches){

        std::vector<std::vector<Individual*>> smap = patch_i->split_by_sex();
        std::vector<Individual*> females = smap[0];
        std::vector<Individual*> males = smap[1];

        Individual* random_male;
        Individual* offspring;
        int num_males = males.size();
        if (num_males > 0){
            for (Individual* female: females){
                int n_off = poisson(mean_off_per_female, main_generator);
                int index = int_uniform(0, males.size()-1, main_generator);
                random_male = males[index];

                bool parent_migrated = (random_male->has_migrated || female->has_migrated);

                for (int i = 0; i < n_off; i++){
                    offspring = new Individual(patch_i, parent_migrated);
                    offspring->gen_haplotype(female, 0);
                    offspring->gen_haplotype(random_male, 1);

                    #if __DEBUG__
                        int n_loci = params["NUM_OF_LOCI"];
                        for (int i = 0; i < n_loci; i++){
                            double al0 = offspring->get_locus(i, 0);
                            double al1 = offspring->get_locus(i, 1);
                            assert(al0 >= 0 && al0 <= 1 && al1 >= 0 && al1 <= 1);
                        }
                    #endif

                    patch_i->add_to_next_gen(offspring);
                }
            }
        }
        patch_i->replace_current_gen();
    }
}

void check_dispersal(){
    for (Patch* patch_i: *patches){
        for (Individual* indiv: patch_i->get_all_individuals()){
            migration_tracker->note_successful_migration(patch_i, indiv->get_patch_born_in());
        }
    }
}

void census(){
    AlleleTracker al_tracker;
    al_tracker.construct_allele_table();
    for (Patch* patch_i: *patches){
        al_tracker.get_ld(patch_i->get_id(), "fitness");
        al_tracker.get_ld(patch_i->get_id(), "pref");
        al_tracker.get_ld(patch_i->get_id(), "neutral");
        al_tracker.get_global_ld("fitness");
        al_tracker.get_global_ld("pref");
        al_tracker.get_global_ld("neutral");
    }
}


std::vector<double> get_fst(){
    std::vector<int> pooled_heterozygos_ct_by_locus;
    int n_loci = params["NUM_OF_LOCI"];

    for (int l = 0; l < n_loci; l++){
        pooled_heterozygos_ct_by_locus.push_back(0);
    }


    std::vector<double> hS_sum;
    for (int l = 0; l < n_loci; l++){
        hS_sum.push_back(0.0);
    }

    int n_patches = 0;

    for (Patch* patch_i : *patches){
        int n_patch = patch_i->get_size();
        if (n_patch > 0){
            n_patches++;

            std::vector<int> patch_heterozygos_ct_by_locus;
            for (int l = 0; l < n_loci; l++){
                patch_heterozygos_ct_by_locus.push_back(0);
            }

            for (Individual* indiv: patch_i->get_all_individuals()){
                for (int l = 0; l < n_loci; l++){
                    if (indiv->get_locus(l, 0) != indiv->get_locus(l,1)){
                        pooled_heterozygos_ct_by_locus[l]++;
                        patch_heterozygos_ct_by_locus[l]++;
                    }
                }
            }

            for (int l = 0; l < n_loci; l++){
                double patch_hS = double(patch_heterozygos_ct_by_locus[l])/double(n_patch);
                hS_sum[l] += patch_hS;
            }
        }
    }

    std::vector<double> f_st;
    int n_total = get_total_population_size();
    if (n_total > 0){
        double hT;
        double hS;


        for (int l = 0; l < n_loci; l++){
            assert(n_patches);
            assert(n_total);

            hS = double(hS_sum[l]) / double(n_patches);
            hT = double(pooled_heterozygos_ct_by_locus[l])/double(n_total);
            if (hT > 0){
                double FST = (hT - hS) / hT;
                f_st.push_back(FST);
            }
            else{
                f_st.push_back(1.0);
            }
        }
    }

    else{
        for (int l = 0; l < n_loci; l++){
            f_st.push_back(0.0);
        }
    }


    return f_st;
}

int get_total_population_size(){
    int ct = 0;
    for (Patch* patch_i : *patches){
        ct += patch_i->get_size();
    }
    return ct;
}


void update_progress_bar(int gen){
    int n_gen = params["NUM_GENERATIONS"];

    double progress = double(gen)/double(n_gen);
    int barWidth = 50;

       std::cout << "[";
       int pos = barWidth * progress;
       for (int i = 0; i < barWidth; ++i) {
           if (i <= pos) std::cout << "=";
           else std::cout << " ";
       }
       std::cout << "] [" << gen << " / " << n_gen << "] \r";
       std::cout.flush();

}
