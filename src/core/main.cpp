#include "include.h"
#include "Patch.h"
#include "Individual.h"
#include "AlleleTracker.h"
#include "MigrationTracker.h"
#include "EnvFactor.h"
#include "GenomeDict.h"

int generation;
MigrationTracker* migration_tracker;
std::map<std::string, double> params;
std::vector<EnvFactor*>* envFactors;
std::vector<Patch*>* patches;
std::vector<std::vector<double>> dist_matrix;
GenomeDict* genome_dict;

int main(int argc, char* argv[]){
    if (argv[1]){
        chdir(argv[1]);
    }


    init();
    int num_gen = params["NUM_GENERATIONS"];
    int census_freq = params["CENSUS_FREQ"];
    int samp_freq = params["MIGRATION_SAMPLE_FREQ"];
    bool static_env_factors = params["ENV_FACTORS_STATIC"];

    for (generation = 0; generation <= num_gen; generation++){
        if (!static_env_factors){
            shift_environment();
        }

        migration();
        selection();
        check_dispersal();

        if (generation % samp_freq == 0 && generation > 0){
            logging();
        }


        mating();

        if (generation % census_freq == 0 && generation > 0){
            census();
            write_efs();
        }

        update_progress_bar(generation);
        migration_tracker->reset_migration_matrix();
    }
    return 0;
}

void migration(){

/*        for (Patch* patch_i: *patches){
            std::vector<Individual*> indivs = patch_i->get_all_individuals();
            for (Individual* indiv: indivs){
                indiv->migrate_old();
            }
        }*/

    double base_mig = params["BASE_MIGRATION_RATE"];

    std::vector<Individual*> indivs;
    for (Patch* patch_i: *patches){
        std::vector<double> row = migration_tracker->get_dispersal_row(patch_i->get_id());


        int n_patches = row.size();
        int n_this_patch = patch_i->get_size();
        int this_patch_id = patch_i->get_id();

        for (int j = 0; j < n_patches; j++){
            if (j != this_patch_id){
                Patch* patch_j = (*patches)[j];
                double prop = row[j];
                int n_i_to_j = int(prop * n_this_patch * base_mig);
                if (n_i_to_j > 0){
                    indivs = patch_i->pick_n_random_indivs(n_i_to_j);
                    for (Individual* indiv_i: indivs){
                        patch_j->add_to_migrant_queue(indiv_i);
                        indiv_i->migrate(patch_j);
                    }
                }
            }
        }
    }

    for (Patch* patch_i: *patches){
        patch_i->add_migrants_to_patch();
    }
}

void selection(){
    for (Patch* patch_i: *patches){
        std::vector<Individual*> indivs = patch_i->get_all_individuals();
        for (Individual* indiv: indivs){
            indiv->calc_fitness();
        }
    }

    for (Patch* patch_i: *patches){
        patch_i->selection();
    }
}

void logging(){
    for (Patch* patch_i: *patches){
        if (patch_i->get_size() > 0){
            int i = patch_i->get_id();
            log_eff_migration(patch_i);
            log_patch(patch_i);
            for (Patch* patch_j: *patches){
                int j = patch_j->get_id();
                if (patch_j->get_size() && i != j){
                    //log_attempted_migration(patch_i, patch_j);
                    log_successful_migration(patch_i, patch_j);
                }
            }
        }
    }
    //get_fst();
}


void write_efs(){
    int size = params["SIDE_LENGTH"];
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            log_env_factors(i,j);
        }
    }
}

void mating(){
    int mean_off_per_female = params["AVG_NUM_OFFSPRING_PER_FEMALE"];

    for (Patch* patch_i: *patches){

        std::vector<std::vector<Individual*>> smap = patch_i->split_by_sex();
        std::vector<Individual*> females = smap[0];
        std::vector<Individual*> males = smap[1];

        Individual* random_male;
        Individual* random_female;


        Individual* offspring;
        int num_males = males.size();
        int num_females = females.size();

        if (num_males > 0 && num_females > 0){
                int m_index, f_index;
                bool parent_migrated;

                for (int i = 0; i < num_females; i++){
                    random_female = females[i];
                    m_index = int_uniform(0, num_males-1, main_generator);
                    random_male = males[m_index];

                    parent_migrated = (random_male->has_migrated || random_female->has_migrated);

                    double n_off_mean = (random_male->get_exp_num_off() + random_female->get_exp_num_off())/double(2);

                    int n_off = poisson(n_off_mean, main_generator);

                    for (int off = 0; off < n_off; off++){
                        offspring = new Individual(patch_i, parent_migrated);
                        offspring->gen_haplotype(random_female, 0);
                        offspring->gen_haplotype(random_male, 1);
                        patch_i->add_to_next_gen(offspring);
                    }


                }

                /*
                for (int i = 0; i < n_off; i++){
                    m_index = int_uniform(0, num_males-1, main_generator);
                    f_index = int_uniform(0, num_females-1, main_generator);
                    random_male = males[m_index];
                    random_female = females[f_index];

                    parent_migrated = (random_male->has_migrated || random_female->has_migrated);

                    offspring = new Individual(patch_i, parent_migrated);
                    offspring->gen_haplotype(random_female, 0);
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
                }*/
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
        int p1 = patch_i->get_id();
        al_tracker.get_ld(p1, "");
        al_tracker.get_ld(p1, "");

        for (Patch* patch_j: *patches){
            int p2 = patch_j->get_id();
            if (p1 > p2){
                al_tracker.get_pairwise_ld(p1, p2, "");
                al_tracker.get_pairwise_ld(p1, p2, "");
            }
        }

    }
    al_tracker.get_global_ld("");
    al_tracker.get_global_ld("");
}

void shift_environment(){
    int num_gen = params["NUM_GENERATIONS"];

    double prop_shifted = 0.5 + (0.5 * double(generation)/double(num_gen));

    for (EnvFactor* ef : *envFactors){
        ef->shift(prop_shifted);
    }
}

void get_fst(){
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

    double hT;
    double hS;

    int n_total = get_total_population_size();

    std::vector<int> neut_loci = genome_dict->neutral_loci;
    std::vector<std::vector<int>> fit_loci = genome_dict->fitness_loci;



    for (int l : neut_loci){
        hS = double(hS_sum[l]) / double(n_patches);
        hT = double(pooled_heterozygos_ct_by_locus[l])/double(n_total);
        if (hT > 0){
            double FST = (hT - hS) / hT;
            log_fst(l, FST, "neutral");
        }
        else{
            log_fst(l, 1.0, "neutral");
        }
    }

    for (std::vector<int> fitness_loci_this_ef : fit_loci){
        for (int l : fitness_loci_this_ef){
            hS = double(hS_sum[l]) / double(n_patches);
            hT = double(pooled_heterozygos_ct_by_locus[l])/double(n_total);
            if (hT > 0){
                double FST = (hT - hS) / hT;
                log_fst(l, FST, "fitness");
            }
            else{
                log_fst(l, 1.0, "fitness");
            }
        }
    }

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
   std::cout << "\t[";
   int pos = barWidth * progress;
   for (int i = 0; i < barWidth; ++i) {
       if (i <= pos) std::cout << "=";
       else std::cout << " ";
   }
   std::cout << "] [" << gen << " / " << n_gen << "] \r";
   std::cout.flush();
}
