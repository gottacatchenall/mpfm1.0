#include "Individual.h"
#include "Patch.h"
#include "GenomeDict.h"
#include "MigrationTracker.h"

int Individual::id_counter = 0;

Individual::Individual(Patch* patch, bool parent_was_migrant){
    int genome_size = params["NUM_OF_LOCI"];
    int s = int_uniform(0,1, patch_generator);

    this->id = this->id_counter++;
    this->sex = s;
    this->patch = patch;
    this->patch_born_in = patch;
    this->has_migrated = false;
    this->haplotype0 = new double[genome_size];
    this->haplotype1 = new double[genome_size];

    this->parent_was_migrant = parent_was_migrant;
    /*for (int i = 0; i < genome_size; i++){
        this->haplotype0[i] = 1.0;
        this->haplotype1[i] = 1.0;
    }*/
}

Individual::~Individual(){
    delete this->haplotype0;
    delete this->haplotype1;
}

int Individual::get_id(){
    return this->id;
}

Patch* Individual::get_patch(){
    return this->patch;
}
Patch* Individual::get_patch_born_in(){
    return this->patch_born_in;
}

void Individual::set_locus(int locus, int haplotype, double val){
    if (haplotype == 0){
        this->haplotype0[locus] = val;
        return;
    }
    else if (haplotype == 1){
        this->haplotype1[locus] = val;
        return;
    }
    assert(0 && "invalid haplo!");
}

double Individual::get_locus(int locus, int haplotype){
    double al = -1;
    if (haplotype == 0){
        al =  this->haplotype0[locus];
    }
    else if (haplotype == 1){
        al = this->haplotype1[locus];
    }
    /*
    if (al > 1 || al < 0){
        assert(0 && "alleles still not working...");
    }*/
    return al;
    assert(0 && "invalid haplotype");
}


// =========================================================================
// Migration
// =========================================================================

void Individual::migrate(Patch* patch){
    this->patch = patch;
    migration_tracker->note_attempted_migration(this->patch_born_in, patch);
    if (patch != this->patch_born_in){
        this->has_migrated = true;
    }
}

void Individual::migrate_old(){
    if (this->has_migrated){
        return;
    }

    //std::vector<Patch*> delta = this->stochastic_foraging();
    //Patch* new_patch = pick_best_patch(delta);
    Patch* best_patch = this->patch;
    double distance_strength = 0.01;
    double pref, dist, decay, score;
    int current_patch_id = this->patch->get_id();


    int mean_num_patches = params["MEAN_NUM_PATCHES_FORAGED"];
    int num_of_patches = poisson(mean_num_patches, main_generator);
    int total_n_patches = patches->size();

    std::vector<Patch*> avail_patches;

    for (int i = 0; i < num_of_patches; i++){
        int index = int_uniform(0, total_n_patches-1, main_generator);
        avail_patches.push_back((*patches)[index]);
    }



    double best_score =  this->calc_pref(this->patch);;
    for (Patch* patch_i: avail_patches){
        dist = abs(dist_matrix[current_patch_id][patch_i->get_id()]);
        decay = exp(-1.0*distance_strength*(dist)*(dist));

        pref = this->calc_pref(patch_i);

        score = pref*decay;

        if (score > best_score){
            best_score = score;
            best_patch = patch_i;
        }
    }

    if (best_patch != this->patch){
        migration_tracker->note_attempted_migration(this->patch, best_patch);
        this->has_migrated = true;
        this->patch->remove_individual(this);
        //best_patch->add_to_migrant_queue(this);
        best_patch->add_individual(this);
    }
    else{
        migration_tracker->note_attempted_migration(this->patch, this->patch);
    }
}

// =========================================================================
// Selection
// =========================================================================
double Individual::calc_pref(Patch* patch_i){
    int n_ef = params["NUM_ENV_FACTORS"];
    int n_loci_per_ef = params["NUM_LOCI_PER_EF"];
    double sigma_s = params["SIGMA_SELECTION"];

    double sel_strength = params["SELECTION_STRENGTH"];

    int locus;
    double theta_i, x_i, s_i, w_i;

    double w = 1.0;
    double s_max_i = -1*sel_strength;

    double draw;
    double locus_weight;

    std::vector<double> theta = patch_i->get_env_factors();
    for (int i = 0; i < n_ef; i++){
        theta_i = theta[i];
        for (int j = 0; j < n_loci_per_ef; j++ ){
            locus = genome_dict->fitness_loci[i][j];

            locus_weight = genome_dict->selection_strengths[locus];

            x_i = this->get_locus(locus, 0);
            draw = abs(locus_weight*normal((x_i - theta_i), sigma_s, main_generator));
            s_i = s_max_i * draw;
            w_i = 1.0 + s_i;
            w = w * w_i;

            x_i = this->get_locus(locus, 1);
            draw = abs(locus_weight*normal((x_i - theta_i), sigma_s, main_generator));
            s_i = s_max_i * draw;
            w_i = 1.0 + s_i;
            w = w * w_i;
        }
    }

    return w;
}



void Individual::calc_fitness(){
    int n_ef = params["NUM_ENV_FACTORS"];
    int n_loci_per_ef = params["NUM_LOCI_PER_EF"];
    double sigma_s = params["SIGMA_SELECTION"];
    double sel_strength = params["SELECTION_STRENGTH"];

    int locus;
    double theta_i, x_i, s_i, w_i;

    double w = 1.0;
    double s_max_i = -1*sel_strength;

    double draw;
    double locus_weight;

    std::vector<double> theta = this->patch->get_env_factors();
    for (int i = 0; i < n_ef; i++){
        theta_i = theta[i];
        for (int j = 0; j < n_loci_per_ef; j++ ){
            locus = genome_dict->fitness_loci[i][j];

            locus_weight = genome_dict->selection_strengths[locus];

            x_i = this->get_locus(locus, 0);
            draw = abs(locus_weight*normal((x_i - theta_i), sigma_s, main_generator));
            s_i = s_max_i * draw;
            w_i = 1.0 + s_i;
            w = w * w_i;

            x_i = this->get_locus(locus, 1);
            draw = abs(locus_weight*normal((x_i - theta_i), sigma_s, main_generator));
            s_i = s_max_i * draw;
            w_i = 1.0 + s_i;
            w = w * w_i;
        }
    }

    this->w = w;
}

double Individual:: get_fitness(){
    return this->w;
}

void Individual::set_exp_num_off(double val){
    this->exp_num_off = val;
}
double Individual::get_exp_num_off(){
    return this->exp_num_off;
}



int Individual::get_sex(){
    return this->sex;
}

void Individual::gen_haplotype(Individual* parent, int offspring_haplo){

    std::vector<int> crossing_over_points = get_crossing_over_points();

    int curr_crossover = 0;
    int curr_chromo = 0;
    int current_haplo = int_uniform(0,1, main_generator);
    int n_loci = params["NUM_OF_LOCI"];

    for (int locus = 0; locus < n_loci; locus++){
        // Check if at start of a new chromosome
        if (locus == genome_dict->chromo_map[curr_chromo]){
            curr_chromo++;
            current_haplo = int_uniform(0,1, main_generator);
        }

        if (crossing_over_points.size() > 0){
            if (locus == crossing_over_points[curr_crossover]){
                curr_crossover++;
                current_haplo = !(current_haplo);
            }
        }


        bool sets_something = false;

        if (real_uniform(0,1, main_generator) < params["MUTATION_RATE"]){
            double new_al = real_uniform(0,1, main_generator);
            assert(new_al > 0 && new_al < 1);
            sets_something = true;
            this->set_locus(locus, offspring_haplo, new_al);
        }

        else{
            double new_al =  parent->get_locus(locus, current_haplo);
            assert(new_al >= 0 && new_al <= 1);
            sets_something = true;
            this->set_locus(locus, offspring_haplo, new_al);
        }

        assert(sets_something);
    }
}

std::vector<int> Individual::get_crossing_over_points(){
    int n_loci = params["NUM_OF_LOCI"];
    int cM = params["GENOME_LENGTH_CENTIMORGANS"];
    int n_crossover_events = poisson(cM, main_generator);

    std::vector<int> crossing_over_points;

    // Generate a random permutation of unique loci to be the crossover points
    int r;
    bool exists;
    for(int i = 0; i < n_crossover_events; i++){
        do{
            r = int_uniform(0,n_loci-1, main_generator);
            exists = (std::find(crossing_over_points.begin(), crossing_over_points.end(), r) != crossing_over_points.end());
        } while(exists);
        crossing_over_points.push_back(r);
    }

    // Sort crossing_over_points
    std::sort(crossing_over_points.begin(), crossing_over_points.end());
    return crossing_over_points;
}
