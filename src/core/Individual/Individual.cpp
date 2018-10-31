#include "Individual.h"
#include "Patch.h"
#include "GenomeDict.h"

int Individual::id_counter = 0;

Individual::Individual(Patch* patch){
    int genome_size = params["NUM_OF_LOCI"];
    int s = int_uniform(0,1, patch_generator);

    this->id = this->id_counter++;
    this->sex = s;
    this->patch = patch;
    this->has_migrated = false;
    this->haplotype0 = new double[genome_size];
    this->haplotype1 = new double[genome_size];

    for (int l = 0; l < genome_size; l++){
        this->haplotype0[l] = -1.0;
        this->haplotype1[l] = -1.0;
    }
}

int Individual::get_id(){
    return this->id;
}

void Individual::set_locus(int locus, int haplotype, double val){
    if (haplotype == 0){
        this->haplotype0[locus] = val;
    }
    else if (haplotype == 1){
        this->haplotype1[locus] = val;
    }
}

double Individual::get_locus(int locus, int haplotype){
    if (haplotype == 0){
        return this->haplotype0[locus];
    }
    else if (haplotype == 1){
        return this->haplotype1[locus];
    }
    assert(0 && "invalid haplotype");
}


// =========================================================================
// Migration
// =========================================================================

void Individual::migrate(){
    if (this->has_migrated){
        return;
    }

    std::vector<Patch*> delta = this->stochastic_foraging();
    Patch* new_patch = pick_best_patch(delta);
    if (new_patch){
        this->has_migrated = true;
        this->patch->remove_individual(this);
        new_patch->add_individual(this);
    }
}

std::vector<Patch*> Individual::stochastic_foraging(){
    int mean_dist = params["MEAN_MIGRATION_DIST"];
    int mean_num_patches = params["MEAN_NUM_PATCHES_FORAGED"];

    std::vector<Patch*> potential_patches;

    potential_patches.push_back(this->patch);

    int num_of_patches = poisson(mean_num_patches, main_generator);

    for (int i = 0; i < num_of_patches; i++){
        double crit_dist = poisson(mean_dist, main_generator);
        potential_patches.push_back(find_nearest_patch_to_crit(crit_dist));
    }

    std::set<Patch*> s( potential_patches.begin(), potential_patches.end() );
    std::vector<Patch*> potential_patches_unique;
    potential_patches_unique.assign( s.begin(), s.end() );

    return potential_patches_unique;
}

Patch* Individual::find_nearest_patch_to_crit(double crit_dist){
    double min = 100;
    double dist_from_crit;
    Patch* best;

    int id = this->patch->get_id();
    int id2;

    for (Patch* patch_j : *patches){
        id2 = patch_j->get_id();
        if (this->patch != patch_j){
            dist_from_crit = dist_matrix[id][id2] - crit_dist;
            if (dist_from_crit < min){
                min = dist_from_crit;
                best = patch_j;
            }
        }
    }

    return best;
}

Patch* Individual::pick_best_patch(std::vector<Patch*> options){
    Patch* best = this->patch;

    double max_pref = 0;
    double pref;

    for (Patch* patch_i : options){
        pref = this->calc_pref(patch_i);
        if (pref > max_pref){
            best = patch_i;
            max_pref = pref;
        }
    }

    return best;
}

double Individual::calc_pref(Patch* patch){
    int n_ef = params["NUM_ENV_FACTORS"];
    int n_loci_per_ef = params["NUM_LOCI_PER_EF"];

    int locus;
    double theta_i, y_i, s_i, p_i;

    double p = 1.0;
    double s_max_i = 0.1;

    std::vector<double> theta = this->patch->get_env_factors();
    for (int i = 0; i < n_ef; i++){
        theta_i = theta[i];

        for (int j = 0; j < n_loci_per_ef; j++ ){
            locus = genome_dict->fitness_loci[i][j];

            y_i = this->get_locus(locus, 0);
            s_i = s_max_i * exp(-(y_i - theta_i)*(y_i - theta_i)); // TODO make this draw from gaussian
            p_i = 1.0 + s_i;
            p = p * p_i;

            y_i = this->get_locus(locus, 1);
            s_i = s_max_i * exp(-(y_i - theta_i)*(y_i - theta_i));  // TODO make this draw from gaussian
            p_i = 1.0 + s_i;
            p = p * p_i;
        }
    }

    return p;
}


// =========================================================================
// Selection
// =========================================================================
