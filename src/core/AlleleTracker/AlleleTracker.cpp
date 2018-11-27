
#include "AlleleTracker.h"
#include "Patch.h"
#include "Individual.h"
#include "GenomeDict.h"

AlleleTracker::AlleleTracker(){
    int n_loci = params["NUM_OF_LOCI"];
    allele_map = new std::vector<allele*>[n_loci];
}


void AlleleTracker::get_ld(int patch_num, std::string type){
    std::vector<allele*> alleles;

    double f_al1, f_al2, f_both, ld;

    Patch* p = (*patches)[patch_num];
    int n_total = p->get_size();

    int n_ef = params["NUM_ENV_FACTORS"];
    int n_loci_per_ef = params["NUM_LOCI_PER_EF"];


    std::vector<int> loci;
    if (type == "fitness"){
        for (int i = 0; i < n_ef; i++){
            for (int j = 0; j < n_loci_per_ef; j++){
                loci.push_back(genome_dict->fitness_loci[i][j]);
            }
        }
    }
    else if (type == "neutral"){
        for (int l : genome_dict->neutral_loci){
            loci.push_back(l);
        }
    }

    std::sort(loci.begin(), loci.end());



    int n_loci = loci.size();

    if (n_total > 0){
        for (int l1 = 0; l1 < n_loci; l1++){
            int l1_i = loci[l1];
            alleles = this->allele_map[l1_i];
            for (allele* al1: alleles){
                f_al1 = double(al1->freq_map[patch_num])/double(2*n_total);
                if (f_al1 > 0.0){
                    log_allele_freq(patch_num, al1->locus, al1->allele_val, f_al1, type);
                }
            }
            for (int l2 = l1+1; l2 < n_loci; l2++){
                int l2_i = loci[l2];
                double sum = 0.0;
                int ct = 0;

                for (allele* al1: alleles){
                    f_al1 = double(al1->freq_map[patch_num])/double(2*n_total);
                    if (f_al1 > 0){
                        for (dependent_allele* al2 : al1->loci[l2_i]){
                            f_both = double(al2->freq_map[patch_num])/double(2*n_total);
                            f_al2 = double(this->find_allele(al2->locus, al2->allele_val)->freq_map[patch_num])/double(2*n_total);

                            if (f_al2 > 0){
                                if (!(f_al1 >= 0.0 && f_al1 <= 1.0) || !((f_al2 >= 0.0 && f_al2 <= 1.0)) || !(f_both >= 0.0 && f_both <= 1.0)){
                                    printf("al1_n = %d\n", al1->freq_map[patch_num]);
                                    printf("al2_n = %d\n", this->find_allele(al2->locus, al2->allele_val)->freq_map[patch_num]);
                                    printf("both_n = %d\n", al2->freq_map[patch_num]);
                                    printf("ld assert\n");
                                    exit(-1);
                                }

                                ld = (f_al1 * f_al2) - f_both;
                                sum += abs(ld);
                                ct++;
                                //log_linkage(patch_num, l1, al1->allele_val, l2, al2->allele_val, ld, type);
                            }
                        }
                    }
                }
                if (ct == 0){
                    // only happens when l2 = l1+1
                    //printf("somehow ct ended up at 0?\n");
                }
                double avg_ld_this_pair = double(sum)/double(ct);
                log_linkage(patch_num, l1, l2, avg_ld_this_pair, type);
            }
        }
    }
}

void AlleleTracker::get_global_ld(std::string type){

    double f_al1, f_al2, f_both, ld;
    int n_ef = params["NUM_ENV_FACTORS"];
    int n_loci_per_ef = params["NUM_LOCI_PER_EF"];

    std::vector<int> loci;
    if (type == "fitness"){
        for (int i = 0; i < n_ef; i++){
            for (int j = 0; j < n_loci_per_ef; j++){
                loci.push_back(genome_dict->fitness_loci[i][j]);
            }
        }
    }
    else if (type == "neutral"){
        for (int l : genome_dict->neutral_loci){
            loci.push_back(l);
        }
    }

    int n_loci = loci.size();
    int n_total = get_total_population_size();
    std::sort(loci.begin(), loci.end());

    std::vector<allele*> alleles;

    for (int l1 = 0; l1 < n_loci; l1++){
        int l1_i = loci[l1];

        for (int l2 = l1+1; l2 < n_loci; l2++){
            alleles = this->allele_map[l1_i];
            int l2_i = loci[l2];
            double sum = 0.0;
            int ct = 0;

            for (allele* al1: alleles){
                f_al1 = double(al1->n_total)/double(2*n_total);
                if (f_al1 > 0){
                    for (dependent_allele* al2 : al1->loci[l2_i]){
                        f_both = double(al2->n_total)/double(2*n_total);
                        f_al2 = double(this->find_allele(al2->locus, al2->allele_val)->n_total)/double(2*n_total);

                        if (f_al2 > 0){
                            if (!(f_al1 >= 0.0 && f_al1 <= 1.0) || !((f_al2 >= 0.0 && f_al2 <= 1.0)) || !(f_both >= 0.0 && f_both <= 1.0)){
                                printf("al1_n = %d\n", al1->n_total);
                                printf("al2_n = %d\n", this->find_allele(al2->locus, al2->allele_val)->n_total);
                                printf("both_n = %d\n", al2->n_total);
                                printf("ld assert\n");
                                exit(-1);
                            }

                            ld = (f_al1 * f_al2) - f_both;
                            sum += abs(ld);
                            ct++;
                            //log_linkage(patch_num, l1, al1->allele_val, l2, al2->allele_val, ld, type);
                        }
                    }
                }
            }
            if (ct == 0){
                // only happens when l2 = l1+1
                //printf("somehow ct ended up at 0?\n");
            }
            double avg_ld_this_pair = double(sum)/double(ct);
            log_global_linkage(l1, l2, avg_ld_this_pair, type);
        }
    }
}

void AlleleTracker::construct_allele_table(){
    assert(allele_map[0].size() == 0);

    int n_loci = params["NUM_OF_LOCI"];
    allele* al1_1;
    allele* al1_2;

    for (Patch* patch_i: *patches){
        int patch_id = patch_i->get_id();
        for (Individual* indiv: patch_i->get_all_individuals()){
            for (int l1 = 0; l1 < n_loci; l1++){
                this->update_tracker(l1, indiv->get_locus(l1, 0), patch_id);
                this->update_tracker(l1, indiv->get_locus(l1, 1), patch_id);
            }
            for (int l1 = 0; l1 < n_loci; l1++){
                al1_1 = find_allele(l1, indiv->get_locus(l1, 0));
                al1_2 = find_allele(l1, indiv->get_locus(l1, 1));
                for (int l2 = l1; l2 < n_loci; l2++){
                    this->add_dependent_allele(al1_1, l2, indiv->get_locus(l2, 0), patch_id);
                    //this->add_dependent_allele(al1_1, l2, indiv->get_locus(l2, 1), patch_id);
                    //this->add_dependent_allele(al1_2, l2, indiv->get_locus(l2, 0), patch_id);
                    this->add_dependent_allele(al1_2, l2, indiv->get_locus(l2, 1), patch_id);
                }
            }
        }
    }
}

void AlleleTracker::update_tracker(int locus, double allele_val, int patch_id){
    int n_loci = params["NUM_OF_LOCI"];
    allele* allele_struct = find_allele(locus, allele_val);
    if (allele_struct != NULL){
        note_allele_seen(allele_struct, patch_id);
    }
    else{
        allele* new_allele = new allele(locus, allele_val, n_loci);
        note_allele_seen(new_allele, patch_id);
        this->allele_map[locus].push_back(new_allele);
    }
}

void AlleleTracker::add_dependent_allele(allele* primary_allele, int dependent_locus, double dependent_allele_val, int patch_id){
    assert(primary_allele != NULL);

    for (dependent_allele* al2 : primary_allele->loci[dependent_locus]){
        if (al2->allele_val == dependent_allele_val){
           note_allele_seen(al2, patch_id);
           return;
        }
    }


    dependent_allele* dep_al = new dependent_allele(dependent_locus, dependent_allele_val);
    primary_allele->loci[dependent_locus].push_back(dep_al);

    note_allele_seen(dep_al, patch_id);
}

void AlleleTracker::note_allele_seen(dependent_allele* al, int patch_id){
    al->n_total++;
    al->freq_map[patch_id]++;
}

void AlleleTracker::note_allele_seen(allele* al, int patch_id){
    al->n_total++;
    al->freq_map[patch_id]++;
}

allele* AlleleTracker::find_allele(int locus, double allele_val){
    std::vector<allele*> allele_vec = this->allele_map[locus];
    for(allele* al: allele_vec) {
        if (al->allele_val == allele_val){
            return al;
        }
    }
    return NULL;
}
