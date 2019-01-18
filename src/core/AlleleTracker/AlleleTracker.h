
#ifndef ALLELE_TRACKER_H
#define ALLELE_TRACKER_H

#include "include.h"

typedef struct dependent_allele{
    dependent_allele(int l, double val){
        locus = l;
        allele_val = val;
        n_total = 0;
    };
    int locus;
    double allele_val;
    int n_total;
    std::unordered_map<int,int> freq_map;
} dependent_allele;

typedef struct allele{
    allele(int l, double val, int n_loci){
        locus = l;
        allele_val = val;
        n_total = 0;

        loci = new std::vector<dependent_allele*>[n_loci];
    };

    int locus;
    double allele_val;
    int n_total;
    std::unordered_map<int,int> freq_map;
    std::vector<dependent_allele*>* loci;
} allele;

class AlleleTracker{
    private:
        std::vector<allele*>* allele_map;
    public:
        AlleleTracker();
        void construct_allele_table();
        void get_ld(int patch_num, std::string type);
        void get_pairwise_ld(int patch1_num, int patch2_num, std::string type);
        void get_global_ld(std::string type);
        void update_tracker(int locus, double allele_val, int patch_id);
        void add_dependent_allele(allele* primary_allele, int dependent_locus, double dependent_allele_val, int patch_id);
        void note_allele_seen(dependent_allele* al, int patch_id);
        void note_allele_seen(allele* al, int patch_id);
        allele* find_allele(int locus, double allele_val);
};

#endif
