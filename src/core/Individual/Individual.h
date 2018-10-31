
#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "include.h"

class Individual{
    private:
        static int id_counter;
        int id;
        int sex;
        Patch* patch;
        double* haplotype0;
        double* haplotype1;
        bool has_migrated;
    public:
        Individual(Patch* patch);
        Patch* get_patch();
        int get_id();
        void set_locus(int locus, int haplotype, double val);
        double get_locus(int locus, int haplotype);

        void migrate();
        std::vector<Patch*> stochastic_foraging();
        Patch* find_nearest_patch_to_crit(double crit_dist);
        Patch* pick_best_patch(std::vector<Patch*> options);

        double calc_pref(Patch* patch);

};

#endif
