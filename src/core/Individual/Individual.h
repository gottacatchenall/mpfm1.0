
#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "include.h"

class Individual{
    private:
        static int id_counter;
        int id;
        int sex;
        Patch* patch_born_in;
        Patch* patch;
        double* haplotype0;
        double* haplotype1;
        double w;
    public:
        bool has_migrated;
        double parent_was_migrant;
        double exp_num_off;

        Individual(Patch* patch, bool parent_was_migrant);
        ~Individual();
        Patch* get_patch();
        Patch* get_patch_born_in();
        int get_id();
        int get_sex();
        double get_fitness();
        void set_locus(int locus, int haplotype, double val);
        double get_locus(int locus, int haplotype);

        void migrate(Patch* patch);
        void migrate_old();

        void set_exp_num_off(double val);
        double get_exp_num_off();

        void calc_fitness();
        double calc_pref(Patch* patch_i);
        void gen_haplotype(Individual* parent, int offspring_haplo);
        std::vector<int> get_crossing_over_points();

};

#endif
