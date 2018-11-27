
#ifndef PATCH_H
#define PATCH_H

#include "include.h"

class Patch{
    private:
        static int id_counter;
        int id;
        int x;
        int y;
        double K;
        std::vector<Individual*>* individuals;
        std::vector<Individual*>* next_gen;
        std::vector<double> env_factors;
    public:
        Patch(int x, int y, double K);
        double get_x();
        double get_y();
        double get_K();
        int get_id();
        int get_size();
        void add_individual(Individual* indiv);
        void remove_individual(Individual* indiv);
        std::vector<Individual*> get_all_individuals();
        std::vector<double> get_env_factors();

        void selection();
        double beverton_holt_prob(int n, double k_prime);
        std::vector<std::vector<Individual*>> split_by_sex();

        void add_to_next_gen(Individual* indiv);
        void replace_current_gen();
};

#endif
