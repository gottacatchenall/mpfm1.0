
#ifndef PATCH_H
#define PATCH_H

#include "include.h"

class Patch{
    private:
        static int id_counter;
        int id;
        double x;
        double y;
        double K;
        std::unordered_map<int,Individual*> individuals;
        std::vector<double> env_factors;
    public:
        Patch(double x, double y, double K);
        double get_x();
        double get_y();
        double get_K();
        int get_id();
        int get_size();
        void add_individual(Individual* indiv);
        void remove_individual(Individual* indiv);
        std::vector<Individual*> get_all_individuals();
        std::vector<double> get_env_factors();
};

#endif
