
#include "Patch.h"
#include "Individual.h"
#include "EnvFactor.h"

int Patch::id_counter = 0;

Patch::Patch(double x, double y, double K){
    this->x = x;
    this->y = y;
    this->K = K;
    this->id = id_counter++;
}

double Patch::get_x(){
    return this->x;
}

double Patch::get_y(){
    return this->y;
}

double Patch::get_K(){
    return this->K;
}

void Patch::add_individual(Individual* indiv){
    this->individuals.insert(std::make_pair(indiv->get_id(), indiv));
}

void Patch::remove_individual(Individual* indiv){
    this->individuals.erase(indiv->get_id());
}

int Patch::get_size(){
    return this->individuals.size();
}

int Patch::get_id(){
    return this->id;
}

std::vector<Individual*> Patch::get_all_individuals(){
   std::vector<Individual*> indiv;
   for (auto ind: this->individuals){
       indiv.push_back(ind.second);
   }

   return indiv;
}

std::vector<double> Patch::get_env_factors(){
    double x = this->x;
    double y = this->y;

    std::vector<double> ef_vector;

    for (EnvFactor* ef_i : *envFactors){
        ef_vector.push_back(ef_i->get_env_factor_value(x,y));
    }

    return ef_vector;
}
