#include "include.h"
#include "Patch.h"
#include "Individual.h"
#include "EnvFactor.h"

int Patch::id_counter = 0;

Patch::Patch(double x, double y, double K){
    this->x = x;
    this->y = y;
    this->K = K;
    this->id = id_counter++;

    this->individuals = new std::vector<Individual*>;
    this->next_gen = new std::vector<Individual*>;
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
    this->individuals->push_back(indiv);
}

void Patch::remove_individual(Individual* indiv){
    int index = std::distance(this->individuals->begin(), std::find(this->individuals->begin(), this->individuals->end(), indiv));
    this->individuals->erase(this->individuals->begin()+index);
}

int Patch::get_size(){
    return this->individuals->size();
}

int Patch::get_id(){
    return this->id;
}

std::vector<Individual*> Patch::get_all_individuals(){
   return (*this->individuals);
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

void Patch::selection(){
    std::vector<Individual*> indivs = this->get_all_individuals();
    double w;
    bool surv;

    double max_fitness = 0.0;
    for (Individual* indiv : indivs){
        w = indiv->get_fitness();
        if (w > max_fitness){
            max_fitness = w;
        }
    }


    int K = this->K;
    int n = indivs.size();

    for (Individual* indiv : indivs){
        double fitness = double(indiv->get_fitness())/double(max_fitness);
        double k_prime = double(K * fitness);


        double prob = this->beverton_holt_prob(n, k_prime);
        surv = false;
        if (real_uniform(0,1, main_generator) < prob){
            surv =  true;
        }

        if (!surv){
            remove_individual(indiv);
            delete indiv;
        }
    }
}

double Patch::beverton_holt_prob(int n, double k_prime){
    int b = params["AVG_NUM_OFFSPRING_PER_FEMALE"];

    double prop_full = double(n)/double(k_prime);
    double prob = double(1.0) / double(1 + (double(b)/double(2) - 1)*(prop_full));
    return prob;
}

std::vector<std::vector<Individual*>> Patch::split_by_sex(){
    std::vector<Individual*> males;
    std::vector<Individual*> females;

    for (Individual* indiv: this->get_all_individuals()){
        if (indiv->get_sex()== 1){
           males.push_back(indiv);
        }
        else if (indiv->get_sex() == 0){
           females.push_back(indiv);
        }
    }
    std::vector<std::vector<Individual*>> res;
    res.push_back(females);
    res.push_back(males);
    return res;
}

void Patch::add_to_next_gen(Individual* indiv){
    this->next_gen->push_back(indiv);
}

void Patch::replace_current_gen(){
    delete this->individuals;
    this->individuals = this->next_gen;
    this->next_gen = new std::vector<Individual*>;
}
