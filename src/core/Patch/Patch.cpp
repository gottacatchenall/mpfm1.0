#include "include.h"
#include "Patch.h"
#include "Individual.h"
#include "EnvFactor.h"

int Patch::id_counter = 0;

Patch::Patch(int x, int y, double K){
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
    bool surv;

    int K = this->K;
    int n = indivs.size();

    int b = params["AVG_NUM_OFFSPRING_PER_FEMALE"];
    double max_fitness = 0.0;

    for (Individual* indiv : indivs){
        double this_indiv_fitness = indiv->get_fitness();
        if (this_indiv_fitness > max_fitness){
            max_fitness = this_indiv_fitness;
        }
        double weight = this_indiv_fitness;
        double num_off = (weight * b);
        indiv->set_exp_num_off(num_off);
    }

    // Beverton Holt Method
    for (Individual* indiv : indivs){
        double rel_fitness = double(indiv->get_fitness())/double(max_fitness);
        double k_prime = double(K * rel_fitness);
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
    int n_now = indivs.size();

    if (n > 0 && n_now <= 1){
        log_extinction(this->get_id());
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

void Patch::add_to_migrant_queue(Individual* indiv){
    this->migrant_queue.push_back(indiv);
}

void Patch::add_migrants_to_patch(){
    for (Individual* indiv: this->migrant_queue){
        this->add_individual(indiv);
    }
    this->migrant_queue.clear();
}

std::vector<Individual*> Patch::pick_n_random_indivs(int n){
    int index;
    int size;
    Individual* indiv;
    std::vector<Individual*> indivs;

    int ct = 0;
    while (ct < n){
        size = this->individuals->size();
        index = int_uniform(0,size-1, main_generator);
        indiv = (*this->individuals)[index];
        if (!indiv->has_migrated){
            indiv->has_migrated = true;
            indivs.push_back(indiv);
            ct++;
            this->remove_individual(indiv);
        }
    }

    return indivs;
}


void Patch::add_to_next_gen(Individual* indiv){
    this->next_gen->push_back(indiv);
}

void Patch::replace_current_gen(){
    delete this->individuals;
    this->individuals = this->next_gen;
    this->next_gen = new std::vector<Individual*>;
}
