#include "include.h"
#include "Patch.h"
#include "Individual.h"


std::map<std::string, double> params;
std::vector<EnvFactor*>* envFactors;
std::vector<Patch*>* patches;
std::vector<std::vector<double>> dist_matrix;
GenomeDict* genome_dict;

int main(int argc, char* argv[]){

    init();

    int num_gen = params["NUM_GENERATIONS"];
    int census_freq = params["CENSUS_FREQ"];

    for (int generation = 0; generation <= num_gen; generation++){
        migration();
        selection();
        //logging();
        //mating();
        if (generation % census_freq == 0){
            //census();
        }
        update_progress_bar(generation);
    }

    return 0;
}

void migration(){
    for (Patch* patch_i: *patches){
        for (Individual* indiv: patch_i->get_all_individuals() ){
            indiv->migrate();
        }
    }
}

void selection(){
    for (Patch* patch_i: *patches){
        for (Individual* indiv: patch_i->get_all_individuals() ){
        //    indiv->calc_fitness();
        }
        //patch_i->selection();
    }
}

void update_progress_bar(int gen){
    int n_gen = params["NUM_GENERATIONS"];

    double progress = double(gen)/double(n_gen);
    int barWidth = 50;

       std::cout << "[";
       int pos = barWidth * progress;
       for (int i = 0; i < barWidth; ++i) {
           if (i <= pos) std::cout << "=";
           else std::cout << " ";
       }
       std::cout << "] [" << gen << " / " << n_gen << "] \r";
       std::cout.flush();

}
