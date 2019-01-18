#include "include.h"
#include "Patch.h"
#include "Individual.h"
#include "EnvFactor.h"
#include "MigrationTracker.h"

void log_patch(Patch* patch_i){
    std::string patch_file = "patches.csv";
    std::ofstream file;
    file.open(patch_file.c_str(), std::ios::app);

    int n_ef = params["NUM_ENV_FACTORS"];


    if(is_file_empty(patch_file)){
        file << "patch_num,generation,x,y,n_indiv,K,prop_of_k,mean_w,se_w,";
        for (int i = 0; i < n_ef; i++){
            file << "ef" << std::to_string(i) << ",";
        }
        file << "\n";
    }

    double w_bar = mean_w(patch_i);
    double w_sd = sd_w(patch_i, w_bar);

    int n = patch_i->get_size();
    int k = patch_i->get_K();

    file << patch_i->get_id() << "," << generation << "," << patch_i->get_x() << "," << patch_i->get_y() << "," << n << "," << k << "," << double(n)/double(k) << "," << w_bar  << "," << w_sd;

    for (double val : patch_i->get_env_factors()){
        file << "," << std::to_string(val);
    }
    file << "\n";
}


double sd_w(Patch* patch_i, double w_bar){
    double se = 0;
    int ct = 0;
    for (Individual* indiv: patch_i->get_all_individuals()){
        se += pow((indiv->get_fitness()-w_bar),2);
        ct++;
    }

    return sqrt(double(se)/double(ct));
}

double mean_w(Patch* patch_i){
    double s = 0;
    int ct = 0;
    for (Individual* indiv: patch_i->get_all_individuals()){
        s += indiv->get_fitness();
        ct++;
    }

    if (ct < 1){
        return 0.0;
    }

    return double(s)/double(ct);
}

void log_fst(int locus, double fst, std::string type){
    std::string fst_file = "f_st.csv";
    std::ofstream file;
    file.open(fst_file.c_str(), std::ios::app);

    if(is_file_empty(fst_file)){
        file << "generation,locus,type,F_st\n";
    }
    file << generation << "," << locus << "," << type << "," << fst << "\n";
}


void log_attempted_migration(Patch* from, Patch* to){
    std::string migration_file = "attempted_migration.csv";
    std::ofstream file;
    file.open(migration_file.c_str(), std::ios::app);
    if(is_file_empty(migration_file)){
        file << "patch_from_num,patch_to_num,generation,num_indiv,prop_of_old_patch_migrants,prop_of_new_patch_migrants\n";
    }

    double em = migration_tracker->get_emigration(from, to);
    double im = migration_tracker->get_immigration(from, to);
    int n_indiv = migration_tracker->get_num_indiv(from, to);
    if (n_indiv > 0){
        file << from->get_id() << "," << to->get_id() << "," << generation << "," << n_indiv << "," << em << "," << im << "\n";
    }
}


void log_successful_migration(Patch* from, Patch* to){
    std::string migration_file = "successful_migration.csv";
    std::ofstream file;
    file.open(migration_file.c_str(), std::ios::app);
    if(is_file_empty(migration_file)){
        file << "patch_from_num,patch_to_num,generation,prop_of_new_patch_migrants\n";
    }

    double im = migration_tracker->get_successful_migration(from, to);
    if (im > 0.0){
        file << from->get_id() << "," << to->get_id() << "," << generation << "," << im <<  "\n";
    }
}

void log_eff_migration(Patch* patch_i){
    std::string eff_migration_file = "eff_migration.csv";
    std::ofstream file;
    file.open(eff_migration_file.c_str(), std::ios::app);
    if(is_file_empty(eff_migration_file)){
        file << "patch_num,generation,eff_migration\n";
    }

    file << patch_i->get_id() << "," << generation << "," << double(migration_tracker->get_eff_migration(patch_i)) << "\n";
 }


void log_population(Patch* patch_i){
    std::string population_file = "population.csv";
    std::ofstream file;
    file.open(population_file.c_str(), std::ios::app);

    if(is_file_empty(population_file)){
        file << "patch_num,generation,num_indiv,proportion_of_k\n";
    }

    int n = patch_i->get_size();

    file << patch_i->get_id() << "," << generation << ","  << n << "," << double(n)/double(patch_i->get_K()) << "\n";

}

void log_linkage(int patch_num, int l1, int l2, double D, std::string type){
    std::string linkage_file = "linkage.csv";

    std::ofstream file;
    file.open(linkage_file.c_str(), std::ios::app);
    if(is_file_empty(linkage_file)){
        file << " patch_num,generation,locus1,locus2,D,type\n";
    }

    file  << patch_num << "," << generation << "," << l1 << "," << l2 << "," << D << "," << type <<  "\n";
}

void log_global_linkage(int l1, int l2, double D, std::string type){
    std::string linkage_file = "global_linkage.csv";

    std::ofstream file;
    file.open(linkage_file.c_str(), std::ios::app);
    if(is_file_empty(linkage_file)){
        file << "generation,locus1,locus2,D,type\n";
    }

    file << generation << "," << l1 << "," << l2 << "," << D << "," << type <<  "\n";
}

void log_pairwise_linkage(int patch_i, int patch_j, double D, std::string type){
    std::string linkage_file = "pairwise_linkage.csv";

    std::ofstream file;
    file.open(linkage_file.c_str(), std::ios::app);
    if(is_file_empty(linkage_file)){
        file << "generation,patch_i,patch_j,D,type\n";
    }

    file << generation << "," << patch_i << "," << patch_j << "," << D << "," << type <<  "\n";
}


void log_allele_freq(int patch_num, int locus, double allele_val, double freq, std::string type){
    std::string allele_freq_file = "allele_freq.csv";
    std::ofstream file;
    file.open(allele_freq_file.c_str(), std::ios::app);
    if(is_file_empty(allele_freq_file)){
        file << "patch_num,generation,locus,allele_val,frequency,type\n";
    }

    file << patch_num << "," << generation << "," << locus << "," << allele_val << "," << freq << "," << type << "\n";
}

void log_locus(int l, int ef, double strength, std::string type){
    std::string locus_file = "loci.csv";
    std::ofstream file;

    file.open(locus_file.c_str(), std::ios::app);
    if(is_file_empty(locus_file)){
        file << "locus,ef,selection_strength,type\n";
    }

    file << l << "," << ef << "," << strength <<  "," << type << "\n";
}

void log_colonization(int patch_num){
    std::string col_file = "colonization.csv";
    std::ofstream file;

    file.open(col_file.c_str(), std::ios::app);
    if(is_file_empty(col_file)){
        file << "generation,patch_num\n";
    }
    file << generation << "," << patch_num << "\n";
}
void log_extinction(int patch_num){
    std::string ext_file = "extinction.csv";
    std::ofstream file;

    file.open(ext_file.c_str(), std::ios::app);
    if(is_file_empty(ext_file)){
        file << "generation,patch_num\n";
    }

    file << generation << "," << patch_num << "\n";
}

void log_env_factors(int x, int y){
    int n_ef = params["NUM_ENV_FACTORS"];

    std::string env_factor_file = "env_factors.csv";
    std::ofstream file;
    file.open(env_factor_file.c_str(), std::ios::app);




    if(is_file_empty(env_factor_file)){
        file << "generation,x,y";

        for (int i = 0; i < n_ef; i++){
            file << "," << "ef" << i;
        }
        file << "\n";
    }

    file << generation << "," << x << "," << y;

    for (EnvFactor* ef : *envFactors){
        file << "," << ef->get_cell_value(x,y);
    }

    file << "\n";
}

bool is_file_empty(std::string path){
    std::ifstream pFile;
    pFile.open(path.c_str());
    return pFile.peek() == std::ifstream::traits_type::eof();
}
