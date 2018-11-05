
#include "MigrationTracker.h"
#include "Individual.h"
#include "Patch.h"

MigrationTracker::MigrationTracker(){

	int n_patches = patches->size();
	this->attempted_migration_matrix = new double*[n_patches];
	this->successful_migration_matrix = new double*[n_patches];
	for (int i = 0; i < n_patches; i++){
		this->attempted_migration_matrix[i] = new double[n_patches];
		this->successful_migration_matrix[i] = new double[n_patches];
	}

	for (int i = 0; i < n_patches; i++){
		for (int j = 0; j < n_patches; j++){
			this->attempted_migration_matrix[i][j] = 0;
			this->successful_migration_matrix[i][j] = 0;
		}
	}

}


void MigrationTracker::note_attempted_migration(Patch* from, Patch* to){
	int x = from->get_id();
	int y = to->get_id();

	this->attempted_migration_matrix[x][y]++;
}


void MigrationTracker::note_successful_migration(Patch* from, Patch* to){
	int x = from->get_id();
	int y = to->get_id();

	this->successful_migration_matrix[x][y]++;
}


double MigrationTracker::get_emigration(Patch* from, Patch* to){
	int x = from->get_id();
	int y = to->get_id();

	return double(attempted_migration_matrix[x][y]) / double(from->get_size());
}


double MigrationTracker::get_immigration(Patch* from, Patch* to){
	int x = from->get_id();
	int y = to->get_id();

	return double(attempted_migration_matrix[x][y]) / double(to->get_size());
}

double MigrationTracker::get_successful_migration(Patch* from, Patch* to){
	int x = from->get_id();
	int y = to->get_id();

	return double(successful_migration_matrix[x][y]) / double(to->get_size());
}

double MigrationTracker::get_eff_migration(Patch* patch_i){
	int eff_migrants = 0;

	for (Individual* indiv : patch_i->get_all_individuals()){
		if (indiv->parent_was_migrant){
			eff_migrants++;
		}
	}


	return double(eff_migrants) / double(patch_i->get_size());
}

int MigrationTracker::get_num_indiv(Patch* from, Patch* to){
	int x = from->get_id();
	int y = to->get_id();
	return attempted_migration_matrix[x][y];
}

void MigrationTracker::reset_migration_matrix(){
	int n_patches = patches->size();

	for (int i = 0; i < n_patches; i++){
		for (int j = 0; j < n_patches; j++){
			this->attempted_migration_matrix[i][j] = 0;
			this->successful_migration_matrix[i][j] = 0;
		}
	}
}
