
#include "MigrationTracker.h"
#include "Individual.h"
#include "Patch.h"

MigrationTracker::MigrationTracker(){
	int n_patches = patches->size();

	double DECAY_STR = params["INCIDENCE_FUNCTION_DECAY"];

	for (int i = 0; i < n_patches; i++){
		this->attempted_migration_matrix.push_back(std::vector<double>());
		this->successful_migration_matrix.push_back(std::vector<double>());
		this->dispersal_kernal.push_back(std::vector<double>());

		for (int j = 0; j < n_patches; j++){
			this->attempted_migration_matrix[i].push_back(0.0);
			this->successful_migration_matrix[i].push_back(0.0);
			this->dispersal_kernal[i].push_back(0.0);
		}
	}

	for (int i = 0; i < n_patches; i++){
		int x1 = (*patches)[i]->get_x();
		int y1 = (*patches)[i]->get_y();

		double row[n_patches];
		double row_sum = 0;
		for (int j = 0; j < n_patches; j++){
			if (i != j){
				int x2 = (*patches)[j]->get_x();
				int y2 = (*patches)[j]->get_y();
				double x_dist = x2 - x1;
	            double y_dist = y2 - y1;

	            double dist = sqrt(pow((x_dist),2)+pow((y_dist),2));
				double potential = exp(-1*DECAY_STR*dist);
				row[j] = potential;
				row_sum += potential;
			}
		}

		for (int j = 0; j < n_patches; j++){
			this->dispersal_kernal[i][j] = double(row[j])/double(row_sum);
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

std::vector<double> MigrationTracker::get_dispersal_row(int num){
	return this->dispersal_kernal[num];
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
