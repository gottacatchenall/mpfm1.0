
#ifndef MIGRATION_TRACKER_H
#define MIGRATION_TRACKER_H

#include "include.h"

class MigrationTracker{
	private:
		std::vector<std::vector<double>> attempted_migration_matrix;
		std::vector<std::vector<double>> successful_migration_matrix;
		std::vector<std::vector<double>> dispersal_kernal;
	public:
		MigrationTracker();
		void note_attempted_migration(Patch* from, Patch* to);
		void note_successful_migration(Patch* from, Patch* to);
		int get_num_indiv(Patch* from, Patch* to);
		double get_eff_migration(Patch* patch_i);
		double get_emigration(Patch* from, Patch* to);
		double get_immigration(Patch* from, Patch* to);
		double get_successful_migration(Patch* from, Patch* to);
		std::vector<double> get_dispersal_row(int num);
		void reset_migration_matrix();
};

#endif
