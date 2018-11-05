
#ifndef MIGRATION_TRACKER_H
#define MIGRATION_TRACKER_H

#include "include.h"

class MigrationTracker{
	private:
		double** attempted_migration_matrix;
		double** successful_migration_matrix;
	public:
		MigrationTracker();
		void note_attempted_migration(Patch* from, Patch* to);
		void note_successful_migration(Patch* from, Patch* to);
		int get_num_indiv(Patch* from, Patch* to);
		double get_eff_migration(Patch* patch_i);
		double get_emigration(Patch* from, Patch* to);
		double get_immigration(Patch* from, Patch* to);
		double get_successful_migration(Patch* from, Patch* to);
		void reset_migration_matrix();
};

#endif
