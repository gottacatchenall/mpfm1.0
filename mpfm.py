#! /usr/bin/env python3

from src.start_run import start_run, print_info, print_run, create_batch_run_file, create_run_directory
from src.parse_args import setup_arg_parser, read_param_table, read_batch_file
import os, copy, multiprocessing, threading
import numpy as np
#from multiprocessing.Dummy import Pool as ThreadPool

def main():
    this_dir = os.path.abspath('.')

    param_table = read_param_table()
    parser = setup_arg_parser(param_table)
    args = vars(parser.parse_args())
    print_info()

    if (args['BATCH']):
        area_decay_space = [1.0, 2.0, 3.0]
        dist_decay_space = [1.0, 2.0, 3.0]
        total_indiv_space = [4000, 6000, 8000]

        n_patches_space = [5, 10, 25]
        n_rep = 30

        n_threads = n_rep * len(n_patches_space) * len(total_indiv_space) * len(dist_decay_space) * len(area_decay_space)

        print('Num of Procs: %d' % n_threads)
        #pool = ThreadPool(n_threads)



        treatment_ct = 0

        for area in area_decay_space:
            for dist_decay in dist_decay_space:
                for n_patches in n_patches_space:
                    for n_indiv in total_indiv_space:

                        treatment_ct += 1
                        path = 'K%d_NP%d_AD%.2f_DE%.2f_rep' % (n_indiv, n_patches, area, dist_decay)

                        for rep in range(n_rep):
                            os.chdir(this_dir)
                            params = {}
                            for param in param_table:
                                params[param] = param_table[param]['default']


                            params["BASE_MIGRATION_RATE"] = 0.05
                            params["MEAN_LOCUS_WEIGHT"] = 0.2

                            params["PATCH_DECAY"] = area
                            params["N_INDIVIDUALS"] = n_indiv
                            params["INCIDENCE_FUNCTION_DECAY"] = dist_decay
                            params["NUM_PATCHES"] = n_patches

                            params["DATA_DIRECTORY"] = path + str(rep)
                            params["RANDOM_SEED"] = np.random.randint(0, 100000000)
                            params["EF_RANDOM_SEED"] = np.random.randint(0, 100000000)
                            params["PATCH_RANDOM_SEED"] = np.random.randint(0, 100000000)
                            params["GENOME_RANDOM_SEED"] = np.random.randint(0, 100000000)


                            # Each treatment goes into unique lb_cmd_file
                            create_batch_run_file(this_dir, params, treatment_ct)

    else:
        params = {}
        for param in param_table:
            if args[param]:
                params[param] = args[param][0]
            else:
                params[param] = param_table[param]['default']
        run = start_run(params)


if __name__ == '__main__':
    main()
