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
        dist_decay_space = [1.0]
        n_ef_space = [1, 2, 3, 4]
        n_patches_space = [50]
        ef_space = [0.2, 0.5, 0.8]
        lw_space = [0.05, 0.2]

        n_rep = 30

        n_threads = n_rep * len(n_patches_space) * len(total_indiv_space) * len(dist_decay_space) * len(area_decay_space)

        print('Num of Procs: %d' % n_threads)
        #pool = ThreadPool(n_threads)

        lpef = [0, 20, 10, 7, 5]

        treatment_ct = 0

        for dist_decay in dist_decay_space:
            for n_patches in n_patches_space:
                for n_ef in num_ef_space:
                    for ef_h in ef_space:
                        for lw in lw_space:
                            treatment_ct += 1
                            path = 'treatment%d_rep' % (treatment_ct)

                            for rep in range(n_rep):
                                os.chdir(this_dir)
                                params = {}
                                for param in param_table:
                                    params[param] = param_table[param]['default']


                                params["BASE_MIGRATION_RATE"] = 0.05

                                params["MEAN_LOCUS_WEIGHT"] = lw
                                params["PATCH_DECAY"] = 1.0
                                params["N_INDIVIDUALS"] = 8000
                                params["INCIDENCE_FUNCTION_DECAY"] = dist_decay
                                params["NUM_PATCHES"] = n_patches
                                params["ENV_FACTOR_H_VALUE"] = ef_h
                                params["NUM_ENV_FACTORS"] = n_ef
                                params["NUM_LOCI_PER_EF"] = lpef[n_ef]

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
