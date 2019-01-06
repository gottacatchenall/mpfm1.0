#! /usr/bin/env python3

from src.start_run import start_run, print_info, print_run, create_batch_run_file
from src.parse_args import setup_arg_parser, read_param_table, read_batch_file
import os, copy, multiprocessing
import numpy as np

def main():
    this_dir = os.path.abspath('.')

    param_table = read_param_table()
    parser = setup_arg_parser(param_table)
    args = vars(parser.parse_args())
    print_info()

    if (args['BATCH']):
        param_list = read_batch_file(args['BATCH'][0])

        for i,params in enumerate(param_list):
            num_rep = int(params['NUM_REPLICATES'])
            procs = []
            for rep in range(num_rep):
                this_rep_params = copy.deepcopy(params)
                this_rep_params["DATA_DIRECTORY"] = params["DATA_DIRECTORY"] + str(rep)
                this_rep_params["RANDOM_SEED"] = str(np.random.randint(0, 100000000))
                this_rep_params["EF_RANDOM_SEED"] = str(np.random.randint(0, 100000000))
                this_rep_params["PATCH_RAND OM_SEED"] = str(np.random.randint(0, 100000000))
                this_rep_params["GENOME_RANDOM_SEED"] = str(np.random.randint(0, 100000000))

                os.chdir(this_dir)
                #create_batch_run_file(this_dir, this_rep_params)
                #print_run(this_rep_params, i, rep)
                p = multiprocessing.Process(target=start_run, args=(this_rep_params,))
                p.start()
                procs.append(p)
            for p in procs:
                p.join()
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
